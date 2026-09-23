# P. vivax support.
#
# fleet replicates malariasimulation's parameter translation rather than calling
# it, so that this package is self-contained and version-robust (same reasoning
# as .back_translations for falciparum). A replicated table is only as good as
# the test that it still agrees with the original, so that is the first thing
# here -- and it is a real risk, not a theoretical one: the vivax table has one
# entry that is not a plain alias (phi_D_min = phi0 * phi1, because a
# malariasimulation vivax list stores phi1 as a ratio), and an alias that
# silently became a rename would change the equilibrium without erroring.

test_that("the replicated vivax translation agrees with malariasimulation's", {
  ms_tr <- get("translate_vivax_parameters", asNamespace("malariasimulation"))
  p <- malariasimulation::get_parameters(parasite = "vivax")

  a <- ms_tr(p)
  b <- translate_vivax_parameters(p)

  expect_setequal(names(a), names(b))
  for (nm in names(a)) {
    expect_equal(b[[nm]], a[[nm]],
                 info = paste("vivax translation differs for", nm))
  }
})

test_that("the vivax translation is additive, and converts phi1 to a level", {
  p <- malariasimulation::get_parameters(parasite = "vivax")
  b <- translate_vivax_parameters(p)

  # additive: every input field survives under its own name
  expect_true(all(names(p) %in% names(b)))
  expect_equal(b$sigma_squared, p$sigma_squared)
  expect_equal(b$n_heterogeneity_groups, p$n_heterogeneity_groups)

  # the one non-alias: a vivax list stores phi1 as phi_D_min/phi_D_max
  expect_equal(b$phi_D_max, p$phi0)
  expect_equal(b$phi_D_min, p$phi0 * p$phi1)
  expect_lt(b$phi_D_min, b$phi_D_max)

  # rates are inverted durations, not durations
  expect_equal(b$r_D, 1 / p$dd)
  expect_equal(b$r_par, 1 / p$ra)
})

test_that("the vivax equilibrium solves on fleet's own age grid", {
  # This is the seed fleet's vivax block will be built on, so it has to work at
  # fleet's resolution, not only at malariasimulation's 1000-point grid. Ages
  # are in YEARS here -- passing days returns NaN or a singular matrix rather
  # than erroring cleanly, which is a cheap mistake to make and an expensive
  # one to diagnose.
  skip_if_not_installed("malariaEquilibriumVivax")
  p <- translate_vivax_parameters(
    malariasimulation::get_parameters(parasite = "vivax"))
  age <- default_age_lower()

  foim <- vapply(c(0.1, 1, 5, 20, 120), function(E) {
    e <- malariaEquilibriumVivax::vivax_equilibrium(EIR = E, ft = 0.2, p = p, age = age)
    states <- e$states[c("S", "U", "A", "D", "T", "P")]
    expect_false(any(vapply(states, function(x) any(is.na(x)), logical(1))))
    # one unit of population, distributed over age x heterogeneity x batch
    expect_equal(sum(vapply(states, sum, numeric(1))), 1, tolerance = 1e-8)
    # the batch dimension is kmax + 1
    expect_equal(dim(e$states$S)[3], p$kmax + 1)
    e$FOIM
  }, numeric(1))

  expect_true(all(is.finite(foim)))
  expect_true(all(foim > 0))
  expect_false(is.unsorted(foim))     # FOIM rises with transmission
})

test_that("vivax shares its heterogeneity nodes with fleet and the IBM", {
  # fleet seeds each heterogeneity stratum separately for falciparum, but the
  # vivax equilibrium returns the stratum dimension already built. That is only
  # usable if the two discretisations are the same one.
  p <- malariasimulation::get_parameters(parasite = "vivax")
  s2 <- p$sigma_squared
  n <- p$n_heterogeneity_groups

  gq <- malariaEquilibrium::gq_normal(n)
  zeta_fleet <- exp(gq$nodes * sqrt(s2) - s2 / 2)

  q <- statmod::gauss.quad.prob(n, dist = "normal", mu = -0.5 * s2, sigma = sqrt(s2))
  zeta_eq <- exp(q$nodes)

  expect_equal(zeta_fleet, zeta_eq)
  expect_equal(gq$weights, q$weights)
})

## ---- the running model -------------------------------------------------------

# A vivax parameter list at EIR 20, optionally with a vivax drug at 50% coverage.
pv_list <- function(drug = NULL, eir = 20, ...) {
  p <- malariasimulation::get_parameters(list(human_population = 10000, ...),
                                         parasite = "vivax")
  if (!is.null(drug)) {
    p <- malariasimulation::set_drugs(p, list(drug))
    p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1,
                                                   coverages = 0.5)
  }
  eqm(p, eir)
}
pop_of <- function(out) with(out, S_count + D_count + A_count + U_count + Tr_count + Ph_count)

test_that("a vivax run conserves population and relaxes gently off the IBM's seed", {
  skip_if_not_installed("malariaEquilibriumVivax")
  out <- run_simulation_ode(365, pv_list())

  expect_equal(out$EIR[1], 20, tolerance = 1e-6)
  pop <- pop_of(out)
  expect_lt(max(abs(pop - 10000)) / 10000, 1e-6)

  # The seed is the IBM's own (malariaEquilibriumVivax on malariasimulation's
  # durations), which neither model holds exactly: both resolve each day in one
  # competing-hazard draw and build up a within-cell spread of immunity the
  # equilibrium does not have. Measured against the IBM they relax together, so
  # the test is only that the relaxation is gentle and bounded.
  pcr <- out$n_detect_pcr_730_3650 / out$n_age_730_3650
  expect_lt(max(abs(pcr / pcr[1] - 1)), 0.02)
})

test_that("with no spread assumed, every node is the cell mean", {
  skip_if_not_installed("malariaEquilibriumVivax")
  # immunity_spread = FALSE must reduce the closure to the plain mean-field
  # curves -- one node, at the mean -- rather than to something else
  p <- pv_list()
  x_on <- build_inputs(p, init_EIR = 20, timesteps = 30)$pars
  p$immunity_spread <- FALSE
  x_off <- build_inputs(p, init_EIR = 20, timesteps = 30)$pars
  expect_equal(x_on$n_q, length(VIVAX_GH_Z))
  expect_equal(sum(x_on$qw), 1)
  expect_equal(x_off[c("spread_on", "n_q", "qz", "qw")],
               list(spread_on = 0, n_q = 1L, qz = 0, qw = 1))
  # the second moments start with no spread: K = J^2 / N in every cell
  N <- x_on$Sv0 + x_on$Dv0 + x_on$Av0 + x_on$Uv0 + x_on$Trv0 + apply(x_on$Phv0, 1:3, sum)
  occ <- N > 1e-12
  expect_equal(x_on$KCv0[occ], x_on$JCv0[occ]^2 / N[occ], tolerance = 1e-10)
  # and the spread builds up, as it does in the IBM: clinical incidence in
  # school-age children rises off the seed only when the spread is modelled
  bands <- list(clinical_incidence_rendering_min_ages = 1825,
                clinical_incidence_rendering_max_ages = 3650)
  on <- run_simulation_ode(730, do.call(pv_list, c(list(drug = NULL), bands)))
  pq <- do.call(pv_list, c(list(drug = NULL), bands)); pq$immunity_spread <- FALSE
  off <- run_simulation_ode(730, pq)
  last <- nrow(on)
  expect_gt(on$n_inc_clinical_1825_3650[last], 1.1 * off$n_inc_clinical_1825_3650[last])
})

test_that("batch clearance alone decays the hypnozoite ladder exactly", {
  skip_if_not_installed("malariaEquilibriumVivax")
  # With no bites (b = 0) and no relapse (f = 0), batches only clear, one at a
  # time at k * gammal, and people die at the age-independent default rate.
  # Summed over the ladder that is linear: K = sum(k * N_k) obeys
  # dK/dt = -(gammal + mu) K exactly, whatever the batch distribution, because
  # births arrive with none. So the solver must reproduce a pure exponential.
  p <- pv_list()
  x <- build_inputs(p, init_EIR = 20, timesteps = 365)$pars
  x$bv <- 0; x$ff <- 0
  mu <- unique(as.vector(x$mu_age_z))
  expect_length(mu, 1)

  sys <- dust2::dust_system_create(malaria_ode, pars = x, n_particles = 1,
    ode_control = dust2::dust_ode_control(atol = 1e-10, rtol = 1e-8))
  dust2::dust_system_set_state_initial(sys)
  idx <- dust2::dust_unpack_index(sys)
  times <- c(0, 100, 365)
  y <- dust2::dust_system_simulate(sys, times)
  dims <- c(x$n_age_v, x$n_het_v, x$n_hyp)
  kk <- c(seq_len(x$n_bat) - 1, rep(0, x$n_hyp - x$n_bat))
  K <- vapply(seq_along(times), function(ti) {
    N <- Reduce(`+`, lapply(c("Sv", "Dv", "Av", "Uv", "Trv", "Trv_slow"),
                            function(nm) array(y[idx[[nm]], ti], dims)))
    N <- N + apply(array(y[idx$Phv, ti], c(dims, x$n_phv)), 1:3, sum)
    sum(apply(N, 3, sum) * kk)
  }, numeric(1))

  expect_gt(K[1], 1)                       # a real ladder to decay
  expect_equal(K, K[1] * exp(-(x$gammal + mu) * times), tolerance = 1e-6)
})

test_that("radical cure clears batches and fills liver-stage protection", {
  skip_if_not_installed("malariaEquilibriumVivax")
  ms <- asNamespace("malariasimulation")
  cq <- pv_list(ms$CQ_params_vivax)
  tq <- pv_list(ms$CQ_TQ_params_vivax)

  # the protected levels exist only when a radical-cure drug is on the schedule
  x_cq <- build_inputs(cq, init_EIR = 20, timesteps = 365)$pars
  x_tq <- build_inputs(tq, init_EIR = 20, timesteps = 365)$pars
  expect_equal(x_cq$n_hyp, x_cq$n_bat)
  expect_equal(x_tq$n_hyp, x_tq$n_bat + VIVAX_MAX_PH)
  expect_equal(x_tq$hyp_vals, 0.713)
  expect_equal(x_tq$eff_hyp_vals, 1 * 0.713)
  # tafenoquine: Weibull(5, 30) liver-stage protection, mean 30 * gamma(1.2)
  expect_equal(x_tq$mean_ls_vals, 30 * gamma(1.2))

  o_cq <- run_simulation_ode(365, cq)
  o_tq <- run_simulation_ode(365, tq)
  for (o in list(o_cq, o_tq)) expect_lt(max(abs(pop_of(o) - 10000)) / 10000, 1e-6)
  # The same carriers at the seed: malariaEquilibriumVivax has no radical cure
  # and lets treated people gain batches, so its batch distribution does not
  # depend on the drug (the IBM starts from it too). Then fewer carriers,
  # relapses and infections under tafenoquine.
  expect_equal(o_tq$n_with_hypnozoites[1], o_cq$n_with_hypnozoites[1])
  last <- nrow(o_cq)
  expect_lt(o_tq$n_with_hypnozoites[last], 0.95 * o_cq$n_with_hypnozoites[last])
  expect_lt(o_tq$n_relapses[last], o_cq$n_relapses[last])
  expect_lt(o_tq$n_detect_pcr_730_3650[last], o_cq$n_detect_pcr_730_3650[last])
})

test_that("vivax output columns are malariasimulation's, and severe is absent", {
  skip_if_not_installed("malariaEquilibriumVivax")
  bands <- list(clinical_incidence_rendering_min_ages = c(0, 1825),
                clinical_incidence_rendering_max_ages = c(1825, 36500),
                severe_incidence_rendering_min_ages = c(0, 1825),
                severe_incidence_rendering_max_ages = c(1825, 36500))
  pv <- do.call(pv_list, c(list(drug = NULL), bands))
  pf <- eqm(malariasimulation::get_parameters(c(list(human_population = 10000), bands)), 20)
  o_pv <- run_simulation_ode(60, pv)
  o_pf <- run_simulation_ode(60, pf)

  # the same columns for either parasite, so runs stay rbind-able
  expect_identical(names(o_pv), names(o_pf))
  # vivax has no severe pathway; falciparum has no relapses or hypnozoites
  sev <- grep("^n_inc_severe_", names(o_pv), value = TRUE)
  expect_true(all(unlist(o_pv[sev]) == 0))
  expect_true(all(o_pf$n_relapses == 0) && all(o_pf$n_with_hypnozoites == 0))
  # relapses are a share of all infections; carriers a share of the population
  expect_true(all(o_pv$n_relapses > 0 & o_pv$n_relapses < o_pv$n_inc_0_36500))
  expect_true(all(o_pv$n_with_hypnozoites > 0 &
                  o_pv$n_with_hypnozoites < o_pv$n_age_0_36500))
  # vivax A is LM-detectable by definition: LM prevalence sits inside PCR
  expect_true(all(o_pv$n_detect_lm_730_3650 < o_pv$n_detect_pcr_730_3650))
})

test_that("vivax output feeds postie, severe bands following the clinical ones", {
  skip_if_not_installed("malariaEquilibriumVivax")
  skip_if_not_installed("postie")
  # site::site_parameters(parasite = "vivax") sets clinical bands and no severe
  # ones (malariasimulation renders no severe for vivax), and postie::get_rates()
  # needs the two families on the same bands -- so the zero severe columns
  # follow the clinical bands.
  pv <- pv_list(clinical_incidence_rendering_min_ages = c(0, 1825, 5475),
                clinical_incidence_rendering_max_ages = c(1825, 5475, 36500))
  o <- run_simulation_ode(60, pv)
  expect_true(all(c("n_inc_severe_0_1825", "n_inc_severe_1825_5475",
                    "n_inc_severe_5475_36500") %in% names(o)))
  r <- suppressWarnings(postie::get_rates(o))
  expect_true(all(r$severe == 0))
  expect_true(all(r$clinical > 0))
})

test_that("under a custom demography vivax sizes mosquitoes as set_equilibrium() does", {
  skip_if_not_installed("malariaEquilibriumVivax")
  # As for falciparum: set_equilibrium() sizes total_M from the equilibrium under
  # the DEFAULT exponential age structure, and the IBM then realises whatever
  # transmission that density supports under the custom mortality. fleet takes
  # the IBM's density and seeds at the EIR its own equilibrium supports for it.
  dr <- c(0.048, 0.007, 0.003, 0.004, 0.008, 0.020, 0.050, 0.120) / 365
  p <- malariasimulation::get_parameters(list(human_population = 10000), parasite = "vivax")
  p <- malariasimulation::set_demography(p,
    agegroups = round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365),
    timesteps = 0, deathrates = matrix(dr, nrow = 1))
  p <- eqm(p, 3)
  inp <- build_inputs(p, 3)
  expect_equal(inp$meta$total_M_ibm, p$total_M, tolerance = 1e-8)
  expect_equal(inp$meta$total_M, p$total_M, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(inp$meta$eir_seed, 3)))
  # the default demography is unaffected
  expect_equal(build_inputs(pv_list(eir = 3), 3)$meta$eir_seed, 3)
})

test_that("vivax refuses the chemoprevention malariasimulation cannot run", {
  p <- malariasimulation::get_parameters(parasite = "vivax")
  p$smc <- TRUE
  expect_error(build_inputs(p, init_EIR = 20), "cannot be combined with SMC")
})
