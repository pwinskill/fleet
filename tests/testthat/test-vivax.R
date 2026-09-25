# P. vivax.
#
# fleet replicates malariasimulation's parameter translation rather than calling
# it, so that this package is self-contained and version-robust (same reasoning
# as .back_translations for falciparum). A replicated table is only as good as
# the test that it still agrees with the original, so that is the first thing
# here -- and it is a real risk, not a theoretical one: the vivax table has one
# entry that is not a plain alias (phi_D_min = phi0 * phi1, because a
# malariasimulation vivax list stores phi1 as a ratio), and an alias that
# silently became a rename would change the equilibrium without erroring.
# The daily mechanics of the vivax block are pinned one by one in
# test-vivax-mechanics.R, and its output columns against the IBM's in
# test-output-parity.R.

test_that("the replicated vivax translation agrees with malariasimulation's", {
  skip_if_not_installed("malariasimulation")
  ms_tr <- get("translate_vivax_parameters", asNamespace("malariasimulation"))
  p <- malariasimulation::get_parameters(parasite = "vivax")
  a <- ms_tr(p)
  b <- translate_vivax_parameters(p)
  expect_setequal(names(a), names(b))
  for (nm in names(a)) {
    expect_equal(b[[nm]], a[[nm]], info = paste("vivax translation differs for", nm))
  }
})

test_that("the vivax translation is additive, and converts phi1 to a level", {
  skip_if_not_installed("malariasimulation")
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

test_that("a vivax list reaches the shared parameters under the falciparum names", {
  skip_if_not_installed("malariasimulation")
  # The model shares one set of progression rates, decays, windows and the
  # clinical curve between the parasites: a vivax list carries da, dd, dt, rc,
  # uc, rm, pcm and phi0/phi1/ic0/kc under the names falciparum uses, so the
  # vivax block reads the vivax values through them. Pin that it does.
  p <- eqm(malariasimulation::get_parameters(parasite = "vivax"), 3)
  x <- build_inputs(p, 3, timesteps = 30)$pars
  expect_equal(x$rA, 1 - exp(-1 / p$da))
  expect_equal(x$rD, 1 - exp(-1 / p$dd))
  expect_equal(x$rT, 1 - exp(-1 / p$dt))
  expect_equal(c(x$d_ica, x$uc_eff, x$PM), c(p$rc, ceiling(p$uc) - 1, p$pcm))
  expect_equal(c(x$phi0, x$phi1, x$ic0, x$kc), c(p$phi0, p$phi1, p$ic0, p$kc))
  expect_equal(c(x$cD, x$cU), c(p$cd, p$cu))
  expect_equal(c(x$d_iaa, x$ua_eff, x$cA_v, x$bv, x$ff, x$gammal),
               c(p$ra, ceiling(p$ua) - 1, p$ca, p$b, p$f, p$gammal))
  # the vivax block at full size, the falciparum one empty
  expect_equal(c(x$pf_on, x$pv_on, x$n_age_v, x$n_het_v, x$n_bat),
               c(0, 1, x$n_age, x$n_het, p$kmax + 1))
  expect_true(all(c(x$S0, x$D0, x$A0, x$U0, x$Tr0) == 0))
  # and the reverse under falciparum: one inert vivax cell
  y <- build_inputs(eqm(malariasimulation::get_parameters(), 3), 3, timesteps = 30)$pars
  expect_equal(c(y$pf_on, y$pv_on, y$n_age_v, y$n_het_v, y$n_hyp), c(1, 0, 1, 1, 1))
  expect_true(all(c(y$Sv0, y$JAv0) == 0))
})

test_that("the vivax equilibrium solves on fleet's own age grid", {
  skip_if_not_installed("malariasimulation")
  # This is the seed fleet's vivax block is built on, so it has to work at
  # fleet's resolution, not only at malariasimulation's 1000-point grid. Ages
  # are in YEARS here -- passing days returns NaN or a singular matrix rather
  # than erroring cleanly, which is a cheap mistake to make and an expensive
  # one to diagnose.
  p <- translate_vivax_parameters(malariasimulation::get_parameters(parasite = "vivax"))
  foim <- vapply(c(0.1, 1, 5, 20, 120), function(E) {
    e <- malariaEquilibriumVivax::vivax_equilibrium(EIR = E, ft = 0.2, p = p,
                                                    age = default_age_lower())
    states <- e$states[c("S", "U", "A", "D", "T", "P")]
    expect_false(any(vapply(states, function(x) any(is.na(x)), logical(1))))
    expect_equal(sum(vapply(states, sum, numeric(1))), 1, tolerance = 1e-8)
    expect_equal(dim(e$states$S)[3], p$kmax + 1)
    e$FOIM
  }, numeric(1))
  expect_true(all(is.finite(foim) & foim > 0))
  expect_false(is.unsorted(foim))     # FOIM rises with transmission
})

test_that("vivax shares its heterogeneity nodes with fleet and the IBM", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("statmod")
  # fleet seeds each heterogeneity stratum separately for falciparum, but the
  # vivax equilibrium returns the stratum dimension already built. That is only
  # usable if the two discretisations are the same one.
  p <- malariasimulation::get_parameters(parasite = "vivax")
  s2 <- p$sigma_squared
  n <- p$n_heterogeneity_groups
  gq <- malariaEquilibrium::gq_normal(n)
  q <- statmod::gauss.quad.prob(n, dist = "normal", mu = -0.5 * s2, sigma = sqrt(s2))
  expect_equal(exp(gq$nodes * sqrt(s2) - s2 / 2), exp(q$nodes))
  expect_equal(gq$weights, q$weights)
})

## ---- the running model -------------------------------------------------------

# A vivax parameter list, optionally with a vivax drug at 50% coverage.
pv_list <- function(drug = NULL, eir = 10, ...) {
  p <- malariasimulation::get_parameters(list(human_population = 10000, ...),
                                         parasite = "vivax")
  if (!is.null(drug)) {
    p <- malariasimulation::set_drugs(p, list(drug))
    p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1,
                                                   coverages = 0.5)
  }
  eqm(p, eir)
}

test_that("the vivax seed is the IBM's: nobody drug-protected, renormalised by stratum", {
  skip_if_not_installed("malariasimulation")
  # create_variables() draws each person's state and batch count from the
  # equilibrium's S, D, A, U and T within their age and heterogeneity group,
  # leaving out its prophylaxis state. With treatment on from timestep 0 the
  # equilibrium holds people in P, and the seed must hold none of them in the
  # prophylaxis chain, and put every age x heterogeneity cell at its share.
  p <- malariasimulation::get_parameters(list(human_population = 10000), parasite = "vivax")
  p <- malariasimulation::set_drugs(p, list(malariasimulation::CQ_params_vivax))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 0.5)
  inp <- build_inputs(eqm(p, 10), 10, timesteps = 30)
  x <- inp$pars
  expect_true(all(x$Phv0 == 0))
  N <- x$Sv0 + x$Dv0 + x$Av0 + x$Uv0 + x$Trv0
  expect_equal(apply(N, c(1, 2), sum), outer(inp$meta$prop, x$het_wt), tolerance = 1e-12)
  expect_gt(sum(x$Trv0), 0)                     # treated at timestep 0: some in Tr
  # the second moments start with no spread: K = J^2 / N in every occupied cell
  occ <- N > 1e-12
  expect_equal(x$KCv0[occ], x$JCv0[occ]^2 / N[occ], tolerance = 1e-10)
  expect_equal(x$KAv0[occ], x$JAv0[occ]^2 / N[occ], tolerance = 1e-10)
})

test_that("a vivax run conserves population and relaxes gently off the IBM's seed", {
  skip_if_not_installed("malariasimulation")
  out <- run_simulation_ode(365, pv_list(eir = 20))
  expect_equal(out$EIR[1], 20, tolerance = 1e-6)
  expect_lt(max(abs(pop_total(out) - 10000)) / 10000, 1e-9)
  # The seed is the IBM's own (malariaEquilibriumVivax on malariasimulation's
  # durations), which neither model holds exactly: both resolve each day in one
  # competing-hazard draw and build up a within-cell spread of immunity the
  # equilibrium does not have. Measured against the IBM they relax together, so
  # the test is only that the relaxation is gentle and bounded.
  pcr <- out$n_detect_pcr_730_3650 / out$n_age_730_3650
  expect_lt(max(abs(pcr / pcr[1] - 1)), 0.02)
})

test_that("with no spread assumed, every node is the cell mean", {
  skip_if_not_installed("malariasimulation")
  # immunity_spread = FALSE must reduce the closure to the plain mean-field
  # curves -- one node, at the mean -- rather than to something else
  p <- pv_list()
  x_on <- build_inputs(p, init_EIR = 10, timesteps = 30)$pars
  p$immunity_spread <- FALSE
  x_off <- build_inputs(p, init_EIR = 10, timesteps = 30)$pars
  expect_equal(x_on$n_q, length(VIVAX_GH_Z))
  expect_equal(sum(x_on$qw), 1)
  expect_equal(x_off[c("spread_on", "n_q", "qz", "qw")],
               list(spread_on = 0, n_q = 1L, qz = 0, qw = 1))
  p$immunity_spread <- "yes"
  expect_error(build_inputs(p, init_EIR = 10), "TRUE or FALSE")
  # and the spread builds up, as it does in the IBM: clinical incidence in
  # school-age children rises off the seed only when the spread is modelled
  bands <- list(clinical_incidence_rendering_min_ages = 1825,
                clinical_incidence_rendering_max_ages = 3650)
  on <- run_simulation_ode(730, do.call(pv_list, c(list(eir = 20), bands)))
  pq <- do.call(pv_list, c(list(eir = 20), bands)); pq$immunity_spread <- FALSE
  off <- run_simulation_ode(730, pq)
  last <- nrow(on)
  expect_gt(on$n_inc_clinical_1825_3650[last], 1.1 * off$n_inc_clinical_1825_3650[last])
})

test_that("radical cure clears batches and fills liver-stage protection", {
  skip_if_not_installed("malariasimulation")
  cq <- pv_list(malariasimulation::CQ_params_vivax)
  tq <- pv_list(malariasimulation::CQ_TQ_params_vivax)
  # the protected levels exist only when a radical-cure drug is on the schedule,
  # a chain matched to the protected days on the daily clock (tafenoquine's
  # Weibull(5, 30): mean 27.0 days, variance 39.9 -> 11 stages)
  x_cq <- build_inputs(cq, init_EIR = 10, timesteps = 365)$pars
  x_tq <- build_inputs(tq, init_EIR = 10, timesteps = 365)$pars
  pd <- .protection_days(5, 30)
  expect_equal(x_cq$n_hyp, x_cq$n_bat)
  expect_equal(x_cq$rc_on, 0)
  expect_equal(x_tq$rc_on, 1)
  expect_equal(x_tq$n_hyp, x_tq$n_bat + erlang_stages(pd$mean, pd$var))
  expect_equal(x_tq$n_hyp - x_tq$n_bat, 11)
  expect_equal(x_tq$mean_ls_vals, pd$mean)
  expect_equal(x_tq$hyp_vals, 0.713)
  expect_equal(x_tq$eff_hyp_vals, 1 * 0.713)
  o_cq <- run_simulation_ode(365, cq)
  o_tq <- run_simulation_ode(365, tq)
  for (o in list(o_cq, o_tq)) expect_lt(max(abs(pop_total(o) - 10000)) / 10000, 1e-9)
  # The same carriers at the seed: malariaEquilibriumVivax has no radical cure,
  # so its batch distribution does not depend on the drug (the IBM starts from
  # it too). Then fewer carriers, relapses and infections under tafenoquine.
  expect_equal(o_tq$n_with_hypnozoites[1], o_cq$n_with_hypnozoites[1])
  last <- nrow(o_cq)
  expect_lt(o_tq$n_with_hypnozoites[last], 0.95 * o_cq$n_with_hypnozoites[last])
  expect_lt(o_tq$n_relapses[last], o_cq$n_relapses[last])
  expect_lt(o_tq$n_detect_pcr_730_3650[last], o_cq$n_detect_pcr_730_3650[last])
})

test_that("vivax has no severe disease, and its own columns hold together", {
  skip_if_not_installed("malariasimulation")
  bands <- list(clinical_incidence_rendering_min_ages = c(0, 1825),
                clinical_incidence_rendering_max_ages = c(1824, 36499),
                incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 36499,
                severe_incidence_rendering_min_ages = 0,
                severe_incidence_rendering_max_ages = 36499)
  o <- run_simulation_ode(60, do.call(pv_list, bands))
  # a severe band set on a vivax list renders, as the IBM's does, and is 0
  expect_true(all(o$n_inc_severe_0_36499 == 0))
  # relapses are a share of all infections; carriers a share of the population
  expect_true(all(o$n_relapses > 0 & o$n_relapses < o$n_infections))
  expect_equal(o$n_infections, o$n_inc_0_36499)
  expect_true(all(o$n_with_hypnozoites > 0 & o$n_with_hypnozoites < pop_total(o)))
  # vivax A is LM-detectable by definition: LM prevalence sits inside PCR, and
  # the IBM renders no p_detect_lm_* for vivax
  expect_true(all(o$n_detect_lm_730_3650 < o$n_detect_pcr_730_3650))
  expect_false(any(grepl("^p_detect_lm_", names(o))))
  # the mean batch count is the carriers' batches over everyone
  expect_true(all(o$hypnozoites_mean > 0 & o$hypnozoites_mean < 10))
})

test_that("under a custom demography vivax sizes mosquitoes as set_equilibrium() does", {
  skip_if_not_installed("malariasimulation")
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
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters(parasite = "vivax")
  p$smc <- TRUE
  expect_error(build_inputs(p, init_EIR = 20), "cannot be combined with SMC")
  p <- malariasimulation::get_parameters(); p$parasite <- "ovale"
  expect_error(build_inputs(p, init_EIR = 20), "'falciparum' or 'vivax'")
})
