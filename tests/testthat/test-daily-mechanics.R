# The daily clock's mechanisms, each pinned on its own so that a failure names
# the mechanism rather than only moving a reference value: when settings and
# deployments take effect, the three delay lines, the immunity update, the
# competing draw, the seed, the chemoprevention treated phase, refractory
# boosting, vaccines that follow their cohort, and the positivity guards.

first_diff <- function(a, b) which(a != b)[1]
state_cols <- c("S_count", "D_count", "A_count", "U_count", "Tr_count", "Ph_count")
first_state_diff <- function(a, b) {
  which(apply(abs(as.matrix(a[state_cols]) - as.matrix(b[state_cols])), 1, max) > 0)[1]
}
nets_at <- function(p, day) {
  malariasimulation::set_bednets(p, timesteps = day, coverages = 0.8, retention = 5 * 365,
    dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1), rnm = matrix(0.24, 1, 1),
    gamman = 2.64 * 365)
}

test_that("settings the IBM reads on the day change on their own timestep", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  # coverage scheduled at timestep 100 treats day 100's cases, which reach Tr at
  # the end of day 100: row 100 (the start of day 100) holds nobody treated
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 100, coverages = 0.6)
  o <- run_simulation_ode(110, eqm(p, 20))
  expect_equal(o$ft[o$timestep %in% 99:101], c(0, 0.6, 0.6))
  expect_equal(o$Tr_count[o$timestep == 100], 0)
  expect_gt(o$Tr_count[o$timestep == 101], 0)
  # a first-line switch on day 200 (n_ph fixed so both runs share a chain)
  mk <- function(switch) {
    q <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params,
                                                 malariasimulation::AL_params))
    if (switch) {
      q <- malariasimulation::set_clinical_treatment(q, drug = 1, timesteps = c(1, 200),
                                                     coverages = c(0.5, 0))
      q <- malariasimulation::set_clinical_treatment(q, drug = 2, timesteps = c(1, 200),
                                                     coverages = c(0, 0.5))
    } else {
      q <- malariasimulation::set_clinical_treatment(q, drug = 1, timesteps = 1, coverages = 0.5)
    }
    eqm(q, 20)
  }
  expect_equal(first_state_diff(run_simulation_ode(210, mk(TRUE), tuning = list(n_ph = 4)),
                                run_simulation_ode(210, mk(FALSE), tuning = list(n_ph = 4))), 201L)
  # a resistance step on day 300
  base_p <- malariasimulation::set_clinical_treatment(
    malariasimulation::set_drugs(gp(), list(malariasimulation::AL_params)),
    drug = 1, timesteps = 1, coverages = 0.6)
  res_p <- malariasimulation::set_antimalarial_resistance(base_p, drug = 1, timesteps = c(0, 300),
    artemisinin_resistance_proportion = c(0, 0.8), partner_drug_resistance_proportion = c(0, 0),
    slow_parasite_clearance_probability = c(0, 0), early_treatment_failure_probability = c(0.5, 0.5),
    late_clinical_failure_probability = c(0, 0), late_parasitological_failure_probability = c(0, 0),
    reinfection_during_prophylaxis_probability = c(0, 0), slow_parasite_clearance_time = 5)
  expect_equal(first_state_diff(run_simulation_ode(310, eqm(res_p, 20)),
                                run_simulation_ode(310, eqm(base_p, 20))), 301L)
})

test_that("each deployment first moves the output on the day after it, through its own channel", {
  skip_if_not_installed("malariasimulation")
  day <- 100
  scen <- list(
    nets = nets_at(gp_bands(), day),
    tbv = malariasimulation::set_tbv(gp_bands(), timesteps = day, coverages = 0.9, ages = c(18, 19, 20)),
    cc = malariasimulation::set_carrying_capacity(gp_bands(), timesteps = day,
                                                  carrying_capacity_scalers = matrix(0.5, 1, 1)))
  # nets act on the biting rate at once, which moves the force of infection on
  # mosquitoes; TBV on infectivity, which mosquitoes read delay_gam = 12.5 days
  # later; a carrying-capacity change passes through the aquatic stages and
  # ceiling(dem) - 1 = 9 days of incubation (day 112), and the extra infectious
  # bites reach humans, and the EIR column, de = 12 days after that
  want <- c(nets = day + 1, tbv = day + 1 + 12, cc = day + 12 + 12)
  base <- run_simulation_ode(150, eqm(gp_bands(), 20))
  num <- setdiff(names(base), "timestep")
  for (nm in names(scen)) {
    o <- run_simulation_ode(150, eqm(scen[[nm]], 20))
    d <- apply(abs(as.matrix(o[num]) - as.matrix(base[num])), 1, max)
    expect_equal(which(d > 0)[1], want[[nm]], label = nm)
  }
})

test_that("the EIR, infectivity and incubation delay lines are the IBM's", {
  skip_if_not_installed("malariasimulation")
  base <- run_simulation_ode(140, eqm(gp_bands(), 20))
  # the EIR column is the EIR biting humans today, which moves the day after the
  # nets arrive; infections follow once that EIR has been lagged by de = 12 days
  o <- run_simulation_ode(140, eqm(nets_at(gp_bands(), 100), 20))
  expect_equal(first_diff(o$FOIM, base$FOIM), 101L)
  expect_equal(first_diff(o$EIR, base$EIR), 101L + 12L)
  expect_equal(first_diff(o$n_inc_0_36500, base$n_inc_0_36500), 101L + 12L)
  # an MDA after day 100 changes infectivity from day 101: FOIM follows it
  # floor(delay_gam) = 12 days later, and the EIR biting humans once those
  # mosquitoes have incubated and the EIR has been lagged
  mda <- malariasimulation::set_mda(
    malariasimulation::set_drugs(gp_bands(), list(malariasimulation::SP_AQ_params)),
    drug = 1, timesteps = 100, coverages = 0.8, min_ages = 0, max_ages = 100 * 365)
  o <- run_simulation_ode(140, eqm(mda, 20))
  expect_equal(first_diff(o$FOIM, base$FOIM), 113L)
  # infected mosquitoes leave incubation ceiling(dem) - 1 = 9 days later, so the
  # infectious bites they deliver change from day 123 and reach humans de = 12
  # days after that
  expect_equal(first_diff(o$EIR, base$EIR), 135L)
  # delay_gam = 12.5 is read halfway between the two saved days either side
  g12 <- function(p) { p$delay_gam <- 12; p }
  o12 <- run_simulation_ode(140, eqm(g12(mda), 20))
  b12 <- run_simulation_ode(140, eqm(g12(gp_bands()), 20))
  expect_equal((o$FOIM[113] - base$FOIM[113]) / (o12$FOIM[113] - b12$FOIM[113]), 0.5,
               tolerance = 1e-4)
})

test_that("a boosted person skips the day's immunity decay", {
  skip_if_not_installed("malariasimulation")
  p <- eqm(malariasimulation::get_parameters(), 20)
  pr <- build_inputs(p, 20, timesteps = 10)$pars
  sys <- dust2::dust_system_create(get_generator(), pars = pr, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_system_run_to_time(sys, 1)
  IB1 <- dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))$IB
  n <- pr$n_age
  # de = 12 whole days, so day 1 reads the seeded EIR; each person's share is
  # zeta psi mean(psi) / mean(zeta psi), mean(zeta) being the quadrature's
  EPS <- outer(pr$psi, pr$zeta) * sum(pr$eir0) / sum(pr$het_wt * pr$zeta)
  qb <- 1 - exp(-EPS)
  pb <- qb / (qb * pr$ub_eff + 1)
  dec <- 1 - exp(-1 / pr$d_ib)
  N0 <- pr$S0 + pr$D0 + pr$A0 + pr$U0 + pr$Tr0 + apply(pr$Ph0, 1:2, sum) + apply(pr$Phc0, 1:2, sum)
  ain <- rbind(sum(pr$mu_age_z[, 1] * N0) * pr$het_wt / N0[1, ],
               pr$r_age[-n] * N0[-n, , drop = FALSE] / N0[-1, , drop = FALSE])
  I0 <- pr$IB_init
  ageing <- ain * (rbind(0, I0[-n, , drop = FALSE]) - I0)
  expect_equal(IB1, I0 + pb - (1 - pb) * dec * I0 + ageing, tolerance = 1e-12, ignore_attr = TRUE)
  # decay-then-boost (everyone decays, the boosted then gain 1) is measurably different
  expect_gt(max(abs(IB1 - (I0 + pb - dec * I0 + ageing))), 1e-3)
})

test_that("infection and progression are one competing draw a day, for the bitten only", {
  skip_if_not_installed("malariasimulation")
  p <- eqm(malariasimulation::get_parameters(), 20)
  pr <- build_inputs(p, 20, timesteps = 10)$pars
  sys <- dust2::dust_system_create(get_generator(), pars = pr, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_system_run_to_time(sys, 1)
  U1 <- dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))$U
  n <- pr$n_age
  b <- pr$b0 * (pr$b1 + (1 - pr$b1) / (1 + ((pr$IB_init + pr$acq_offset) / pr$ib0)^pr$kb))
  qb <- 1 - exp(-outer(pr$psi, pr$zeta) * sum(pr$eir0) / sum(pr$het_wt * pr$zeta))
  # the bitten face infection (hazard -log(1 - b)) and progression together; the
  # unbitten face progression alone
  hb <- -log(1 - b); hA <- -log(1 - pr$rA); hU <- -log(1 - pr$rU)
  pAU <- qb * (1 - (1 - b) * (1 - pr$rA)) * hA / (hb + hA) + (1 - qb) * pr$rA
  leaveU <- qb * (1 - (1 - b) * (1 - pr$rU)) + (1 - qb) * pr$rU      # infected or recovered
  re <- pr$r_age + pr$mu_age_z[, 1]
  ageU <- rbind(0, pr$r_age[-n] * pr$U0[-n, , drop = FALSE])
  expect_equal(U1, pr$U0 + ageU - re * pr$U0 + pAU * pr$A0 - leaveU * pr$U0,
               tolerance = 1e-12, ignore_attr = TRUE)
  # one draw with the bite-averaged infection probability instead is measurably different
  p_inf <- qb * b; h <- -log(1 - p_inf)
  pAU_avg <- (1 - (1 - p_inf) * (1 - pr$rA)) * hA / (h + hA)
  expect_gt(max(abs(pAU_avg - pAU) / pAU), 1e-3)
})

test_that("row 1 is the seed: the state at the start of day 1", {
  skip_if_not_installed("malariasimulation")
  p <- eqm(malariasimulation::get_parameters(), 20)
  pr <- build_inputs(p, 20, timesteps = 5)$pars
  o <- run_simulation_ode(5, p)
  hp <- p$human_population
  expect_equal(o$S_count[1], sum(pr$S0) * hp, tolerance = 1e-12)
  expect_equal(o$A_count[1], sum(pr$A0) * hp, tolerance = 1e-12)
  expect_gt(abs(o$S_count[2] / (sum(pr$S0) * hp) - 1), 1e-6)   # and row 2 is not
})

test_that("a treated run starts where the IBM starts it", {
  skip_if_not_installed("malariasimulation")
  al <- function(ts) {
    q <- malariasimulation::set_drugs(gp_bands(), list(malariasimulation::AL_params))
    eqm(malariasimulation::set_clinical_treatment(q, drug = 1, timesteps = ts, coverages = 0.4), 20)
  }
  untreated <- build_inputs(eqm(gp_bands(), 20), 20, timesteps = 10)
  from1 <- build_inputs(al(1), 20, timesteps = 10)
  from0 <- build_inputs(al(0), 20, timesteps = 10)
  # the IBM draws its humans from the equilibrium under the coverage at timestep
  # 0 and sizes its mosquitoes under the coverage at timestep 1
  # (set_equilibrium()): treatment from timestep 1 seeds untreated humans among
  # mosquitoes sized for a treated population
  expect_equal(from1$pars$S0, untreated$pars$S0)
  expect_equal(sum(from1$pars$Tr0), 0)
  expect_gt(sum(from0$pars$Tr0), 0)
  # (the infectious stock is fixed by the EIR target; the susceptible stock, and
  # the incubation queue, by the force of infection the mosquitoes were sized at)
  expect_equal(from1$pars$Sm0, from0$pars$Sm0)
  expect_equal(from1$pars$inc0, from0$pars$inc0)
  expect_false(isTRUE(all.equal(from1$pars$Sm0, untreated$pars$Sm0)))
  # with treatment in force from timestep 0 both clocks agree and the seed holds
  o <- run_simulation_ode(365, al(0))
  expect_lt(max(abs(pfpr(o) / pfpr(o)[1] - 1)), 0.01)
})

test_that("a chemoprevention round sends the detectable it clears through Tr first", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters(list(prevalence_rendering_min_ages = 91,
                                              prevalence_rendering_max_ages = 1825))
  p <- malariasimulation::set_drugs(p, list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_smc(p, drug = 1, timesteps = 100, coverages = 0.9,
                                  min_ages = 91, max_ages = 1825)
  p <- eqm(p, 120)
  inp <- build_inputs(p, 120, timesteps = 200)
  o <- run_simulation_ode(200, p)
  expect_equal(pop_total(o), rep(p$human_population, nrow(o)), tolerance = 1e-10)
  expect_true(all(as.matrix(o[state_cols]) >= -1e-9))
  # nobody is in Tr without treatment, until the round: then the dosed clinical
  # and LM-detectable asymptomatic are, and they leave at the treated rate
  expect_equal(o$Tr_count[o$timestep <= 100], rep(0, 100))
  tr <- o$Tr_count[o$timestep %in% 101:110]
  expect_gt(tr[1], 0)
  expect_equal(tr[-1] / tr[-length(tr)], rep(1 - inp$pars$rT, 9), tolerance = 0.02)
  # they stay LM-detectable: prevalence in the dosed band falls over the Tr stay,
  # not overnight
  pr <- pfpr(o, "91_1825")
  expect_gt(pr[102], 0.5 * pr[100])
  expect_lt(pr[110], 0.6 * pr[100])
  # the protection chain after Tr carries the dose's protection less the Tr stay
  expect_gt(inp$meta$n_phct, 1L)
  expect_lt(1 / inp$pars$rP_ct, 1 / inp$pars$rP_c)
})

test_that("refractory boosting follows where the boosted have got to", {
  skip_if_not_installed("malariasimulation")
  p <- eqm(malariasimulation::get_parameters(), 50)
  pr <- build_inputs(p, 50, timesteps = 10)$pars
  sys <- dust2::dust_system_create(get_generator(), pars = pr, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_system_run_to_time(sys, 1)
  ICA1 <- dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))$ICA
  ## day 1's rates, written out independently of the model
  n <- pr$n_age
  qb <- 1 - exp(-outer(pr$psi, pr$zeta) * sum(pr$eir0) / sum(pr$het_wt * pr$zeta))
  b <- pr$b0 * (pr$b1 + (1 - pr$b1) / (1 + ((pr$IB_init + pr$acq_offset) / pr$ib0)^pr$kb))
  ICA0 <- pr$ICA_init
  ICM <- pr$PM * outer(pr$icm_factor, colSums(ICA0 * pr$mask20))
  phi <- pr$phi0 * (pr$phi1 + (1 - pr$phi1) / (1 + (((ICA0 + pr$acq_offset) + ICM) / pr$ic0)^pr$kc))
  hb <- -log(1 - b); hA <- -log(1 - pr$rA); hU <- -log(1 - pr$rU)
  iA <- qb * (1 - (1 - b) * (1 - pr$rA)) * hb / (hb + hA)
  pAU <- qb * (1 - (1 - b) * (1 - pr$rA)) * hA / (hb + hA) + (1 - qb) * pr$rA
  iU <- qb * (1 - (1 - b) * (1 - pr$rU)) * hb / (hb + hU)
  pUS <- qb * (1 - (1 - b) * (1 - pr$rU)) * hU / (hb + hU) + (1 - qb) * pr$rU
  p_inf <- qb * b
  inf <- p_inf * pr$S0 + iA * pr$A0 + iU * pr$U0
  mu <- pr$mu_age_z[, 1]
  ## where someone boosted today is on each of the next uc_eff days (no treatment
  ## here, so every clinical case is in D), and the refractory days that makes
  aA <- 1 - phi; aD <- phi; aU <- 0 * phi; KA <- aA; KU <- aU
  for (d in seq_len(pr$wc - 1)) {
    nA <- aA * (1 - pAU - phi * iA - mu) + aD * pr$rD + (1 - phi) * aU * iU
    nD <- aD * (1 - pr$rD - mu) + phi * (aA * iA + aU * iU)
    nU <- aU * (1 - iU - pUS - mu) + aA * pAU
    aA <- nA; aD <- nD; aU <- nU; KA <- KA + aA; KU <- KU + aU
  }
  N0 <- pr$S0 + pr$D0 + pr$A0 + pr$U0 + pr$Tr0 + apply(pr$Ph0, 1:2, sum) + apply(pr$Phc0, 1:2, sum)
  pbC <- inf / (1 + KA * iA + KU * iU) / N0
  dec <- 1 - exp(-1 / pr$d_ica)
  ain <- rbind(sum(mu * N0) * pr$het_wt / N0[1, ],
               pr$r_age[-n] * N0[-n, , drop = FALSE] / N0[-1, , drop = FALSE])
  ageing <- ain * (rbind(0, ICA0[-n, , drop = FALSE]) - ICA0)
  expect_equal(ICA1, ICA0 + pbC - (1 - pbC) * dec * ICA0 + ageing, tolerance = 1e-12,
               ignore_attr = TRUE)
  ## the renewal form, the same p / (p u + 1) in every state, under-boosts
  w <- pr$uc_eff
  pb_renewal <- (p_inf * pr$S0 / (p_inf * w + 1) + iA * pr$A0 / (iA * w + 1) +
                   iU * pr$U0 / (iU * w + 1)) / N0
  expect_gt(sum(pbC * N0), sum(pb_renewal * N0))
  expect_gt(max((pbC - pb_renewal) / pbC), 0.05)
})

test_that("a mass PEV campaign protects the cohort it vaccinated as it ages", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_mass_pev(malariasimulation::get_parameters(),
    profile = malariasimulation::rtss_profile, timesteps = 100, coverages = 0.8,
    min_ages = 150, max_ages = 510, min_wait = 0, booster_spacing = 365,
    booster_coverage = matrix(0), booster_profile = list(malariasimulation::rtss_booster_profile))
  inp <- build_inputs(eqm(p, 20), 20, timesteps = 1200)
  pv <- inp$pars$pev_vals; tt <- inp$pars$pev_times
  lo <- inp$meta$age_lo; hi <- inp$meta$age_hi
  # 360 days in, the people vaccinated on day 100, aged 150 to 510 days then, are
  # 260 days older: nobody younger than 410 days is protected (the band it was
  # vaccinated in has been refilled from births), and the cohort itself still is
  t1 <- 360
  prot <- 1 - pv[, tt == t1]
  expect_true(all(prot[hi <= 150 + 260 - 1] < 1e-12))
  expect_true(all(prot[lo >= 150 + 260 & hi <= 510 + 260] > 0))
  # and before its efficacy started, nobody was
  expect_true(all(pv[, tt == 100 + max(p$pev_doses) - 1] == 1))
})

test_that("the positivity guards leave room for the day's ageing and deaths", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                    list(malariasimulation::AL_params))
  p <- eqm(malariasimulation::set_mda(p, drug = 1, timesteps = 100, coverages = 0.8,
                                      min_ages = 0, max_ages = 100 * 365), 20)
  fine <- default_age_lower(n_group = 209)
  inp <- build_inputs(p, 20, age_lower = fine, timesteps = 200)
  room <- 1 - max(inp$pars$r_age + apply(inp$pars$mu_age_z, 1, max))
  expect_lte(inp$meta$n_phc * inp$pars$rP_c, room)
  expect_error(build_inputs(p, 20, age_lower = fine, n_phc = inp$meta$n_phc + 1, timesteps = 200),
               "more than it holds")
  o <- run_simulation_ode(200, p, tuning = list(age_lower = fine))
  expect_true(all(as.matrix(o[state_cols]) >= -1e-12))
  # a chain with no drug behind it is not checked
  expect_silent(build_inputs(eqm(malariasimulation::get_parameters(), 20), 20,
                             n_ph = 40, timesteps = 10))
})

test_that("tuning is validated however it arrives", {
  expect_warning(tn <- as_ode_tuning(list(n_sub = 4, rtol = 1e-6)), "ignored")
  expect_equal(tn$n_sub, 4L)
  expect_error(ode_tuning(default_age_lower(), NULL, NULL, 32L, NULL, 1e-6), "named")
  expect_error(ode_tuning(n_sub = Inf), "whole number")
  stale <- structure(list(age_lower = default_age_lower(), n_ph = NULL, n_phc = NULL,
                          n_sub = 0L, odin_file = NULL), class = "fleet_ode_tuning")
  expect_error(as_ode_tuning(stale), "n_sub")
  skip_if_not_installed("malariasimulation")
  p <- eqm(malariasimulation::get_parameters(), 20)
  expect_error(run_simulation_ode(10, p, list(n_sub = 4, rtol = 1e-6)), "fourth")
})
