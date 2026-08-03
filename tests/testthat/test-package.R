# Package-level: warnings, error paths, edge cases, output contract.

test_that("all intervention modules build silently (no unsupported warnings)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  # custom demography and time-varying EPI coverage are supported (no warning)
  expect_silent(build_inputs(malariasimulation::set_demography(
    gp(), agegroups = c(2920, 36500), timesteps = 0,
    deathrates = matrix(c(1e-4, 1e-3), nrow = 1)), 20))
  epi <- malariasimulation::set_pev_epi(gp(), profile = malariasimulation::rtss_profile,
    timesteps = c(1, 1000), coverages = c(0.5, 0.9), min_wait = 0, age = 180,
    booster_spacing = 360, booster_coverage = matrix(0.8, 2, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  expect_silent(build_inputs(epi, 20))
  expect_silent(build_inputs(gp(overrides = list(model_seasonality = TRUE)), 20))
})

test_that("time-varying EPI coverage: protection tracks each cohort's vaccination-date coverage", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  age_mid <- default_age_lower() * 365 + 30
  # coverage steps up from 0.2 to 0.9 at year 5
  epi <- malariasimulation::set_pev_epi(gp(), profile = malariasimulation::rtss_profile,
    timesteps = c(1, 5 * 365), coverages = c(0.2, 0.9), min_wait = 0, age = 180,
    booster_spacing = 360, booster_coverage = matrix(0, 2, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  ps <- pev_series(epi, age_mid, 10 * 365)
  # pick a protected age group and look along calendar time: cohorts vaccinated
  # before the year-5 step-up carry 0.2 coverage, later cohorts carry 0.9, so the
  # FOI multiplier for that age group must take >1 distinct value over time
  vax_complete <- 180 + max(epi$pev_doses)
  i <- which(age_mid > vax_complete + 200)[1]
  expect_gt(length(unique(round(ps$vals[i, ], 6))), 1)     # genuinely time-varying
  # the coverage step must land at THIS cohort's vaccination date (grid - tsince =
  # 5*365 => grid = 5*365 + tsince), not the calendar step-up (grid = 5*365). The
  # largest drop in the multiplier is the 0.2->0.9 step; pin its calendar time.
  tsince <- age_mid[i] - vax_complete
  last_dose <- max(epi$pev_doses)
  step_g <- ps$times[which.min(diff(ps$vals[i, ])) + 1]
  # ms samples EPI coverage at the FIRST-dose date, so the 0.2->0.9 step lands when this
  # cohort's first dose (last_dose before efficacy onset) crosses year 5:
  # grid - tsince - last_dose = 5*365  =>  grid = 5*365 + tsince + last_dose.
  expect_gt(step_g, 5 * 365 + tsince + last_dose - 45)   # after the cohort's first-dose date
  expect_lt(step_g, 5 * 365 + tsince + last_dose + 45)   # (a current-time bug would step at 5*365)
  # and the deepest protection (highest coverage) must beat a constant-0.2 program
  epi_lo <- malariasimulation::set_pev_epi(gp(), profile = malariasimulation::rtss_profile,
    timesteps = 1, coverages = 0.2, min_wait = 0, age = 180, booster_spacing = 360,
    booster_coverage = matrix(0, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  lo <- pev_series(epi_lo, age_mid, 10 * 365)
  expect_lt(min(ps$vals[i, ]), min(lo$vals[i, ]))
})

test_that("PEV boosters are scaled by primary coverage (not the whole population)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  age_mid <- default_age_lower() * 365 + 30
  epi <- malariasimulation::set_pev_epi(gp(), profile = malariasimulation::rtss_profile,
    timesteps = 1, coverages = 0.5, min_wait = 0, age = 180, booster_spacing = 360,
    booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  ps <- pev_series(epi, age_mid, 8 * 365)
  # the total FOI reduction can never exceed the primary coverage (0.5): only the
  # vaccinated can be protected, boosted or not. The old whole-population booster
  # layer let it reach ~0.52 (> 0.5); verify the coverage ceiling now holds.
  expect_lte(max(1 - ps$vals), 0.5 + 1e-9)
})

test_that("multi-drug resistance is coverage-weighted and drives dynamics", {
  skip_if_not_installed("malariasimulation")
  mk <- function(p, drug, art, etf, spc, tclear)
    malariasimulation::set_antimalarial_resistance(p, drug = drug, timesteps = c(0, 1000),
      artemisinin_resistance_proportion = c(0, art),
      partner_drug_resistance_proportion = c(0, 0),
      slow_parasite_clearance_probability = c(spc, spc),
      early_treatment_failure_probability = c(etf, etf),
      late_clinical_failure_probability = c(0, 0),
      late_parasitological_failure_probability = c(0, 0),
      reinfection_during_prophylaxis_probability = c(0, 0),
      slow_parasite_clearance_time = tclear)
  base_p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
    list(malariasimulation::AL_params, malariasimulation::SP_AQ_params))
  # coverages deliberately DON'T sum to 1, so the normalised share wn = cov/sum(cov)
  # differs from raw coverage (pins the /sum(cw) normalisation)
  base_p <- malariasimulation::set_clinical_treatment(base_p, drug = 1, timesteps = 1, coverages = 0.6)
  base_p <- malariasimulation::set_clinical_treatment(base_p, drug = 2, timesteps = 1, coverages = 0.3)
  w1 <- 0.6 / 0.9; w2 <- 0.3 / 0.9                     # normalised treatment shares
  # resistance on BOTH drugs, different arms -> the population value is the sum
  p <- mk(base_p, drug = 1, art = 0.5, etf = 0.4, spc = 0.2, tclear = 10)
  p <- mk(p, drug = 2, art = 1.0, etf = 0.1, spc = 0.3, tclear = 20)
  r <- resistance_series(p, translate_parameters(p))
  expect_equal(tail(r$etf, 1), w1 * 0.5 * 0.4 + w2 * 1.0 * 0.1, tolerance = 1e-9)
  expect_equal(tail(r$spc, 1), w1 * 0.5 * 0.2 + w2 * 1.0 * 0.3, tolerance = 1e-9)
  # rT_slow = 1 / coverage*spc-weighted mean clearance time (the rate feeding the ODE)
  s1 <- 0.5 * 0.2; s2 <- 1.0 * 0.3
  num <- w1 * s1 * 10 + w2 * s2 * 20; den <- w1 * s1 + w2 * s2
  expect_equal(r$rT_slow, 1 / (num / den), tolerance = 1e-9)
  # end-to-end: resistance raises prevalence vs the SAME treatment without resistance
  base <- run_simulation_ode(1500, base_p, init_EIR = 30)
  o <- run_simulation_ode(1500, p, init_EIR = 30)
  expect_true(all(is.finite(o$p_detect_lm_730_3650)))
  expect_gt(o$p_detect_lm_730_3650[1500], base$p_detect_lm_730_3650[1500])
})

test_that("resistance rT_slow uses the PEAK timestep, not the last", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                    list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.5)
  # resistance rises to 0.8 then is controlled back to 0 -> spc = (0, 0.4, 0)
  p <- malariasimulation::set_antimalarial_resistance(p, drug = 1, timesteps = c(0, 100, 1000),
    artemisinin_resistance_proportion = c(0, 0.8, 0),
    partner_drug_resistance_proportion = c(0, 0, 0),
    slow_parasite_clearance_probability = c(0.5, 0.5, 0.5),
    early_treatment_failure_probability = c(0, 0, 0),
    late_clinical_failure_probability = c(0, 0, 0),
    late_parasitological_failure_probability = c(0, 0, 0),
    reinfection_during_prophylaxis_probability = c(0, 0, 0),
    slow_parasite_clearance_time = 10)
  r <- resistance_series(p, translate_parameters(p))
  expect_equal(tail(r$spc, 1), 0)               # schedule ends at zero resistance
  expect_equal(r$rT_slow, 1 / 10, tolerance = 1e-9)   # from the peak (dt=10), not 1/rT_base
})

test_that("custom demography changes the equilibrium age structure", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  dr <- c(0.05, 0.008, 0.003, 0.004, 0.008, 0.02, 0.05, 0.12) / 365
  p <- malariasimulation::set_demography(gp(),
    agegroups = round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365),
    timesteps = 0, deathrates = matrix(dr, nrow = 1))
  inp <- build_inputs(p, 20)
  am <- inp$meta$age_mid; prop <- inp$meta$prop
  expect_equal(sum(prop), 1, tolerance = 1e-10)               # normalised
  # high infant/elderly mortality => fewer under-5s than constant-hazard (~0.21)
  expect_lt(sum(prop[am < 5 * 365]), 0.15)
  # default (constant hazard) still holds flat at equilibrium (offset 0 = exact seed)
  pflat <- gp(); pflat$acquired_immunity_offset <- 0
  o <- run_simulation_ode(400, pflat, init_EIR = 20)
  expect_lt(max(abs(o$EIR - o$EIR[1])), 1e-6)
})

test_that("invalid inputs error clearly", {
  skip_if_not_installed("malariasimulation")
  expect_error(run_simulation_ode(10, malariasimulation::get_parameters(), init_EIR = 0),
               "positive")
  expect_error(run_simulation_ode(10, malariasimulation::get_parameters(parasite = "vivax"),
                                  init_EIR = 10), "falciparum")
})

test_that("IRS reduces prevalence; zero-coverage IRS == baseline", {
  skip_if_not_installed("malariasimulation")
  base <- run_simulation_ode(1200, malariasimulation::get_parameters(), init_EIR = 20)
  irs <- function(cov) malariasimulation::set_spraying(
    malariasimulation::get_parameters(), timesteps = 365, coverages = cov,
    ls_theta = matrix(2.025, 1, 1), ls_gamma = matrix(-0.009, 1, 1),
    ks_theta = matrix(-2.222, 1, 1), ks_gamma = matrix(0.008, 1, 1),
    ms_theta = matrix(-1.232, 1, 1), ms_gamma = matrix(-0.009, 1, 1))
  z <- run_simulation_ode(1200, irs(0), init_EIR = 20)
  expect_lt(max(abs(z$p_detect_lm_730_3650 - base$p_detect_lm_730_3650)), 1e-9)
  o <- run_simulation_ode(1200, irs(0.8), init_EIR = 20)
  expect_lt(min(o$p_detect_lm_730_3650), base$p_detect_lm_730_3650[1])
})

test_that("carrying-capacity scaler = 1 is a no-op (flat)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_carrying_capacity(
    malariasimulation::get_parameters(), timesteps = 365,
    carrying_capacity_scalers = matrix(1, 1, 1))
  p$acquired_immunity_offset <- 0                        # flat-machinery check
  o <- run_simulation_ode(800, p, init_EIR = 20)
  expect_lt(max(abs(o$EIR - o$EIR[1])), 1e-6)
})

test_that("seasonal cycle is annually periodic (settled)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters(overrides = list(
    model_seasonality = TRUE, g0 = 2, g = c(0.3, 0.6, 0.9), h = c(0.1, 0.4, 0.7)))
  o <- run_simulation_ode(3650, p, init_EIR = 20)
  y3 <- mean(o$EIR[1096:1460]); y4 <- mean(o$EIR[1461:1825])
  expect_lt(abs(y3 - y4), 0.5)   # year-on-year mean stable (settled cycle)
})

test_that("kitchen-sink run is finite, conserved, and postie-consumable", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("postie")
  gp <- malariasimulation::get_parameters
  p <- gp(overrides = list(model_seasonality = TRUE))
  p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params,
                                            malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
  p <- malariasimulation::set_bednets(p, timesteps = 365, coverages = 0.6, retention = 3 * 365,
    dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1), rnm = matrix(0.24, 1, 1),
    gamman = 2.64 * 365)
  p <- malariasimulation::set_smc(p, drug = 2, timesteps = c(365, 395, 425),
    coverages = rep(0.9, 3), min_ages = rep(round(0.25 * 365), 3),
    max_ages = rep(round(5 * 365), 3))
  p <- malariasimulation::set_pev_epi(p, profile = malariasimulation::rtss_profile,
    timesteps = 1, coverages = 0.8, min_wait = 0, age = 180, booster_spacing = 360,
    booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  base <- run_simulation_ode(730, gp(), init_EIR = 20)
  o <- run_simulation_ode(730, p, init_EIR = 20)
  expect_true(all(is.finite(as.matrix(o[sapply(o, is.numeric)]))))
  expect_identical(names(o), names(base))               # column stability
  pop <- with(o, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-5)
  expect_true(all(o$p_detect_lm_730_3650 >= 0 & o$p_detect_lm_730_3650 <= 1))
  epi <- get_epi_outputs(o)
  expect_true(all(epi$prevalence$lm_prevalence_2_10 >= 0 &
                  epi$prevalence$lm_prevalence_2_10 <= 1))
  expect_true(all(is.finite(epi$rates$clinical) & epi$rates$clinical >= 0))
  expect_true(all(is.finite(epi$rates$dalys)))
})

test_that("timesteps = 1 returns two rows", {
  skip_if_not_installed("malariasimulation")
  o <- run_simulation_ode(1, malariasimulation::get_parameters(), init_EIR = 20)
  expect_equal(nrow(o), 2)
})
