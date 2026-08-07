test_that("no intervention (multi-species) holds flat at equilibrium", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  p <- malariasimulation::set_species(
    p, list(malariasimulation::arab_params, malariasimulation::fun_params),
    c(0.6, 0.4)
  )
  p$acquired_immunity_offset <- 0; p$bite_dedup <- 0                    # flat-machinery check (exact seed)
  o <- run_simulation_ode(1000, p, init_EIR = 20)
  expect_lt(max(abs(o$EIR - o$EIR[1])) / o$EIR[1], 1e-2)
  expect_lt(max(abs(o$p_detect_lm_730_3650 - o$p_detect_lm_730_3650[1])) / o$p_detect_lm_730_3650[1], 1e-2)
})

test_that("zero-coverage bednets == no intervention", {
  skip_if_not_installed("malariasimulation")
  base <- run_simulation_ode(1000, malariasimulation::get_parameters(), init_EIR = 20)
  p <- malariasimulation::set_bednets(
    malariasimulation::get_parameters(), timesteps = 100, coverages = 0,
    retention = 3 * 365, dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
    rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365
  )
  nets0 <- run_simulation_ode(1000, p, init_EIR = 20)
  expect_lt(max(abs(nets0$p_detect_lm_730_3650 - base$p_detect_lm_730_3650)), 1e-9)
})

test_that("bednets reduce prevalence then revert as they decay", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_bednets(
    malariasimulation::get_parameters(), timesteps = 365, coverages = 0.8,
    retention = 2 * 365, dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
    rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365
  )
  o <- run_simulation_ode(2000, p, init_EIR = 20)
  base <- o$p_detect_lm_730_3650[1]
  expect_lt(min(o$p_detect_lm_730_3650), base)          # nets reduce prevalence
  expect_gt(o$p_detect_lm_730_3650[2001], min(o$p_detect_lm_730_3650))  # reverts
})

test_that("higher net coverage gives lower prevalence (monotone)", {
  skip_if_not_installed("malariasimulation")
  prev <- vapply(c(0, 0.4, 0.8), function(cov) {
    p <- malariasimulation::set_bednets(
      malariasimulation::get_parameters(), timesteps = 200, coverages = cov,
      retention = 5 * 365, dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
      rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365
    )
    run_simulation_ode(1200, p, init_EIR = 20)$p_detect_lm_730_3650[1201]
  }, numeric(1))
  expect_true(all(diff(prev) < 0))
})

test_that("time-varying treatment lowers prevalence and ft column tracks schedule", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                    list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1000,
                                                 coverages = 0.6)
  o <- run_simulation_ode(2000, p, init_EIR = 20)
  expect_equal(o$ft[500], 0)          # before treatment
  expect_equal(o$ft[2001], 0.6)       # after treatment
  expect_lt(o$p_detect_lm_730_3650[2001], o$p_detect_lm_730_3650[500])
})

test_that("antimalarial resistance raises prevalence at fixed coverage", {
  skip_if_not_installed("malariasimulation")
  base_p <- malariasimulation::set_clinical_treatment(
    malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                 list(malariasimulation::AL_params)),
    drug = 1, timesteps = 1, coverages = 0.6)
  res_p <- malariasimulation::set_antimalarial_resistance(
    base_p, drug = 1, timesteps = 1,
    artemisinin_resistance_proportion = 0.8, partner_drug_resistance_proportion = 0,
    slow_parasite_clearance_probability = 0.5, early_treatment_failure_probability = 0.5,
    late_clinical_failure_probability = 0, late_parasitological_failure_probability = 0,
    reinfection_during_prophylaxis_probability = 0, slow_parasite_clearance_time = 20)
  base <- run_simulation_ode(1825, base_p, init_EIR = 20)$p_detect_lm_730_3650[1826]
  res <- run_simulation_ode(1825, res_p, init_EIR = 20)$p_detect_lm_730_3650[1826]
  expect_gt(res, base)
})

test_that("MDA clears infection (prevalence dips, population conserved)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                    list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_mda(p, drug = 1, timesteps = c(365, 730),
    coverages = c(0.8, 0.8), min_ages = rep(round(0.25 * 365), 2),
    max_ages = rep(round(5 * 365), 2))
  o <- run_simulation_ode(1200, p, init_EIR = 20)
  pop <- with(o, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)         # conserved
  expect_lt(min(o$p_detect_lm_730_3650), o$p_detect_lm_730_3650[1])  # MDA reduces
})

test_that("single-round MDA runs (regression: no rep() self-recursion)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_mda(
    malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                 list(malariasimulation::SP_AQ_params)),
    drug = 1, timesteps = 365, coverages = 0.8,
    min_ages = round(0.25 * 365), max_ages = round(5 * 365))
  o <- run_simulation_ode(800, p, init_EIR = 20)
  pop <- with(o, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
  expect_lt(min(o$p_detect_lm_730_3650), o$p_detect_lm_730_3650[1])
})

test_that("SMC reduces transmission and conserves population", {
  skip_if_not_installed("malariasimulation")
  smc <- malariasimulation::set_smc(
    malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                 list(malariasimulation::SP_AQ_params)),
    drug = 1, timesteps = c(365, 395, 425), coverages = rep(0.9, 3),
    min_ages = rep(round(0.25 * 365), 3), max_ages = rep(round(5 * 365), 3))
  o <- run_simulation_ode(800, smc, init_EIR = 20)
  pop <- with(o, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
  expect_lt(min(o$p_detect_lm_730_3650), o$p_detect_lm_730_3650[1])
})

test_that("PMC generates continuous-cadence events and runs finite", {
  skip_if_not_installed("malariasimulation")
  pmc <- malariasimulation::set_pmc(
    malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                 list(malariasimulation::SP_AQ_params)),
    drug = 1, timesteps = c(1, 365), coverages = c(0.9, 0.9),
    ages = round(c(2.5, 3.5, 9) * 30))
  ev <- chemoprevention_events(pmc, 730)
  expect_gt(length(ev), 20)          # monthly cadence, not just 2 timesteps
  o <- run_simulation_ode(730, pmc, init_EIR = 20)
  expect_true(all(is.finite(o$p_detect_lm_730_3650)))
  pop <- with(o, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
})

test_that("zero-coverage MDA == no intervention", {
  skip_if_not_installed("malariasimulation")
  base <- run_simulation_ode(600, malariasimulation::get_parameters(), init_EIR = 20)
  p <- malariasimulation::set_mda(
    malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                 list(malariasimulation::SP_AQ_params)),
    drug = 1, timesteps = 200, coverages = 0,
    min_ages = round(0.25 * 365), max_ages = round(5 * 365))
  z <- run_simulation_ode(600, p, init_EIR = 20)
  expect_lt(max(abs(z$p_detect_lm_730_3650 - base$p_detect_lm_730_3650)), 1e-9)
})

test_that("PEV (EPI and mass) reduces prevalence", {
  skip_if_not_installed("malariasimulation")
  epi <- malariasimulation::set_pev_epi(
    malariasimulation::get_parameters(), profile = malariasimulation::rtss_profile,
    timesteps = 1, coverages = 0.9, min_wait = 0, age = round(6 * 30),
    booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  o <- run_simulation_ode(2000, epi, init_EIR = 20)
  expect_lt(o$p_detect_lm_730_3650[2001], o$p_detect_lm_730_3650[1])

  mass <- malariasimulation::set_mass_pev(
    malariasimulation::get_parameters(), profile = malariasimulation::rtss_profile,
    timesteps = 365, coverages = 0.9, min_ages = round(0.5 * 365),
    max_ages = round(15 * 365), min_wait = 0, booster_spacing = round(12 * 30),
    booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  om <- run_simulation_ode(1200, mass, init_EIR = 20)
  expect_lt(min(om$p_detect_lm_730_3650), om$p_detect_lm_730_3650[1] - 0.05)
  pop <- with(om, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
})

test_that("TBV reduces onward transmission (EIR and prevalence fall)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_tbv(malariasimulation::get_parameters(),
                                  timesteps = 365, coverages = 0.9, ages = 5:15)
  o <- run_simulation_ode(1200, p, init_EIR = 20)
  expect_lt(min(o$EIR), o$EIR[1])
  expect_lt(min(o$p_detect_lm_730_3650), o$p_detect_lm_730_3650[1])
})

test_that("zero-coverage vaccines == no intervention (PEV and TBV)", {
  skip_if_not_installed("malariasimulation")
  base <- run_simulation_ode(800, malariasimulation::get_parameters(), init_EIR = 20)
  pe <- malariasimulation::set_mass_pev(
    malariasimulation::get_parameters(), profile = malariasimulation::rtss_profile,
    timesteps = 200, coverages = 0, min_ages = round(2 * 365), max_ages = round(10 * 365),
    min_wait = 0, booster_spacing = round(12 * 30), booster_coverage = matrix(0, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  z <- run_simulation_ode(800, pe, init_EIR = 20)
  expect_lt(max(abs(z$p_detect_lm_730_3650 - base$p_detect_lm_730_3650)), 1e-9)
  tb <- malariasimulation::set_tbv(malariasimulation::get_parameters(),
                                   timesteps = 200, coverages = 0, ages = 5:15)
  zt <- run_simulation_ode(800, tb, init_EIR = 20)
  expect_lt(max(abs(zt$p_detect_lm_730_3650 - base$p_detect_lm_730_3650)), 1e-9)
})

test_that("mass PEV: higher coverage -> larger reduction; only target ages affected", {
  skip_if_not_installed("malariasimulation")
  mk <- function(cov) malariasimulation::set_mass_pev(
    malariasimulation::get_parameters(), profile = malariasimulation::rtss_profile,
    timesteps = 200, coverages = cov, min_ages = round(2 * 365), max_ages = round(10 * 365),
    min_wait = 0, booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  dip <- vapply(c(0, 0.4, 0.8), function(cov)
    min(run_simulation_ode(900, mk(cov), init_EIR = 20)$p_detect_lm_730_3650), numeric(1))
  expect_true(all(diff(dip) < 0))
  # age-targeting: mass PEV on 2-10y reduces FOI only in that band, ever
  am <- build_inputs(mk(0.8), 20)$meta$age_mid
  pv <- pev_series(mk(0.8), am, 900)$vals
  in_band <- am >= 2 * 365 & am < 10 * 365
  expect_true(all(pv[!in_band, ] == 1))       # out-of-band never reduced
  expect_true(any(pv[in_band, ] < 1))         # in-band reduced at some time
})

test_that("seasonality: aseasonal is flat, seasonal oscillates around the mean", {
  skip_if_not_installed("malariasimulation")
  p0 <- malariasimulation::get_parameters(); p0$acquired_immunity_offset <- 0; p0$bite_dedup <- 0
  o0 <- run_simulation_ode(1000, p0, init_EIR = 20)
  expect_lt(max(abs(o0$EIR - o0$EIR[1])) / o0$EIR[1], 1e-2)          # aseasonal flat (machinery)
  p <- malariasimulation::get_parameters(overrides = list(
    model_seasonality = TRUE, g0 = 2, g = c(0.3, 0.6, 0.9), h = c(0.1, 0.4, 0.7)))
  o <- run_simulation_ode(1825, p, init_EIR = 20)
  eir_y4 <- o$EIR[1461:1825]
  expect_gt(max(eir_y4) - min(eir_y4), 5)               # genuine seasonal swing
  expect_true(abs(mean(eir_y4) - 20) < 4)               # mean near target EIR
})

test_that("flexible carrying capacity scales transmission", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_carrying_capacity(
    malariasimulation::get_parameters(), timesteps = 365,
    carrying_capacity_scalers = matrix(0.5, 1, 1))
  o <- run_simulation_ode(1460, p, init_EIR = 20)
  expect_lt(o$EIR[1460], o$EIR[1] - 1)                  # halving K lowers EIR
})
