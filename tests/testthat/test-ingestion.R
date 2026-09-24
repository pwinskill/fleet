# Parameter-ingestion completeness fixes: PEV multi-booster + multi-band, chemo-
# prevention age-band overlap + resistance, IRS accumulation, all-infection
# incidence output, eq_params, and the surfaced approximation warnings.

test_that("PEV: multiple boosters each contribute (not just the first)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  am <- default_age_lower() * 365 + 30
  epi <- function(nb) {
    bs <- if (nb == 2) c(360, 720) else 360
    malariasimulation::set_pev_epi(gp(), profile = malariasimulation::rtss_profile,
      timesteps = 1, coverages = 0.9, min_wait = 0, age = 180,
      booster_spacing = bs, booster_coverage = matrix(0.8, 1, nb),
      booster_profile = rep(list(malariasimulation::rtss_booster_profile), nb))
  }
  v1 <- pev_series(epi(1), am, 8 * 365)$vals
  v2 <- pev_series(epi(2), am, 8 * 365)$vals
  expect_gt(max(abs(v1 - v2)), 0.02)               # the 2nd booster changes protection
})

test_that("mass PEV protects every target age band, at every campaign", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  am <- default_age_lower() * 365 + 30
  mp <- malariasimulation::set_mass_pev(gp(), profile = malariasimulation::rtss_profile,
    timesteps = 365, coverages = 0.9, min_ages = c(0.5, 6) * 365, max_ages = c(5, 10) * 365,
    min_wait = 0, booster_spacing = 365, booster_coverage = matrix(0.8, 1, 1),
    booster_profile = list(malariasimulation::rtss_booster_profile))
  ps <- pev_series(mp, am, 5 * 365)
  # a week after efficacy starts, each band's cohort, now that much older
  g <- which(ps$times >= 365 + max(mp$pev_doses) + 7)[1]
  shift <- ps$times[g] - 365
  b1 <- which(am >= 0.5 * 365 + shift + 30 & am < 5 * 365 + shift)[1]
  b2 <- which(am >= 6 * 365 + shift + 30 & am < 10 * 365 + shift)[1]
  expect_gt(1 - ps$vals[b1, g], 0.02)              # band 1 protected
  expect_gt(1 - ps$vals[b2, g], 0.02)              # band 2 protected (was dropped before)
})

test_that("chemoprevention age band uses fractional overlap (narrow band not dropped)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  # PMC dose at 15 months -> a 30-day band that contains no age-group midpoint on
  # the default grid; overlap weighting must still clear a positive fraction.
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_pmc(p, drug = 1, timesteps = 1, coverages = 0.9,
                                  ages = round(1.25 * 365))
  inp <- build_inputs(p, 20, timesteps = 3 * 365)
  ev <- chemoprevention_events(p, 3 * 365)
  e <- ev[[1]]
  ov <- pmax(0, pmin(e$hi, inp$meta$age_hi) - pmax(e$lo, inp$meta$age_lo)) /
        (inp$meta$age_hi - inp$meta$age_lo)
  expect_gt(max(ov), 0)                            # some group overlaps (not silently dropped)
  expect_lte(max(ov), 1 + 1e-9)
})

test_that("IRS accumulates over rounds (closely-spaced rounds add protection)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  irs <- function(ts, cov) malariasimulation::set_spraying(gp(), timesteps = ts,
    coverages = cov,
    ls_theta = matrix(2.025, length(ts), 1), ls_gamma = matrix(-0.009, length(ts), 1),
    ks_theta = matrix(-2.222, length(ts), 1), ks_gamma = matrix(0.008, length(ts), 1),
    ms_theta = matrix(-1.232, length(ts), 1), ms_gamma = matrix(-0.009, length(ts), 1))
  a2 <- vector_control_series(irs(c(200, 260), c(0.9, 0.3)), 400)$a[1, ]
  a1 <- vector_control_series(irs(260, 0.3), 400)$a[1, ]
  expect_lt(min(a2), min(a1))   # a still-fresh high-coverage 1st round lowers biting further
})

test_that("all-infection incidence is output and exceeds clinical incidence", {
  skip_if_not_installed("malariasimulation")
  o <- run_simulation_ode(400, eqm(gp_bands(), 30))
  expect_true("n_inc_730_3650" %in% names(o))
  expect_true(all(o$n_inc_730_3650 >= o$n_inc_clinical_730_3650 - 1e-9))   # all >= clinical
  expect_true(all(is.finite(o$n_inc_0_36500)))
})

test_that("resistance on a chemoprevention drug reduces its clearance", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_smc(p, drug = 1, timesteps = 365, coverages = 0.9,
                                  min_ages = round(0.25 * 365), max_ages = round(5 * 365))
  p <- malariasimulation::set_antimalarial_resistance(p, drug = 1, timesteps = c(0, 1),
    artemisinin_resistance_proportion = c(0, 0.5),
    partner_drug_resistance_proportion = c(0, 0),
    slow_parasite_clearance_probability = c(0, 0),
    early_treatment_failure_probability = c(0.4, 0.4),
    late_clinical_failure_probability = c(0, 0),
    late_parasitological_failure_probability = c(0, 0),
    reinfection_during_prophylaxis_probability = c(0, 0),
    slow_parasite_clearance_time = 10)
  ev <- chemoprevention_events(p, 2 * 365)
  # cleared fraction = cov * eff * (1 - art*etf) = 0.9*eff*(1 - 0.5*0.4) = 0.9*eff*0.8
  expect_equal(ev[[1]]$frac, 0.9 * p$drug_efficacy[1] * 0.8, tolerance = 1e-9)
})

test_that("multi-drug first-line switch is modelled as a time-varying blend", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  # AL, then switch to SP-AQ at t=1000. The drug blend (cT/drug_eff/rP) must vary
  # in time, not freeze at peak-coverage weights (the old approximation, which used
  # to warn). No "peak-coverage" warning should be emitted now that it is modelled.
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::AL_params,
                                               malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = c(1, 1000),
                                                 coverages = c(0.6, 0))
  p <- malariasimulation::set_clinical_treatment(p, drug = 2, timesteps = c(1, 1000),
                                                 coverages = c(0, 0.6))
  wmsgs <- character(0)
  inp <- withCallingHandlers(
    build_inputs(p, 20, timesteps = 1200),
    warning = function(w) { wmsgs <<- c(wmsgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("peak-coverage", wmsgs)))
  # the odin model receives a genuine time-varying drug series spanning the switch
  expect_true(inp$pars$n_dmix > 1)
  # AL and SP-AQ differ in prophylaxis / efficacy, so the blended series changes
  expect_gt(length(unique(round(c(inp$pars$rP_vals, inp$pars$drug_eff_vals), 8))), 1)
})

test_that("custom eq_params are honoured by the seed", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  p$eq_params <- list(phi0 = 0.999)                # extreme, easy to detect
  inp <- build_inputs(p, 20)
  expect_equal(inp$meta$eqp$phi0, 0.999)
})
