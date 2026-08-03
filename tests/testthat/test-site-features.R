# Features needed to run real site::site_parameters() inputs: logistic net
# retention, multi-distribution net usage, and time-varying demography.

test_that("logistic net-retention survival is correct (S(0)=1, S(half_life)=0.5, monotone)", {
  p <- list(bednet_logistic_half_life = 825, bednet_logistic_k = 20)
  expect_equal(.net_survival(0, p), 1)
  expect_equal(.net_survival(825, p), 0.5, tolerance = 1e-6)   # half-life => 50% retained
  sn <- seq(0, 2000, by = 25)
  s <- .net_survival(sn, p)
  expect_true(all(diff(s) <= 1e-12))                            # monotone non-increasing
  expect_true(all(s >= 0 & s <= 1))
  # hard cap on net lifetime: l = half_life / sqrt(1 - k/(k - log 0.5))
  l <- 825 / sqrt(1 - 20 / (20 - log(0.5)))
  expect_equal(.net_survival(l + 100, p), 0)                    # beyond max lifetime => 0
  expect_lt(.net_survival(0.8 * l, p), 1e-4)                    # deep tail is ~0
})

test_that("logistic-retention bed nets build and run (no exponential-only error)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_bednets(malariasimulation::get_parameters(),
    timesteps = 365, coverages = 0.6, retention = 3 * 365,
    dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
    rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)
  p$bednet_retention <- NULL                    # switch to logistic (as site params supply)
  p$bednet_logistic_half_life <- 3 * 365
  p$bednet_logistic_k <- 20
  expect_silent(vc <- vector_control_series(p, 2000))
  o <- run_simulation_ode(1200, p, init_EIR = 20)
  expect_true(all(is.finite(o$p_detect_lm_730_3650)))
  expect_lt(min(o$p_detect_lm_730_3650), o$p_detect_lm_730_3650[1])   # nets reduce prevalence
})

test_that("repeated net distributions accumulate more usage than a single one", {
  skip_if_not_installed("malariasimulation")
  base <- malariasimulation::set_bednets(malariasimulation::get_parameters(),
    timesteps = 365, coverages = 0.3, retention = 5 * 365,
    dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
    rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)
  rep3 <- malariasimulation::set_bednets(malariasimulation::get_parameters(),
    timesteps = c(365, 730, 1095), coverages = rep(0.3, 3), retention = 5 * 365,
    dn0 = matrix(0.387, 3, 1), rn = matrix(0.563, 3, 1),
    rnm = matrix(0.24, 3, 1), gamman = rep(2.64 * 365, 3))
  # at day 1200, three 30%-distributions have protected more people (lower biting
  # rate a) than one 30%-distribution back at day 365
  a1 <- vector_control_series(base, 1300)$a[1, ]
  a3 <- vector_control_series(rep3, 1300)$a[1, ]
  expect_lt(min(a3), min(a1))
})

test_that("time-varying demography: mu_age series varies over time and runs", {
  skip_if_not_installed("malariasimulation")
  dr <- matrix(c(1e-4, 1e-3,        # t = 0 baseline
                 5e-5, 5e-4), nrow = 2, byrow = TRUE)   # t = 10y (mortality halved)
  p <- malariasimulation::set_demography(malariasimulation::get_parameters(),
    agegroups = c(1825, 36500), timesteps = c(0, 10 * 365), deathrates = dr)
  inp <- build_inputs(p, 20, timesteps = 20 * 365)
  expect_equal(inp$pars$n_mut, 2)                                   # two demography knots
  expect_true(any(inp$pars$mu_age_z[, 1] != inp$pars$mu_age_z[, 2]))# genuinely time-varying
  expect_equal(inp$pars$mu_age_t, c(0, 10 * 365))
  o <- run_simulation_ode(15 * 365, p, init_EIR = 20)
  expect_true(all(is.finite(o$p_detect_lm_730_3650)))
})
