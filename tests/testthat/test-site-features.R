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
  o <- run_simulation_ode(1200, eqm(p, 20))
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
  o <- run_simulation_ode(15 * 365, eqm(p, 20))
  expect_true(all(is.finite(o$p_detect_lm_730_3650)))
})

test_that("the live ageing and FOIM coefficients equal the frozen ones at the seed", {
  skip_if_not_installed("malariasimulation")
  # Immunity ageing uses the per-capita inflow into an age cell,
  # ain = r_age[i-1] * N[i-1,j] / N[i,j], read off the LIVE population, and the
  # FOIM denominator is the live sum(psi*N)/sum(N). Both were once frozen at the
  # seed as re[i] = r_age[i] + mu_age[i] and a `mean_psi` parameter, which is
  # exact only while the age structure is stationary.
  #
  # This is the compatibility half of that change: at the seed the two forms
  # must agree to machine precision, or every validated result moves. It is
  # what lets the fix be a no-op under constant death rates.
  inp <- build_inputs(eqm(malariasimulation::get_parameters(
    list(human_population = 10000)), 20), 20)
  pr <- inp$pars
  N <- pr$S0 + pr$D0 + pr$A0 + pr$U0 + pr$Tr0 +
    apply(pr$Ph0, c(1, 2), sum) + apply(pr$Phc0, c(1, 2), sum)

  mu0 <- pr$mu_age_z[, 1]
  re <- pr$r_age + mu0
  ain <- rbind(sum(mu0 * rowSums(N)) * pr$het_wt / N[1, ],       # births at i = 1
               pr$r_age[-pr$n_age] * N[-pr$n_age, ] / N[-1, ])
  expect_equal(ain, matrix(re, pr$n_age, pr$n_het), tolerance = 1e-12,
               ignore_attr = TRUE)

  # and the FOIM normalisation, which is the population mean of psi. Weighted by
  # psi alone: zeta averages 1 within every age band, so it cancels.
  expect_equal(sum(pr$psi * rowSums(N)) / sum(N),
               sum(inp$meta$prop * pr$psi), tolerance = 1e-12)

  # `mean_psi` is no longer a parameter; the model computes it from Npop
  expect_false("mean_psi" %in% names(pr))
})

test_that("time-varying mortality drives the immunity ageing coefficient", {
  skip_if_not_installed("malariasimulation")
  # The test above pins the two forms AGREEING at the seed, which is a property
  # of a stationary age structure and would still hold if the coefficient were
  # reverted to the frozen re[i]. This one pins the difference: the same
  # mortality step as the test above, run past it, where the structure is still
  # moving and the two forms genuinely disagree.
  #
  # Reverting the coefficient to the frozen re[i] = r_age[i] + mu_age[i], and
  # recompiling, gives a final prevalence of 0.6443813 and a final EIR of
  # 36.57542 instead of the values below -- 2.5e-4 and 7.6e-4 relative, which
  # is small but is thousands of times the tolerance asserted here. Measured by
  # actually reverting, not assumed: restoring the coefficient returns these
  # numbers bit for bit. The gap is modest on this two-band fixture; on a real
  # site-file mortality series it reaches several per cent.
  dr <- matrix(c(1e-4, 1e-3, 5e-5, 5e-4), nrow = 2, byrow = TRUE)
  p <- malariasimulation::set_demography(malariasimulation::get_parameters(),
    agegroups = c(1825, 36500), timesteps = c(0, 10 * 365), deathrates = dr)
  o <- run_simulation_ode(15 * 365, eqm(p, 20))
  n <- nrow(o)

  # the step really does move the age structure, or the test proves nothing
  expect_gt(abs(o$EIR[n] / o$EIR[1] - 1), 0.05)
  expect_equal(o$p_detect_lm_730_3650[1], 0.6476852095, tolerance = 1e-7)
  expect_equal(o$p_detect_lm_730_3650[n], 0.6442193184, tolerance = 1e-7)
  expect_equal(o$EIR[n], 36.60322457, tolerance = 1e-7)
})
