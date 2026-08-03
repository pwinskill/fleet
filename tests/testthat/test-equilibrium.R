test_that("model holds flat at the malariaEquilibrium equilibrium (ft = 0)", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("malariaEquilibrium")

  p <- malariasimulation::get_parameters()
  # default acq_offset is 0 (the validated best mean-field match), for which the seed IS
  # the exact fixed point; this verifies the numerical machinery (equilibrium seeding,
  # immunity flux-aging, Erlang-lag seeding) holds flat. Set explicitly for clarity.
  p$acquired_immunity_offset <- 0
  out <- run_simulation_ode(3650, p, init_EIR = 20)

  # adult EIR is reproduced and held flat
  expect_equal(out$EIR[1], 20, tolerance = 1e-3)
  expect_lt(max(abs(out$EIR - out$EIR[1])), 1e-5)

  # LM prevalence 2-10 is held flat
  prev <- out$p_detect_lm_730_3650
  expect_lt(max(abs(prev - prev[1])), 1e-6)

  # population conserved
  pop <- with(out, S_count + D_count + A_count + U_count + Tr_count + Ph_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
})

test_that("model holds flat at equilibrium with treatment (ft > 0)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  p$acquired_immunity_offset <- 0                    # test the machinery at the fixed point
  p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1,
                                                 coverages = 0.4)
  out <- run_simulation_ode(3650, p, init_EIR = 20)
  expect_lt(max(abs(out$EIR - out$EIR[1])), 1e-5)
  prev <- out$p_detect_lm_730_3650
  expect_lt(max(abs(prev - prev[1])), 1e-6)
  expect_equal(out$ft[1], 0.4)
})

test_that("acquired-immunity offset toggle: default 0 holds flat, 0.5 shifts like the IBM Hill", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  off <- run_simulation_ode(3650, p, init_EIR = 20)             # default offset 0 -> flat seed
  p_on <- p; p_on$acquired_immunity_offset <- 0.5
  on  <- run_simulation_ode(3650, p_on, init_EIR = 20)          # +0.5 -> literal IBM Hill calls
  # default (0) holds flat at the malariaEquilibrium seed
  prev_off <- off$p_detect_lm_730_3650
  expect_lt(max(abs(prev_off - prev_off[1])), 1e-6)
  # the +0.5 toggle is active: it relaxes off the (offset-free) seed, bounded and stable
  prev_on <- on$p_detect_lm_730_3650
  expect_gt(max(abs(prev_on - prev_on[1])), 1e-6)
  expect_true(all(is.finite(prev_on)) && all(prev_on > 0 & prev_on < 1))
  # +0.5 raises effective acquired immunity => lower severe susceptibility (theta)
  sev_col <- grep("^n_inc_severe_", names(on))[1]
  expect_lt(sum(on[[sev_col]]), sum(off[[sev_col]]))
})

test_that("PfPR increases monotonically with EIR (simulated, not seed)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  prev <- vapply(c(1, 10, 50, 200), function(e) {
    o <- run_simulation_ode(365, p, init_EIR = e)
    o$p_detect_lm_730_3650[nrow(o)]   # final (simulated) timestep
  }, numeric(1))
  expect_true(all(diff(prev) > 0))
})

test_that("output feeds postie without error", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("postie")
  p <- malariasimulation::get_parameters()
  out <- run_simulation_ode(730, p, init_EIR = 20)
  epi <- get_epi_outputs(out)
  expect_true(all(c("rates", "prevalence") %in% names(epi)))
  expect_true("lm_prevalence_2_10" %in% names(epi$prevalence))
  expect_true(all(c("clinical", "severe", "dalys") %in% names(epi$rates)))
})
