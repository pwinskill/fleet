test_that("model holds flat at the malariaEquilibrium equilibrium (ft = 0)", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("malariaEquilibrium")

  p <- malariasimulation::get_parameters()
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
  p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1,
                                                 coverages = 0.4)
  out <- run_simulation_ode(3650, p, init_EIR = 20)
  expect_lt(max(abs(out$EIR - out$EIR[1])), 1e-5)
  prev <- out$p_detect_lm_730_3650
  expect_lt(max(abs(prev - prev[1])), 1e-6)
  expect_equal(out$ft[1], 0.4)
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
