test_that("every parameter a draw sets reaches fleet's dynamics, at its live value", {
  skip_if_not_installed("malariasimulation")
  # set_parameter_draw() swaps in one draw of the core parameters from the model
  # fit's joint posterior. fleet has to follow each as the IBM does: move one at a
  # time to its value in a draw far from the median, and fleet's dynamic inputs
  # must move. A parameter fleet takes from somewhere else leaves them unchanged
  # -- as iv0 did, read from set_equilibrium()'s frozen copy, which holds
  # malariaEquilibrium's default whatever the draw.
  base <- malariasimulation::get_parameters()
  far <- malariasimulation::set_parameter_draw(base, 256)
  drawn <- setdiff(names(base)[!mapply(identical, base, far[names(base)])], "parameter_draw")
  expect_setequal(drawn, names(malariasimulation::parameter_draws_pf[[256]]))
  ref <- build_inputs(eqm(base, 20), 20, timesteps = 10)$pars
  # The seeded states move with any parameter that moves the equilibrium, so
  # they cannot show whether the dynamics follow one. Leave out everything that
  # moves with the seed EIR alone.
  seeded <- build_inputs(eqm(base, 21), 21, timesteps = 10)$pars
  dyn <- names(ref)[mapply(identical, ref, seeded[names(ref)])]
  # the inputs fleet takes unchanged from the parameter list, as the IBM reads them
  as_read <- c(b0 = "b0", ib0 = "ib0", kb = "kb", fd0 = "fd0", ad = "ad0", gammad = "gd",
               d1 = "d1", id0 = "id0", kd = "kd", phi0 = "phi0", phi1 = "phi1", ic0 = "ic0",
               kc = "kc", theta0 = "theta0", theta1 = "theta1", iv0 = "iv0", kv = "kv",
               fv0 = "fv0", av = "av", gammav = "gammav", cd = "cD", gamma1 = "g_inf",
               cu = "cU")
  for (nm in drawn) {
    p <- base; p[[nm]] <- far[[nm]]
    got <- build_inputs(eqm(p, 20), 20, timesteps = 10)$pars
    expect_false(identical(got[dyn], ref[dyn]),
                 info = paste(nm, "does not reach fleet's dynamic inputs"))
    if (nm %in% names(as_read))
      expect_equal(got[[as_read[[nm]]]], far[[nm]], info = paste(nm, "is not its live value"))
  }
})

test_that("a run under a parameter draw holds its seed, with no frozen-copy warning", {
  skip_if_not_installed("malariasimulation")
  # Draws 143 and 256 sit near the 5th and 95th percentiles of iv0. A parameter
  # the seed and the dynamics disagreed on would show as drift off the seed; a
  # live value overridden by the frozen copy would warn.
  for (d in c(143, 256)) {
    p <- eqm(malariasimulation::set_parameter_draw(malariasimulation::get_parameters(), d), 20)
    expect_no_warning(o <- run_simulation_ode(5 * 365, p))
    prev <- pfpr(o)
    expect_lt(max(abs(prev / prev[1] - 1)), 1e-2)
  }
})

test_that("draws with b0 near 1 run on the default grid, and nothing in them goes negative", {
  skip_if_not_installed("malariasimulation")
  # Draws 410 and 464 have b0 0.959 and 0.990. The default grid's 16.6-day infant
  # groups lose 0.06 of their people a day to ageing and death, so a check that
  # bounded a day's infections by b0, as if everyone were bitten every day,
  # refused both; but the people ageing into a group every day make good what its
  # infections take. They run, and the whole state stays non-negative every day.
  for (d in c(410, 464)) {
    for (eir in c(20, 120)) {
      p <- eqm(malariasimulation::set_parameter_draw(malariasimulation::get_parameters(), d), eir)
      expect_gt(p$b0, 0.95)
      expect_no_error(run_simulation_ode(365, p))
      expect_gte(scan_positivity(365, p)$lowest, 0)
    }
  }
})
