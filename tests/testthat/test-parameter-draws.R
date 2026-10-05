test_that("every parameter a draw sets reaches fleet's inputs", {
  skip_if_not_installed("malariasimulation")
  # set_parameter_draw() swaps in one draw of the core parameters from the model
  # fit's joint posterior. fleet has to follow each as the IBM does: move one at a
  # time to its value in a draw far from the median, and fleet's inputs must move.
  # A parameter fleet takes from somewhere else leaves them unchanged -- as iv0
  # did, read from set_equilibrium()'s frozen copy, which never carries it.
  base <- malariasimulation::get_parameters()
  far <- malariasimulation::set_parameter_draw(base, 256)
  drawn <- setdiff(names(base)[!mapply(identical, base, far[names(base)])], "parameter_draw")
  expect_gt(length(drawn), 20)
  ref <- build_inputs(eqm(base, 20), 20, timesteps = 10)$pars
  for (nm in drawn) {
    p <- base; p[[nm]] <- far[[nm]]
    got <- build_inputs(eqm(p, 20), 20, timesteps = 10)$pars
    expect_false(identical(got, ref), info = paste(nm, "does not reach fleet's inputs"))
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
