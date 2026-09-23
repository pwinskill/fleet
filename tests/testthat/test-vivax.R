# P. vivax support.
#
# fleet replicates malariasimulation's parameter translation rather than calling
# it, so that this package is self-contained and version-robust (same reasoning
# as .back_translations for falciparum). A replicated table is only as good as
# the test that it still agrees with the original, so that is the first thing
# here -- and it is a real risk, not a theoretical one: the vivax table has one
# entry that is not a plain alias (phi_D_min = phi0 * phi1, because a
# malariasimulation vivax list stores phi1 as a ratio), and an alias that
# silently became a rename would change the equilibrium without erroring.

test_that("the replicated vivax translation agrees with malariasimulation's", {
  ms_tr <- get("translate_vivax_parameters", asNamespace("malariasimulation"))
  p <- malariasimulation::get_parameters(parasite = "vivax")

  a <- ms_tr(p)
  b <- translate_vivax_parameters(p)

  expect_setequal(names(a), names(b))
  for (nm in names(a)) {
    expect_equal(b[[nm]], a[[nm]],
                 info = paste("vivax translation differs for", nm))
  }
})

test_that("the vivax translation is additive, and converts phi1 to a level", {
  p <- malariasimulation::get_parameters(parasite = "vivax")
  b <- translate_vivax_parameters(p)

  # additive: every input field survives under its own name
  expect_true(all(names(p) %in% names(b)))
  expect_equal(b$sigma_squared, p$sigma_squared)
  expect_equal(b$n_heterogeneity_groups, p$n_heterogeneity_groups)

  # the one non-alias: a vivax list stores phi1 as phi_D_min/phi_D_max
  expect_equal(b$phi_D_max, p$phi0)
  expect_equal(b$phi_D_min, p$phi0 * p$phi1)
  expect_lt(b$phi_D_min, b$phi_D_max)

  # rates are inverted durations, not durations
  expect_equal(b$r_D, 1 / p$dd)
  expect_equal(b$r_par, 1 / p$ra)
})

test_that("the vivax equilibrium solves on fleet's own age grid", {
  # This is the seed fleet's vivax block will be built on, so it has to work at
  # fleet's resolution, not only at malariasimulation's 1000-point grid. Ages
  # are in YEARS here -- passing days returns NaN or a singular matrix rather
  # than erroring cleanly, which is a cheap mistake to make and an expensive
  # one to diagnose.
  skip_if_not_installed("malariaEquilibriumVivax")
  p <- translate_vivax_parameters(
    malariasimulation::get_parameters(parasite = "vivax"))
  age <- default_age_lower()

  foim <- vapply(c(0.1, 1, 5, 20, 120), function(E) {
    e <- malariaEquilibriumVivax::vivax_equilibrium(EIR = E, ft = 0.2, p = p, age = age)
    states <- e$states[c("S", "U", "A", "D", "T", "P")]
    expect_false(any(vapply(states, function(x) any(is.na(x)), logical(1))))
    # one unit of population, distributed over age x heterogeneity x batch
    expect_equal(sum(vapply(states, sum, numeric(1))), 1, tolerance = 1e-8)
    # the batch dimension is kmax + 1
    expect_equal(dim(e$states$S)[3], p$kmax + 1)
    e$FOIM
  }, numeric(1))

  expect_true(all(is.finite(foim)))
  expect_true(all(foim > 0))
  expect_false(is.unsorted(foim))     # FOIM rises with transmission
})

test_that("vivax shares its heterogeneity nodes with fleet and the IBM", {
  # fleet seeds each heterogeneity stratum separately for falciparum, but the
  # vivax equilibrium returns the stratum dimension already built. That is only
  # usable if the two discretisations are the same one.
  p <- malariasimulation::get_parameters(parasite = "vivax")
  s2 <- p$sigma_squared
  n <- p$n_heterogeneity_groups

  gq <- malariaEquilibrium::gq_normal(n)
  zeta_fleet <- exp(gq$nodes * sqrt(s2) - s2 / 2)

  q <- statmod::gauss.quad.prob(n, dist = "normal", mu = -0.5 * s2, sigma = sqrt(s2))
  zeta_eq <- exp(q$nodes)

  expect_equal(zeta_fleet, zeta_eq)
  expect_equal(gq$weights, q$weights)
})
