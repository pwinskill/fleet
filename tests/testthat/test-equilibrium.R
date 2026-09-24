test_that("model stays near the malariaEquilibrium seed (ft = 0)", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("malariaEquilibrium")

  p <- malariasimulation::get_parameters()
  # default acq_offset is 0 (the validated best mean-field match), for which the seed IS
  # the least-drifting configuration; this verifies the numerical machinery (seeding,
  # immunity transport with ageing, the delay lines) holds flat. Set explicitly for
  # clarity.
  p$acquired_immunity_offset <- 0; p$bite_dedup <- 0
  out <- run_simulation_ode(3650, eqm(p, 20))

  # adult EIR is reproduced and held flat
  expect_equal(out$EIR[1], 20, tolerance = 1e-3)
  expect_lt(max(abs(out$EIR - out$EIR[1])) / out$EIR[1], 1e-2)

  # LM prevalence 2-10 is held flat
  prev <- pfpr(out)
  expect_lt(max(abs(prev - prev[1])) / prev[1], 1e-2)

  # population conserved
  pop <- with(out, S_count + D_count + A_count + U_count + Tr_count)
  expect_equal(max(pop), min(pop), tolerance = 1e-6)
})

test_that("model stays near the malariaEquilibrium seed with treatment in force at the seed", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  p$acquired_immunity_offset <- 0; p$bite_dedup <- 0                    # test the machinery at the fixed point
  p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
  # from timestep 0, so the IBM's two seeding clocks agree (treatment from timestep
  # 1 starts the humans untreated, as the IBM does: see test-daily-mechanics.R)
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 0,
                                                 coverages = 0.4)
  out <- run_simulation_ode(3650, eqm(p, 20))
  expect_lt(max(abs(out$EIR - out$EIR[1])) / out$EIR[1], 1e-2)
  prev <- pfpr(out)
  expect_lt(max(abs(prev - prev[1])) / prev[1], 1e-2)
  expect_equal(out$ft[1], 0.4)
})

test_that("acquired-immunity offset toggle: default 0 holds flat, 0.5 shifts like the IBM Hill", {
  skip_if_not_installed("malariasimulation")
  p <- gp_bands()
  p$bite_dedup <- 0                          # isolate the offset: linear FOI, least drift
  off <- run_simulation_ode(3650, eqm(p, 20))             # default offset 0 -> flat seed
  p_on <- p; p_on$acquired_immunity_offset <- 0.5
  on  <- run_simulation_ode(3650, eqm(p_on, 20))          # +0.5 -> literal IBM Hill calls
  # default (0) holds flat at the malariaEquilibrium seed
  prev_off <- pfpr(off)
  expect_lt(max(abs(prev_off - prev_off[1])) / prev_off[1], 1e-2)
  # the +0.5 toggle is active: it relaxes off the (offset-free) seed, bounded and stable
  prev_on <- pfpr(on)
  expect_gt(max(abs(prev_on - prev_on[1])), 1e-6)
  expect_true(all(is.finite(prev_on)) && all(prev_on > 0 & prev_on < 1))
  # +0.5 raises effective acquired immunity => lower severe susceptibility (theta)
  sev_col <- grep("^n_inc_severe_", names(on))[1]
  expect_lt(sum(on[[sev_col]]), sum(off[[sev_col]]))
})

test_that("bite deduplication saturates the infection hazard at high exposure", {
  skip_if_not_installed("malariasimulation")
  gp <- gp_bands
  # High EIR => expected bites/person/day approaches 1, where the IBM's one-infection-
  # per-timestep cap binds: the saturating hazard must give LESS infection (and so less
  # clinical incidence) than the unbounded linear form.
  p_sat <- gp(); p_sat$bite_dedup <- 1
  p_lin <- gp(); p_lin$bite_dedup <- 0
  sat <- run_simulation_ode(1825, eqm(p_sat, 200))
  lin <- run_simulation_ode(1825, eqm(p_lin, 200))
  cc <- grep("^n_inc_clinical_", names(sat))[1]
  rows <- 1000:1825                                  # settled window
  expect_lt(sum(sat[[cc]][rows]), sum(lin[[cc]][rows]))
  # ... and negligibly different at low exposure, where 1 - exp(-EPS) ~ EPS
  sat_lo <- run_simulation_ode(1095, eqm(p_sat, 0.5))
  lin_lo <- run_simulation_ode(1095, eqm(p_lin, 0.5))
  r <- sum(sat_lo[[cc]]) / sum(lin_lo[[cc]])
  expect_gt(r, 0.98); expect_lte(r, 1.001)
})

test_that("PfPR increases monotonically with EIR (simulated, not seed)", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  prev <- vapply(c(1, 10, 50, 200), function(e) {
    o <- run_simulation_ode(365, eqm(p, e))
    pfpr(o)[nrow(o)]   # final (simulated) timestep
  }, numeric(1))
  expect_true(all(diff(prev) > 0))
})

test_that("output feeds postie without error", {
  skip_if_not_installed("malariasimulation")
  skip_if_not_installed("postie")
  # postie derives clinical and severe rates band by band, so the list renders
  # both over the same bands, as a site file's does
  lo <- c(0, 5, 15) * 365; hi <- c(5, 15, 100) * 365 - 1
  p <- malariasimulation::get_parameters(list(
    clinical_incidence_rendering_min_ages = lo, clinical_incidence_rendering_max_ages = hi,
    severe_incidence_rendering_min_ages = lo, severe_incidence_rendering_max_ages = hi))
  out <- run_simulation_ode(730, eqm(p, 20))
  # the count table is malariasimulation-shaped, so postie consumes it directly,
  # exactly as it consumes an IBM run
  prevalence <- postie::get_prevalence(out, diagnostic = "lm")
  rates <- postie::get_rates(out)
  expect_true("lm_prevalence_2_10" %in% names(prevalence))
  expect_true(all(c("clinical", "severe", "dalys") %in% names(rates)))
})
