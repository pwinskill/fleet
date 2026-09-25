# Absolute reference values -- the regression net under the rest of this suite.
#
# WHAT this is for. Every other test here asserts a *relation*: flat at the
# equilibrium seed, monotone in EIR, intervention-on below intervention-off. A
# relation is blind to scale, so a mutation that changed the model's answer could
# pass all of them. Two did: doubling `sev_inc_a` (severe episodes exactly
# twice as many) and reverting the A-term correction in `inc_a` (+2.7% on
# all-infection incidence) each failed zero assertions in the whole suite.
#
# This file closes that hole. It pins the model's ABSOLUTE epidemiological output
# levels -- one quantity from every output family, at low, medium and high
# transmission -- to 1e-6, so any change that moves a number announces itself
# here instead of passing silently. It asserts nothing about whether the numbers
# are *right*; agreement with the IBM is fleetcheck's job. Its only job is
# to make a change visible.
#
# HOW to regenerate, when a model change is intended and the new numbers are the
# ones to keep:
#
#   FLEET_REGENERATE_REFERENCE=1 Rscript -e 'testthat::test_file(
#     "tests/testthat/test-reference.R", package = "fleet", load_package = "installed")'
#
# against a fresh R CMD INSTALL (devtools::test() compiles a debug build into
# src/, see CONTRIBUTING.md), which rewrites the reference CSVs from the current
# model and then trivially passes. Do it deliberately, and read the CSV diff: it is the
# record of exactly what the change did to the model's output. Never regenerate
# to make a red test go green.

## The pinned scenario. Nothing here may drift, or the CSV stops meaning anything:
## default parameters, default `ode_tuning()`, three decades of transmission, and a
## horizon long enough for the run to leave the seed and settle again.
REF_DAYS <- 730
REF_EIRS <- c(1, 20, 120)

## One quantity per output family, chosen so that no family can move unnoticed:
## PfPR 2-10 and its underlying LM count, PCR detection, clinical incidence in
## under-5s and all-ages, severe incidence in 2-10s and all-ages, all-infection
## incidence in 2-10s and all-ages, the two transmission scalars (FOIM onto
## mosquitoes, EIR back), the population infection-state counts, and the age-band
## denominators everything above is a count over.
##
## Why S/D/A/U are here. They were the last family in the run with no absolute
## assertion anywhere in the suite. The only other test that touches them is a
## conservation check, S + D + A + U + Tr + Ph = population, and that sum is
## scale-invariant: every possible redistribution *among* the states passes it. So
## the states could all move and nothing in the suite would say so -- exactly the
## hole this file exists to close.
##
## Why the other three state columns are not here. ft, Tr_count and Ph_count are
## absent or identically 0 across this scenario (it treats nobody, so nothing
## enters Tr or Ph). A pinned 0 would assert nothing about the model and would
## read as coverage the file does not have. They belong in a scenario with
## treatment switched on, or nowhere.
##
## Why only two of the n_age bands. n_age_* is emitted over the union of every
## rendering band. The 2-10 one is already pinned to 1e-6 in effect -- it is
## n_detect_lm_730_3650 over PfPR, and it is constant in time so its sum is
## implied too -- and the others are not implied by anything: n_age_0_1825 pins
## the demography's under-5 share and n_age_0_36500 pins the population scaling
## that turns every fraction in the model into the counts above.
REF_QUANTITIES <- c(
  "n_detect_lm_730_3650", "n_detect_pcr_730_3650",
  "n_inc_clinical_0_1825", "n_inc_clinical_0_36500",
  "n_inc_severe_730_3650", "n_inc_severe_0_36500",
  "n_inc_0_36500", "FOIM", "EIR",
  "S_count", "D_count", "A_count", "U_count",
  "n_age_0_1825", "n_age_0_36500", "n_inc_730_3650")

## Default parameters, rendering clinical, severe and all-infection incidence as
## well as the default 2-10 prevalence: malariasimulation renders only what the
## list asks for, and so does fleet.
reference_parameters <- function() gp_bands()

## Two deterministic reductions per series: the final day pins the LEVEL the run
## settles to, the sum over every output day pins the TRAJECTORY that got there.
## A change to either -- a shifted equilibrium or a different approach to it -- is
## caught, and a change that moved the level while preserving the sum cannot hide.
reference_summary <- function() {
  rows <- list()
  for (eir in REF_EIRS) {
    out <- run_simulation_ode(REF_DAYS, eqm(reference_parameters(), eir))
    expect_equal(nrow(out), REF_DAYS)                    # one row per day, 1..REF_DAYS
    expect_true(all(REF_QUANTITIES %in% names(out)))
    for (q in REF_QUANTITIES) {
      v <- out[[q]]
      rows[[length(rows) + 1]] <- data.frame(init_EIR = eir, quantity = q,
                                             statistic = "final", value = v[length(v)])
      rows[[length(rows) + 1]] <- data.frame(init_EIR = eir, quantity = q,
                                             statistic = "sum", value = sum(v))
    }
  }
  do.call(rbind, rows)
}

test_that("model output matches the committed reference values", {
  skip_if_not_installed("malariasimulation")

  got <- reference_summary()

  if (nzchar(Sys.getenv("FLEET_REGENERATE_REFERENCE"))) {
    # Deliberate regeneration (see the header). Full double precision: the test
    # asserts to 1e-6, so the file must not be the thing that loses digits.
    write.csv(transform(got, value = vapply(got$value, format, character(1), digits = 17)),
              test_path("reference-values.csv"), row.names = FALSE, quote = FALSE)
    message("FLEET_REGENERATE_REFERENCE: rewrote reference-values.csv from the current model")
  }

  ref <- utils::read.csv(test_path("reference-values.csv"))
  # the CSV must cover exactly the scenario above -- a quantity quietly dropped
  # from the file would otherwise silently stop being pinned
  key <- function(d) paste(d$init_EIR, d$quantity, d$statistic)
  expect_setequal(key(ref), key(got))
  ref <- ref[match(key(got), key(ref)), ]

  for (i in seq_len(nrow(got))) {
    expect_equal(got$value[i], ref$value[i], tolerance = 1e-6,
                 label = sprintf("%s [%s] at init_EIR = %g",
                                 got$quantity[i], got$statistic[i], got$init_EIR[i]))
  }
})


## ---- Intervention reference -------------------------------------------------
#
# WHY a second file. The scenario above has no interventions, so it is blind to
# the two defects that moved the most: vector control acted up to 10 days BEFORE
# deployment (linear interpolation with no knot at onset - 1, EIR 66% down the day
# before any net existed), and TBV / mass-PEV age selection matched the age-group
# MIDPOINT, so an age set that hit no midpoint was a silent no-op -- `ages = 18:20`
# returned a run bit-identical to no TBV at all. Neither could ever move a value
# in reference-values.csv. This block is the net under them.
#
# It is a separate CSV, not extra rows in the first one, so the transmission
# baseline stays byte-identical to what it was committed as and the two concerns
# can be regenerated independently. Same regeneration switch.

REF_INT_DAYS <- 400
NET_DAY <- 100

## The days are the point. `before` is the deployment day itself: a round on day
## 100 lands at the end of it, so row 100 MUST equal the no-intervention run, and
## it is what goes red if a deployment ever acts a day early. Each scenario is then
## sampled on the first day its effect reaches the sampled quantities and again
## later, so the samples pin that the effect arrives, and at the right size. Nets
## cut the biting rate from day 101, and humans feel it de = 12 days later (row
## 113). TBV cuts onward infectivity from day 101; mosquitoes read it delay_gam
## later, incubate for 9 days, and their bites reach humans de after that (row
## 135). Mass PEV protects nobody until its last primary dose, pev_doses = 0, 45,
## 90 days, so from row 191.
REF_INT_DAYS_SAMPLED <- list(
  nets = c(before = NET_DAY, onset = NET_DAY + 13, after = NET_DAY + 18),
  tbv_offgrid = c(before = NET_DAY, onset = NET_DAY + 35, after = NET_DAY + 65),
  mass_pev_split = c(before = NET_DAY, onset = NET_DAY + 91, after = NET_DAY + 136))

reference_int_scenarios <- function() {
  gp <- function() gp_bands()
  list(
    ## Onset timing. Nets are the cleanest probe: one campaign, a hard start date.
    nets = malariasimulation::set_bednets(
      gp(), timesteps = NET_DAY, coverages = 0.8, retention = 5 * 365,
      dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
      rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365),

    ## Fractional age weighting. 18:20 hits no age-group midpoint on the default
    ## grid (above 14 the midpoints are 14.4, 15.8, 17.5, 19.2, 21.3, ...), so under the old rule
    ## this scenario was indistinguishable from no TBV.
    tbv_offgrid = malariasimulation::set_tbv(
      gp(), timesteps = NET_DAY, coverages = 0.9, ages = c(18, 19, 20)),

    ## Band folding. Two adjacent bands of ONE campaign that together tile
    ## [0, 1825) must give exactly what the single band [0, 1825) gives: the
    ## straddling age group's two weights sum to 1, so splitting is a no-op.
    ## Folding them through .combine() instead of adding cost 6.3% here.
    mass_pev_split = malariasimulation::set_mass_pev(
      gp(), profile = malariasimulation::rtss_profile,
      timesteps = NET_DAY, coverages = 1,
      min_ages = c(0, 1195), max_ages = c(1195, 1825), min_wait = 0,
      booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
      booster_profile = list(malariasimulation::rtss_booster_profile)))
}

REF_INT_QUANTITIES <- c("EIR", "n_detect_lm_730_3650", "n_inc_clinical_730_3650")

reference_int_summary <- function() {
  rows <- list()
  for (nm in names(reference_int_scenarios())) {
    out <- run_simulation_ode(REF_INT_DAYS, eqm(reference_int_scenarios()[[nm]], 20))
    expect_true(all(REF_INT_QUANTITIES %in% names(out)))
    for (q in REF_INT_QUANTITIES) {
      days <- REF_INT_DAYS_SAMPLED[[nm]]
      for (s in names(days)) {
        rows[[length(rows) + 1]] <- data.frame(
          scenario = nm, quantity = q, statistic = s,
          value = out[[q]][out$timestep == days[[s]]])
      }
      rows[[length(rows) + 1]] <- data.frame(
        scenario = nm, quantity = q, statistic = "final",
        value = out[[q]][nrow(out)])
      rows[[length(rows) + 1]] <- data.frame(
        scenario = nm, quantity = q, statistic = "sum", value = sum(out[[q]]))
    }
  }
  do.call(rbind, rows)
}

test_that("intervention output matches the committed reference values", {
  skip_if_not_installed("malariasimulation")

  got <- reference_int_summary()

  if (nzchar(Sys.getenv("FLEET_REGENERATE_REFERENCE"))) {
    write.csv(transform(got, value = vapply(got$value, format, character(1), digits = 17)),
              test_path("reference-interventions.csv"), row.names = FALSE, quote = FALSE)
    message("FLEET_REGENERATE_REFERENCE: rewrote reference-interventions.csv")
  }

  ref <- utils::read.csv(test_path("reference-interventions.csv"))
  key <- function(d) paste(d$scenario, d$quantity, d$statistic)
  expect_setequal(key(ref), key(got))
  ref <- ref[match(key(got), key(ref)), ]

  for (i in seq_len(nrow(got))) {
    expect_equal(got$value[i], ref$value[i], tolerance = 1e-6,
                 label = sprintf("%s / %s [%s]", got$scenario[i], got$quantity[i],
                                 got$statistic[i]))
  }
})

test_that("no intervention acts before the day it is deployed", {
  skip_if_not_installed("malariasimulation")

  ## The onset-timing invariant stated directly, independent of any committed
  ## number: up to and including the deployment day, every scenario must still be
  ## sitting on the no-intervention trajectory, in every column -- the round lands
  ## at the end of its day. This is what caught the 10-day lead, and it stays
  ## meaningful even if the CSV is regenerated for an unrelated reason.
  base <- run_simulation_ode(REF_INT_DAYS, eqm(gp_bands(), 20))
  num <- setdiff(names(base), "timestep")
  for (nm in names(reference_int_scenarios())) {
    out <- run_simulation_ode(REF_INT_DAYS, eqm(reference_int_scenarios()[[nm]], 20))
    i <- which(out$timestep <= NET_DAY)
    expect_equal(as.matrix(out[i, num]), as.matrix(base[i, num]), tolerance = 1e-9,
                 label = sprintf("%s up to the deployment day", nm))
  }
  ## and the nets act the day after: the biting rate falls at once, so the force
  ## of infection on mosquitoes moves on row NET_DAY + 1
  nets <- run_simulation_ode(REF_INT_DAYS, eqm(reference_int_scenarios()$nets, 20))
  expect_false(isTRUE(all.equal(nets$FOIM[NET_DAY + 1], base$FOIM[NET_DAY + 1])))
  ## nor late: each probe's first effect on people lands on the onset day it is
  ## sampled on
  for (nm in names(reference_int_scenarios())) {
    out <- run_simulation_ode(REF_INT_DAYS, eqm(reference_int_scenarios()[[nm]], 20))
    q <- out$n_inc_clinical_730_3650; b <- base$n_inc_clinical_730_3650
    expect_equal(which(abs(q - b) > 1e-9 * abs(b))[1], REF_INT_DAYS_SAMPLED[[nm]][["onset"]],
                 label = sprintf("%s: first day clinical incidence moves", nm))
  }
})


## ---- P. vivax reference -----------------------------------------------------
#
# The same net for the vivax block, which neither CSV above can see: a vivax
# list at low and moderate transmission, untreated, and two moving to radical
# cure mid-run, by primaquine (four liver-stage levels) and by tafenoquine
# (eleven), pinned the same way (final day and whole-run sum, to 1e-6), in a CSV
# of its own. Same regeneration switch.

REF_PV_DAYS <- 730

reference_pv_scenarios <- function() {
  pv <- function() {
    malariasimulation::get_parameters(list(
      human_population = 10000,
      clinical_incidence_rendering_min_ages = c(0, 0),
      clinical_incidence_rendering_max_ages = c(1824, 36499),
      age_group_rendering_min_ages = 0, age_group_rendering_max_ages = 36499),
      parasite = "vivax")
  }
  to_radical_cure <- function(rc_drug) {
    rc <- pv()
    rc <- malariasimulation::set_drugs(rc, list(malariasimulation::CQ_params_vivax, rc_drug))
    rc <- malariasimulation::set_clinical_treatment(rc, drug = 1, timesteps = c(1, 365),
                                                    coverages = c(0.3, 0))
    malariasimulation::set_clinical_treatment(rc, drug = 2, timesteps = c(1, 365),
                                              coverages = c(0, 0.5))
  }
  list(eir1 = eqm(pv(), 1), eir10 = eqm(pv(), 10),
       radical_cure = eqm(to_radical_cure(malariasimulation::CQ_PQ_params_vivax), 10),
       radical_cure_tq = eqm(to_radical_cure(malariasimulation::CQ_TQ_params_vivax), 10))
}

REF_PV_QUANTITIES <- c(
  "n_detect_lm_730_3650", "n_detect_pcr_730_3650",
  "n_inc_clinical_0_1824", "n_inc_clinical_0_36499",
  "n_infections", "n_relapses", "n_with_hypnozoites", "FOIM", "EIR",
  "S_count", "D_count", "A_count", "U_count", "Tr_count",
  "iaa_mean", "ica_mean", "hypnozoites_mean", "n_age_0_36499")

reference_pv_summary <- function() {
  rows <- list()
  scen <- reference_pv_scenarios()
  for (nm in names(scen)) {
    out <- run_simulation_ode(REF_PV_DAYS, scen[[nm]])
    expect_true(all(REF_PV_QUANTITIES %in% names(out)))
    for (q in REF_PV_QUANTITIES) {
      v <- out[[q]]
      rows[[length(rows) + 1]] <- data.frame(scenario = nm, quantity = q,
                                             statistic = "final", value = v[length(v)])
      rows[[length(rows) + 1]] <- data.frame(scenario = nm, quantity = q,
                                             statistic = "sum", value = sum(v))
    }
  }
  do.call(rbind, rows)
}

test_that("vivax output matches the committed reference values", {
  skip_if_not_installed("malariasimulation")

  got <- reference_pv_summary()

  if (nzchar(Sys.getenv("FLEET_REGENERATE_REFERENCE"))) {
    write.csv(transform(got, value = vapply(got$value, format, character(1), digits = 17)),
              test_path("reference-vivax.csv"), row.names = FALSE, quote = FALSE)
    message("FLEET_REGENERATE_REFERENCE: rewrote reference-vivax.csv")
  }

  ref <- utils::read.csv(test_path("reference-vivax.csv"))
  key <- function(d) paste(d$scenario, d$quantity, d$statistic)
  expect_setequal(key(ref), key(got))
  ref <- ref[match(key(got), key(ref)), ]

  for (i in seq_len(nrow(got))) {
    expect_equal(got$value[i], ref$value[i], tolerance = 1e-6,
                 label = sprintf("%s / %s [%s]", got$scenario[i], got$quantity[i],
                                 got$statistic[i]))
  }
})

## ---- P. vivax intervention reference ------------------------------------------
#
# The vivax counterpart of reference-interventions.csv, in a CSV of its own: the
# same three probes -- one net campaign, TBV on ages that straddle age groups, a
# mass PEV campaign split across two bands -- on a vivax list at EIR 10, sampled
# on the deployment day and on the first day each effect reaches the sampled
# quantities. The days are vivax's own, because its delay lines are: humans feel
# the EIR de = 10 days on, and mosquitoes feel human infectivity the same day
# (delay_gam = 0). Nets cut the biting rate from day 101, which humans feel on
# row 111; TBV's fewer infected mosquitoes leave incubation 9 days after day 101
# and bite 10 days after that (row 120); mass PEV protects from its last primary
# dose (row 191), and cuts relapses as well as bites.

REF_PV_INT_DAYS_SAMPLED <- list(
  nets = c(before = NET_DAY, onset = NET_DAY + 11, after = NET_DAY + 16),
  tbv_offgrid = c(before = NET_DAY, onset = NET_DAY + 20, after = NET_DAY + 50),
  mass_pev_split = c(before = NET_DAY, onset = NET_DAY + 91, after = NET_DAY + 136))

pv_int_bands <- function() {
  malariasimulation::get_parameters(list(
    human_population = 10000,
    clinical_incidence_rendering_min_ages = 730,
    clinical_incidence_rendering_max_ages = 3650), parasite = "vivax")
}

reference_pv_int_scenarios <- function() {
  list(
    nets = malariasimulation::set_bednets(
      pv_int_bands(), timesteps = NET_DAY, coverages = 0.8, retention = 5 * 365,
      dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
      rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365),
    tbv_offgrid = malariasimulation::set_tbv(
      pv_int_bands(), timesteps = NET_DAY, coverages = 0.9, ages = c(18, 19, 20)),
    mass_pev_split = malariasimulation::set_mass_pev(
      pv_int_bands(), profile = malariasimulation::rtss_profile,
      timesteps = NET_DAY, coverages = 1,
      min_ages = c(0, 1195), max_ages = c(1195, 1825), min_wait = 0,
      booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
      booster_profile = list(malariasimulation::rtss_booster_profile)))
}

REF_PV_INT_QUANTITIES <- c("EIR", "n_detect_lm_730_3650", "n_inc_clinical_730_3650",
                           "n_relapses")

reference_pv_int_summary <- function() {
  rows <- list()
  scen <- reference_pv_int_scenarios()
  for (nm in names(scen)) {
    out <- run_simulation_ode(REF_INT_DAYS, eqm(scen[[nm]], 10))
    expect_true(all(REF_PV_INT_QUANTITIES %in% names(out)))
    for (q in REF_PV_INT_QUANTITIES) {
      days <- REF_PV_INT_DAYS_SAMPLED[[nm]]
      for (s in names(days)) {
        rows[[length(rows) + 1]] <- data.frame(
          scenario = nm, quantity = q, statistic = s,
          value = out[[q]][out$timestep == days[[s]]])
      }
      rows[[length(rows) + 1]] <- data.frame(
        scenario = nm, quantity = q, statistic = "final", value = out[[q]][nrow(out)])
      rows[[length(rows) + 1]] <- data.frame(
        scenario = nm, quantity = q, statistic = "sum", value = sum(out[[q]]))
    }
  }
  do.call(rbind, rows)
}

test_that("vivax intervention output matches the committed reference values", {
  skip_if_not_installed("malariasimulation")

  got <- reference_pv_int_summary()

  if (nzchar(Sys.getenv("FLEET_REGENERATE_REFERENCE"))) {
    write.csv(transform(got, value = vapply(got$value, format, character(1), digits = 17)),
              test_path("reference-vivax-interventions.csv"), row.names = FALSE, quote = FALSE)
    message("FLEET_REGENERATE_REFERENCE: rewrote reference-vivax-interventions.csv")
  }

  ref <- utils::read.csv(test_path("reference-vivax-interventions.csv"))
  key <- function(d) paste(d$scenario, d$quantity, d$statistic)
  expect_setequal(key(ref), key(got))
  ref <- ref[match(key(got), key(ref)), ]

  for (i in seq_len(nrow(got))) {
    expect_equal(got$value[i], ref$value[i], tolerance = 1e-6,
                 label = sprintf("%s / %s [%s]", got$scenario[i], got$quantity[i],
                                 got$statistic[i]))
  }
})

test_that("no vivax intervention acts before the day it is deployed, nor late", {
  skip_if_not_installed("malariasimulation")

  ## As for falciparum: up to and including the deployment day every probe sits
  ## on the no-intervention trajectory in every column. And each probe's first
  ## effect on people -- the day's clinical cases and relapses -- lands on its
  ## onset day, so the days above are the onsets they are said to be. (The EIR
  ## follows PEV 20 days later, once fewer infections have reached mosquitoes and
  ## come back.)
  base <- run_simulation_ode(REF_INT_DAYS, eqm(pv_int_bands(), 10))
  num <- setdiff(names(base), "timestep")
  scen <- reference_pv_int_scenarios()
  for (nm in names(scen)) {
    out <- run_simulation_ode(REF_INT_DAYS, eqm(scen[[nm]], 10))
    i <- which(out$timestep <= NET_DAY)
    expect_equal(as.matrix(out[i, num]), as.matrix(base[i, num]), tolerance = 1e-9,
                 label = sprintf("vivax %s up to the deployment day", nm))
    onset <- REF_PV_INT_DAYS_SAMPLED[[nm]][["onset"]]
    for (q in c("n_inc_clinical_730_3650", "n_relapses")) {
      moved <- which(abs(out[[q]] - base[[q]]) > 1e-9 * abs(base[[q]]))
      expect_equal(moved[1], onset, label = sprintf("vivax %s: first day %s moves", nm, q))
    }
  }
})
