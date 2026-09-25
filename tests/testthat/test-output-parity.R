# Raw-output parity with malariasimulation.
#
# fleet returns malariasimulation's daily table: for a parameter list, the
# columns the IBM would render for it, under the same names and with the same
# meanings. These tests run the IBM itself beside fleet across rendering
# configurations chosen to exercise the band rules -- the defaults, a site-style
# partition, bands that share an edge, fractional and very large edges, single
# days, overlapping families, the immunity-mean bands, a treated run, and the
# vivax families and treatments -- and check the column sets exactly,
# fleet's own identities exactly, and the values against a larger IBM run
# loosely: tight enough to catch a column that means something else, not a test
# of level agreement, which is fleetcheck's job.

## the age-banded families, the part of the table the rendering lists control
BAND_FAMILIES <- paste0("^(n_age|n_detect_lm|p_detect_lm|n_detect_pcr|n_inc|p_inc|",
                        "n_inc_clinical|p_inc_clinical|n_inc_severe|p_inc_severe)_")

## malariasimulation columns fleet does not render, and fleet's three of its own
IBM_ONLY <- c("n_bitten", "infectivity", "EIR_gamb", "FOIM_gamb", "mu_gamb",
              "E_gamb_count", "L_gamb_count", "P_gamb_count", "Sm_gamb_count",
              "Pm_gamb_count", "Im_gamb_count", "total_M_gamb", "natural_deaths")
IBM_ONLY_TREATED <- c("n_treated", "n_drug_efficacy_failures", "n_successfully_treated")
FLEET_ONLY <- c("EIR", "FOIM", "Ph_count")

render_configs <- function() list(
  defaults = list(),
  site_partition = list(
    prevalence_rendering_min_ages = c(2, 0) * 365,
    prevalence_rendering_max_ages = c(10, 100) * 365 - 1,
    clinical_incidence_rendering_min_ages = c(0, 5, 15) * 365,
    clinical_incidence_rendering_max_ages = c(5, 15, 100) * 365 - 1,
    severe_incidence_rendering_min_ages = c(0, 5) * 365,
    severe_incidence_rendering_max_ages = c(5, 100) * 365 - 1,
    incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 100 * 365 - 1,
    age_group_rendering_min_ages = c(0, 1, 2) * 365,
    age_group_rendering_max_ages = c(1, 2, 3) * 365 - 1),
  shared_edges = list(
    clinical_incidence_rendering_min_ages = c(0, 1825),
    clinical_incidence_rendering_max_ages = c(1825, 5475)),
  fractional_and_large = list(
    prevalence_rendering_min_ages = 0.5 * 365, prevalence_rendering_max_ages = 5 * 365 + 0.5,
    clinical_incidence_rendering_min_ages = 1e5, clinical_incidence_rendering_max_ages = 2e5),
  single_days = list(
    age_group_rendering_min_ages = c(0, 1, 364, 365),
    age_group_rendering_max_ages = c(0, 1, 364, 729)),
  immunity_bands = list(
    ica_rendering_min_ages = c(0, 1825), ica_rendering_max_ages = c(1824, 36499),
    ib_rendering_min_ages = 0, ib_rendering_max_ages = 1824,
    id_rendering_min_ages = c(730, 1825), id_rendering_max_ages = c(3650, 36499),
    icm_rendering_min_ages = 0, icm_rendering_max_ages = 364,
    iva_rendering_min_ages = 0, iva_rendering_max_ages = 36499,
    ivm_rendering_min_ages = 0, ivm_rendering_max_ages = 364),
  overlapping = list(
    prevalence_rendering_min_ages = c(0, 730), prevalence_rendering_max_ages = c(1825, 3650),
    clinical_incidence_rendering_min_ages = c(0, 0, 1825),
    clinical_incidence_rendering_max_ages = c(36500, 1825, 5475),
    severe_incidence_rendering_min_ages = c(0, 0), severe_incidence_rendering_max_ages = c(36500, 1825)))

test_that("fleet renders the IBM's columns, and only those, for every rendering configuration", {
  skip_if_not_installed("malariasimulation")
  cfgs <- render_configs()
  cfgs$treated <- cfgs$site_partition
  for (nm in names(cfgs)) {
    p <- malariasimulation::get_parameters(c(list(human_population = 200), cfgs[[nm]]))
    if (nm == "treated") {
      p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
      p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
    }
    p <- eqm(p, 20)
    ibm <- malariasimulation::run_simulation(3, p)
    fl <- suppressWarnings(run_simulation_ode(3, p))
    expect_equal(fl$timestep, ibm$timestep, label = nm)
    expect_setequal(grep(BAND_FAMILIES, names(fl), value = TRUE),
                    grep(BAND_FAMILIES, names(ibm), value = TRUE))
    expect_setequal(setdiff(names(ibm), names(fl)),
                    c(IBM_ONLY, if (nm == "treated") IBM_ONLY_TREATED))
    expect_setequal(setdiff(names(fl), names(ibm)), FLEET_ONLY)
  }
})

## The vivax families: relapses and hypnozoite carriage over their own lists,
## the vivax immunity means, and the shared families on a vivax list.
render_configs_pv <- function() list(
  defaults = list(),
  families = list(
    prevalence_rendering_min_ages = c(0, 730), prevalence_rendering_max_ages = c(729, 3650),
    clinical_incidence_rendering_min_ages = c(0, 1825),
    clinical_incidence_rendering_max_ages = c(1824, 36499),
    incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 36499,
    age_group_rendering_min_ages = c(0, 365), age_group_rendering_max_ages = c(364, 1824),
    incidence_relapse_rendering_min_ages = c(0, 1825),
    incidence_relapse_rendering_max_ages = c(1824, 36499),
    n_with_hypnozoites_rendering_min_ages = c(0, 1825),
    n_with_hypnozoites_rendering_max_ages = c(1824, 36499),
    iaa_rendering_min_ages = c(0, 1825), iaa_rendering_max_ages = c(1824, 36499),
    iam_rendering_min_ages = 0, iam_rendering_max_ages = 364,
    ica_rendering_min_ages = 1825, ica_rendering_max_ages = 36499,
    hypnozoites_rendering_min_ages = c(0, 730), hypnozoites_rendering_max_ages = c(729, 36499)),
  severe_band = list(
    severe_incidence_rendering_min_ages = 0, severe_incidence_rendering_max_ages = 36499))

test_that("fleet renders the IBM's vivax columns, and only those", {
  skip_if_not_installed("malariasimulation")
  cfgs <- render_configs_pv()
  cfgs$treated <- cfgs$families
  cfgs$radical_cure <- cfgs$families
  for (nm in names(cfgs)) {
    # large enough that every relapse band sees a relapse on some day: the IBM
    # renders n_inc_relapse_* only on such days, so a band that never sees one
    # would have no column
    # set on the list rather than passed to get_parameters(), which does not
    # know incidence_relapse_rendering_*: the IBM reads it all the same
    p <- malariasimulation::get_parameters(list(human_population = 5000), parasite = "vivax")
    p[names(cfgs[[nm]])] <- cfgs[[nm]]
    if (nm %in% c("treated", "radical_cure")) {
      drug <- if (nm == "treated") malariasimulation::CQ_params_vivax else
        malariasimulation::CQ_PQ_params_vivax
      p <- malariasimulation::set_drugs(p, list(drug))
      p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
    }
    p <- eqm(p, 10)
    ibm <- malariasimulation::run_simulation(5, p)
    fl <- run_simulation_ode(5, p)
    expect_equal(fl$timestep, ibm$timestep, label = nm)
    expect_setequal(setdiff(names(ibm), names(fl)),
                    c(IBM_ONLY, if (nm %in% c("treated", "radical_cure")) IBM_ONLY_TREATED))
    expect_setequal(setdiff(names(fl), names(ibm)), FLEET_ONLY)
  }
})

test_that("each p_* column is the expected count behind its n_* partner", {
  skip_if_not_installed("malariasimulation")
  o <- run_simulation_ode(60, eqm(gp_bands(), 20))
  for (tag in c("730_3650")) {
    # the IBM samples n_detect_lm_* and sums the probabilities into p_detect_lm_*;
    # a mean field's expected count is both, and prevalence is the count over n_age
    expect_equal(o[[paste0("p_detect_lm_", tag)]], o[[paste0("n_detect_lm_", tag)]])
    expect_true(all(pfpr(o, tag) > 0 & pfpr(o, tag) < 1))
    expect_gt(min(o[[paste0("p_detect_lm_", tag)]]), 1)     # a count, not a proportion
  }
  for (tag in c("0_1825", "730_3650", "0_36500")) {
    expect_equal(o[[paste0("p_inc_clinical_", tag)]], o[[paste0("n_inc_clinical_", tag)]])
  }
  for (tag in c("730_3650", "0_36500")) {
    expect_equal(o[[paste0("p_inc_severe_", tag)]], o[[paste0("n_inc_severe_", tag)]])
    # malariasimulation 3.0.0 renders p_inc_* from a field no list sets: always 0
    expect_true(all(o[[paste0("p_inc_", tag)]] == 0))
  }
  expect_equal(o$n_infections, o$n_inc_0_36500)
})

test_that("a band holds the whole-day ages lower to upper, both ends included", {
  skip_if_not_installed("malariasimulation")
  bands <- list(
    c(0, 364), c(365, 729), c(730, 36499),                   # tile every age
    c(0, 1824), c(1825, 5474),                               # tile 0..5474
    c(0, 1825), c(1825, 5475),                               # share day 1825
    c(1825, 1825), c(5475, 5475),                            # single days
    c(0, 1000), c(0, 999), c(1000, 1000),                    # a day as a difference
    c(182.5, 1825.5), c(183, 1825),                          # fractional edges
    c(0, 0))                                                 # nobody is under a day old
  p <- malariasimulation::get_parameters(list(
    age_group_rendering_min_ages = vapply(bands, `[`, 0, 1),
    age_group_rendering_max_ages = vapply(bands, `[`, 0, 2)))
  expect_warning(o <- run_simulation_ode(30, eqm(p, 20)), "0_0")
  n <- function(lo, hi) o[[make.names(paste0("n_age_", lo, "_", hi))]]
  # the last band reaches above the absorbing top group's lower edge (80 years
  # on the default grid), so it holds the whole group
  N <- pop_total(o)
  expect_equal(n(0, 364) + n(365, 729) + n(730, 36499), N, tolerance = 1e-10)
  # bands sharing an edge double-count its day, as the IBM's do
  expect_equal(n(0, 1825) + n(1825, 5475),
               n(0, 1824) + n(1825, 5474) + n(1825, 1825) + n(5475, 5475), tolerance = 1e-12)
  expect_equal(n(1000, 1000), n(0, 1000) - n(0, 999), tolerance = 1e-10)
  expect_true(all(n(1000, 1000) > 0))
  expect_equal(n(182.5, 1825.5), n(183, 1825))
  expect_true(all(n(0, 0) == 0))
})

test_that("population counts are the IBM's: the drug-protected are in S_count", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::set_drugs(gp_bands(), list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.6)
  p <- malariasimulation::set_smc(p, drug = 1, timesteps = 100, coverages = 0.8,
                                  min_ages = 91, max_ages = 1825)
  o <- run_simulation_ode(200, eqm(p, 20))
  expect_equal(pop_total(o), rep(p$human_population, nrow(o)), tolerance = 1e-10)
  expect_true(all(o$Ph_count <= o$S_count))
  expect_gt(max(o$Ph_count), 0)
  # ft appears with treatment, as the IBM renders it, and not without
  expect_true("ft" %in% names(o))
  expect_false("ft" %in% names(run_simulation_ode(5, eqm(gp_bands(), 20))))
})

test_that("the banded counts agree with the IBM's in meaning", {
  skip_if_not_installed("malariasimulation")
  skip_on_cran()
  p <- malariasimulation::get_parameters(list(
    human_population = 5e4,
    prevalence_rendering_min_ages = c(365, 730, 0), prevalence_rendering_max_ages = c(1824, 3650, 36499),
    clinical_incidence_rendering_min_ages = c(0, 1825, 5475),
    clinical_incidence_rendering_max_ages = c(1824, 5474, 36499),
    severe_incidence_rendering_min_ages = c(0, 1825), severe_incidence_rendering_max_ages = c(1825, 36500),
    incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 36499,
    age_group_rendering_min_ages = c(365, 7300), age_group_rendering_max_ages = c(729, 7664)))
  p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
  p <- eqm(p, 20)
  set.seed(1)
  ibm <- malariasimulation::run_simulation(30, p)
  fl <- run_simulation_ode(30, p)
  ## 30-day totals. The tolerances bound fleet against the IBM's expected counts
  ## (its p_* columns where it has them) closely enough that a column meaning
  ## something else -- a proportion, a band off by a year, a family rendered over
  ## the wrong list -- fails by a wide margin; severe allows for its small count
  ## and for the IBM's own severe incidence relaxing from the seed.
  tot <- function(d, col) sum(d[[col]])
  rel <- function(col, ibm_col = col) tot(fl, col) / tot(ibm, ibm_col) - 1
  for (col in grep("^n_age_", names(ibm), value = TRUE)) {
    expect_lt(abs(rel(col)), 0.03, label = col)
  }
  for (tag in c("365_1824", "730_3650", "0_36499")) {
    expect_lt(abs(rel(paste0("p_detect_lm_", tag))), 0.05, label = tag)
    expect_lt(abs(rel(paste0("n_detect_pcr_", tag))), 0.05, label = tag)
  }
  for (tag in c("0_1824", "1825_5474", "5475_36499")) {
    expect_lt(abs(rel(paste0("n_inc_clinical_", tag), paste0("p_inc_clinical_", tag))), 0.10,
              label = tag)
  }
  for (tag in c("0_1825", "1825_36500")) {
    expect_lt(abs(rel(paste0("n_inc_severe_", tag), paste0("p_inc_severe_", tag))), 0.35,
              label = tag)
  }
  expect_lt(abs(rel("n_infections")), 0.05)
  for (st in c("S_count", "A_count", "D_count", "U_count")) expect_lt(abs(rel(st)), 0.05, label = st)
  expect_lt(abs(rel("Tr_count")), 0.10)
  expect_equal(fl$ft, ibm$ft)
  ## the acquired immunity means start where the IBM's do, from the same
  ## equilibrium, and barely move in a month. Not the maternal ones: the IBM
  ## seeds falciparum from the equilibrium bin one step older than each person
  ## (initial_immunity() takes which.max(age < lower edges)), so its infants start
  ## with the maternal immunity of 0.1 years later, ~0.58 of their own. That
  ## passes once the seeded infants age out; its births take their mothers'.
  for (v in c("ica", "ib", "iva", "id"))
    expect_lt(abs(rel(paste0(v, "_mean"))), 0.03, label = v)
})

test_that("the vivax counts agree with the IBM's in meaning", {
  skip_if_not_installed("malariasimulation")
  skip_on_cran()
  p <- malariasimulation::get_parameters(list(
    human_population = 5e4,
    prevalence_rendering_min_ages = c(730, 0), prevalence_rendering_max_ages = c(3650, 36499),
    clinical_incidence_rendering_min_ages = c(0, 1825),
    clinical_incidence_rendering_max_ages = c(1824, 36499),
    incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 36499,
    n_with_hypnozoites_rendering_min_ages = c(0, 1825),
    n_with_hypnozoites_rendering_max_ages = c(1824, 36499)), parasite = "vivax")
  p$incidence_relapse_rendering_min_ages <- c(0, 1825)
  p$incidence_relapse_rendering_max_ages <- c(1824, 36499)
  p <- malariasimulation::set_drugs(p, list(malariasimulation::CQ_PQ_params_vivax))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
  p <- eqm(p, 10)
  set.seed(1)
  ibm <- malariasimulation::run_simulation(30, p)
  fl <- run_simulation_ode(30, p)
  ## 30-day totals, as for falciparum. The IBM renders n_inc_relapse_* only on
  ## days with a relapse in the band, NA otherwise: a day without one.
  tot <- function(d, col) sum(d[[col]], na.rm = TRUE)
  rel <- function(col, ibm_col = col) tot(fl, col) / tot(ibm, ibm_col) - 1
  for (col in grep("^n_age_", names(ibm), value = TRUE)) expect_lt(abs(rel(col)), 0.03, label = col)
  for (tag in c("730_3650", "0_36499")) {
    expect_lt(abs(rel(paste0("n_detect_lm_", tag))), 0.05, label = tag)
    expect_lt(abs(rel(paste0("n_detect_pcr_", tag))), 0.05, label = tag)
  }
  for (tag in c("0_1824", "1825_36499")) {
    expect_lt(abs(rel(paste0("n_inc_clinical_", tag), paste0("p_inc_clinical_", tag))), 0.10,
              label = tag)
    expect_lt(abs(rel(paste0("n_inc_relapse_", tag))), 0.10, label = tag)
    expect_lt(abs(rel(paste0("n_with_hypnozoites_", tag))), 0.05, label = tag)
  }
  for (col in c("n_infections", "n_relapses", "n_with_hypnozoites"))
    expect_lt(abs(rel(col)), 0.05, label = col)
  for (st in c("S_count", "A_count", "D_count", "U_count")) expect_lt(abs(rel(st)), 0.05, label = st)
  expect_lt(abs(rel("Tr_count")), 0.10)
  ## vivax seeds each person from their own bin (initial_state_vivax() compares
  ## with the upper edges), so the maternal means agree from the start too, to
  ## within the few mothers who set them
  for (v in c("ica", "iaa", "hypnozoites"))
    expect_lt(abs(rel(paste0(v, "_mean"))), 0.03, label = v)
  for (v in c("icm", "iam")) expect_lt(abs(rel(paste0(v, "_mean"))), 0.08, label = v)
  expect_equal(fl$ft, ibm$ft)
})
