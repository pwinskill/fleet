# Scenario definitions and the shared per-run summariser.
#
# Sourced by BOTH run_replicates.R (which runs them through both models and
# writes the CSVs) and check_drift.R (which re-runs fleet only and compares
# against the committed reference). They live here so a drift check cannot
# silently test a different set of scenarios from the one the reference was
# built on -- and so sourcing the scenarios does not start a 25-minute run.
#
# Expects ROOT and theme.R's constants to be in scope already.

if (!exists("AGE_EDGES")) stop("source comparison/constants.R before this file.")
if (!exists("SMOKE")) SMOKE <- nzchar(Sys.getenv("CMP_SMOKE"))

## ---- shared parameter scaffolding -------------------------------------------
band_lo <- head(AGE_EDGES, -1) * 365; band_hi <- tail(AGE_EDGES, -1) * 365
tag_of  <- function(lo, hi) paste0(round(lo), "_", round(hi))
AGE_TAGS <- tag_of(band_lo, band_hi)

## rendering ranges: every scenario gets the three defaults (LM prevalence 2-10y,
## clinical and severe incidence 0-5y); the reference scenario adds the
## age-profile bands to all three families.
set_bands <- function(p, age_profile = FALSE) {
  add <- function(lo, hi) if (age_profile) list(c(band_lo, lo), c(band_hi, hi)) else list(lo, hi)
  r <- add(730, 3650)
  p$prevalence_rendering_min_ages <- r[[1]]; p$prevalence_rendering_max_ages <- r[[2]]
  ## both incidence families over under-5 AND all ages: the impact figure reports
  ## the burden a programme actually carries as well as the young-child burden
  r <- add(c(0, 0), c(1825, 36500))
  p$clinical_incidence_rendering_min_ages <- r[[1]]; p$clinical_incidence_rendering_max_ages <- r[[2]]
  p$severe_incidence_rendering_min_ages   <- r[[1]]; p$severe_incidence_rendering_max_ages   <- r[[2]]
  p
}
base_params <- function(seasonal = FALSE, age_profile = FALSE) {
  ov <- list(human_population = POP)
  if (seasonal) ov <- c(ov, list(model_seasonality = TRUE), SEASON)
  set_bands(get_parameters(ov), age_profile)
}

## ---- scenarios ---------------------------------------------------------------
## Each returns list(p = parameters, eir = init_EIR, years = horizon).
Y_INT <- BURN_Y * 365                       # intervention start (day)
scenarios <- list()

for (E in EIR_GRID) scenarios[[paste0("eir_", E)]] <- list(
  p = set_equilibrium(base_params(age_profile = (E == EIR_REF)), init_EIR = E),
  eir = E, years = BURN_Y + 3L)

scenarios$seasonal <- list(
  p = set_equilibrium(base_params(seasonal = TRUE), init_EIR = EIR_REF),
  eir = EIR_REF, years = BURN_Y + 3L)

scenarios$nets <- local({
  p <- set_bednets(base_params(), timesteps = Y_INT, coverages = 0.8,
                   retention = 5 * 365, dn0 = matrix(0.387), rn = matrix(0.563),
                   rnm = matrix(0.24), gamman = 2.64 * 365)
  list(p = set_equilibrium(p, init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + 6L)
})

scenarios$irs <- local({
  rounds <- Y_INT + c(0, 1, 2) * 365
  m <- function(v) matrix(v, nrow = length(rounds), ncol = 1)
  p <- set_spraying(base_params(), timesteps = rounds, coverages = rep(0.8, 3),
                    ls_theta = m(2.025), ls_gamma = m(-0.009),
                    ks_theta = m(-2.222), ks_gamma = m(0.008),
                    ms_theta = m(-1.232), ms_gamma = m(-0.009))
  list(p = set_equilibrium(p, init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + 6L)
})

scenarios$smc <- local({
  p <- base_params(seasonal = TRUE)
  p <- set_drugs(p, list(SP_AQ_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.45)
  rounds <- as.vector(sapply(0:2, function(y) Y_INT + y * 365 + c(0, 30, 60, 90) + 200))
  p <- set_smc(p, drug = 1, timesteps = rounds, coverages = rep(0.9, length(rounds)),
               min_ages = rep(round(0.25 * 365), length(rounds)),
               max_ages = rep(round(5 * 365), length(rounds)))
  list(p = set_equilibrium(p, init_EIR = 15), eir = 15, years = BURN_Y + 3L)
})

scenarios$pev <- local({
  p <- set_pev_epi(base_params(), profile = rtss_profile, timesteps = Y_INT,
                   coverages = 0.9, min_wait = 0, age = 5 * 30,
                   booster_spacing = 12 * 30, booster_coverage = matrix(0.8),
                   booster_profile = list(rtss_booster_profile))
  list(p = set_equilibrium(p, init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + 6L)
})

scenarios$treatment <- local({
  p <- set_drugs(base_params(), list(AL_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = c(1, Y_INT), coverages = c(0.2, 0.6))
  list(p = set_equilibrium(p, init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + 6L)
})

## custom demography: high infant and elderly mortality, so the equilibrium age
## structure departs strongly from the default exponential one
scenarios$demography <- local({
  dr <- c(0.048, 0.007, 0.003, 0.004, 0.008, 0.020, 0.050, 0.120) / 365
  ag <- round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365)
  p <- set_demography(base_params(age_profile = TRUE), agegroups = ag, timesteps = 0,
                      deathrates = matrix(dr, nrow = 1))
  list(p = set_equilibrium(p, init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + 3L)
})

## ---- ts_*: long-horizon programme scenarios ---------------------------------
## The scenarios above each isolate ONE builder over 6 years, which is the right
## shape for attributing a difference but not for showing what a programme looks
## like. These five run 15 years past deployment in a seasonal setting so the
## repeated-campaign dynamics are visible -- five net distributions decaying and
## being replaced, SMC pulsing four times a year for fifteen years -- and they
## share a common no-intervention reference so the three tiers (nothing, one
## thing, everything) can be read against each other. All at EIR 20 seasonal, so
## every row of the figure is the same setting with more added to it.
TS_Y <- 15L
ts_base <- function() set_bands(get_parameters(c(list(human_population = POP),
                                                 list(model_seasonality = TRUE), SEASON)))
## nets every 3 years: 5 campaigns over the 15-year window
ts_net_rounds <- Y_INT + seq(0, by = 3 * 365, length.out = 5L)
ts_nets_on <- function(p) {
  n <- length(ts_net_rounds)
  set_bednets(p, timesteps = ts_net_rounds, coverages = rep(0.8, n),
              retention = 5 * 365, dn0 = matrix(rep(0.387, n)), rn = matrix(rep(0.563, n)),
              rnm = matrix(rep(0.24, n)), gamman = rep(2.64 * 365, n))
}
## SMC: 4 monthly rounds a year, every year of the window, peak-season aligned
ts_smc_rounds <- as.vector(sapply(seq_len(TS_Y) - 1L,
                                  function(y) Y_INT + y * 365 + c(0, 30, 60, 90) + 200))
ts_smc_on <- function(p) {
  n <- length(ts_smc_rounds)
  set_smc(p, drug = 1, timesteps = ts_smc_rounds, coverages = rep(0.9, n),
          min_ages = rep(round(0.25 * 365), n), max_ages = rep(round(5 * 365), n))
}
## case management: SP-AQ throughout (SMC needs a drug), scaled up at deployment
ts_treat_on <- function(p, hi = 0.6) set_clinical_treatment(
  set_drugs(p, list(SP_AQ_params)), drug = 1, timesteps = c(1, Y_INT), coverages = c(0.2, hi))
ts_drug_only <- function(p) set_clinical_treatment(
  set_drugs(p, list(SP_AQ_params)), drug = 1, timesteps = 1, coverages = 0.2)

ts_scen <- list(
  ts_none  = function() ts_drug_only(ts_base()),
  ts_nets  = function() ts_nets_on(ts_drug_only(ts_base())),
  ts_smc   = function() ts_smc_on(ts_drug_only(ts_base())),
  ts_treat = function() ts_treat_on(ts_base()),
  ts_all   = function() ts_smc_on(ts_nets_on(ts_treat_on(ts_base()))))
for (nm in names(ts_scen)) scenarios[[nm]] <- local({
  f <- ts_scen[[nm]]
  list(p = set_equilibrium(f(), init_EIR = EIR_REF), eir = EIR_REF, years = BURN_Y + TS_Y)
})

## CMP_ONLY=a,b -> run just those scenarios and merge their rows into the existing CSVs
ONLY <- Filter(nzchar, strsplit(Sys.getenv("CMP_ONLY"), ",")[[1]])
## the digest must describe the WHOLE scenario set, not whichever subset this
## invocation happens to run -- otherwise every CMP_ONLY run reports the
## reference as stale when nothing has changed
SCENARIOS_ALL <- scenarios
if (length(ONLY)) { stopifnot(all(ONLY %in% names(scenarios))); scenarios <- scenarios[ONLY] }

if (SMOKE) scenarios <- lapply(scenarios, function(s) { s$years <- 4L; s })

## ---- one summariser for BOTH models ------------------------------------------
## Takes a wide daily output table and returns compact tidy pieces. The ODE table
## carries a day-0 seed row the IBM lacks; callers drop it first so both bin on
## days 1..N and monthly bins align. The age profile is only returned when the
## table carries the age-profile bands (the reference scenario).
summarise_run <- function(df, years, tags = AGE_TAGS) {
  n <- nrow(df); day <- seq_len(n)
  pooled <- function(num, den, rows) sum(num[rows]) / sum(den[rows])
  rate   <- function(kind, tag, rows, per = 365)
    pooled(df[[paste0("n_inc_", kind, "_", tag)]], df[[paste0("n_age_", tag)]], rows) * per
  prev   <- function(tag, rows) pooled(df[[paste0("n_detect_lm_", tag)]], df[[paste0("n_age_", tag)]], rows)
  obs <- (n - 3 * 365 + 1):n                    # final three years
  ## equilibrium quantities over the observation window
  eq <- data.frame(
    pfpr_2_10 = prev("730_3650", obs), clin_0_5 = rate("clinical", "0_1825", obs),
    sev_0_5 = rate("severe", "0_1825", obs, 365 * 1000),
    clin_all = rate("clinical", "0_36500", obs),
    sev_all = rate("severe", "0_36500", obs, 365 * 1000))
  ## age profiles over the observation window (reference scenario only)
  age <- if (all(paste0("n_detect_lm_", tags) %in% names(df))) {
    n_band <- vapply(tags, function(tg) sum(df[[paste0("n_age_", tg)]][obs]), numeric(1))
    data.frame(
      age_lo = band_lo / 365, age_hi = band_hi / 365, age_mid = (band_lo + band_hi) / 2 / 365,
      pop_frac = n_band / sum(n_band),            # share of the 0-85 population in each band
      prev = vapply(tags, prev, numeric(1), rows = obs),
      clin = vapply(tags, rate, numeric(1), kind = "clinical", rows = obs),
      sev  = vapply(tags, rate, numeric(1), kind = "severe", rows = obs, per = 365 * 1000))
  }
  ## monthly series (30-day bins), pooled within each bin; `year` = bin start
  mon <- (day - 1) %/% 30
  mo <- function(f) as.numeric(tapply(day, mon, f))
  monthly <- data.frame(
    year = as.numeric(names(tapply(day, mon, length))) * 30 / 365,
    pfpr_2_10 = mo(function(ix) prev("730_3650", ix)),
    clin_0_5  = mo(function(ix) rate("clinical", "0_1825", ix)),
    sev_0_5   = mo(function(ix) rate("severe", "0_1825", ix, 365 * 1000)),
    clin_all  = mo(function(ix) rate("clinical", "0_36500", ix)),
    sev_all   = mo(function(ix) rate("severe", "0_36500", ix, 365 * 1000)))
  ## final-year day-of-year series (for the seasonal cycle), 7-day pooled bins
  fy <- (n - 365 + 1):n; wk <- (seq_along(fy) - 1) %/% 7
  doy <- data.frame(
    doy = as.numeric(tapply(seq_along(fy), wk, function(ix) mean(ix))),
    pfpr_2_10 = as.numeric(tapply(seq_along(fy), wk, function(ix) prev("730_3650", fy[ix]))),
    clin_0_5  = as.numeric(tapply(seq_along(fy), wk, function(ix) rate("clinical", "0_1825", fy[ix]))))
  ## realised EIR (IBM emits EIR_<species> as total bites; fleet emits per-adult-per-year)
  eir <- if (any(grepl("^EIR_", names(df)))) {
    ec <- grep("^EIR_", names(df), value = TRUE)
    sum(rowSums(df[obs, ec, drop = FALSE])) / length(obs) / POP * 365
  } else mean(df$EIR[obs])
  eq$eir_realised <- eir
  list(eq = eq, age = age, monthly = monthly, doy = doy)
}
tag_parts <- function(r, nm, model, rep)
  lapply(r, function(x) if (is.null(x)) NULL else cbind(scenario = nm, model = model, rep = rep, x))

## ---- provenance ---------------------------------------------------------------
## A hash of what the IBM was actually run on. The committed IBM rows stay valid
## for any fleet-side change -- the IBM does not depend on fleet -- but they go
## stale the moment a scenario definition or malariasimulation itself changes.
## Without this there is no way to notice that, and a drift check would keep
## comparing against a reference built for different scenarios.
##
## Hashes the fields that actually reach the IBM: every parameter that differs
## from a bare get_parameters(), plus the EIR and horizon. Deliberately not the
## whole list, which carries defaults that shift between malariasimulation
## versions for reasons unrelated to these scenarios.
scenario_digest <- function(scen = SCENARIOS_ALL) {
  base <- get_parameters()
  canon <- function(s) {
    changed <- Filter(function(k) !identical(s$p[[k]], base[[k]]), names(s$p))
    list(eir = s$eir, years = s$years,
         p = lapply(sort(changed), function(k) list(k, s$p[[k]])))
  }
  digest::digest(lapply(scen[sort(names(scen))], canon), algo = "xxhash64")
}

ibm_reference <- function() list(
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  malariasimulation = as.character(utils::packageVersion("malariasimulation")),
  R = paste0(R.version$major, ".", R.version$minor),
  n_rep = N_REP, population = POP, burn_in_years = BURN_Y,
  scenarios = sort(names(SCENARIOS_ALL)), scenario_digest = scenario_digest())

## ---- running fleet -------------------------------------------------------------
## Shared so a drift check cannot accidentally run fleet differently from the way
## the committed reference was produced: same tuning, same horizons, same
## summariser, same day-0 row dropped.
run_fleet <- function(scen = scenarios) {
  suppressMessages(library(fleet))
  log_msg("fleet: %d scenarios", length(scen))
  out <- lapply(names(scen), function(nm) {
    s <- scen[[nm]]
    el <- system.time(
      ## s$p already carries init_EIR: every scenario is built through
      ## set_equilibrium(), which is where fleet reads the target EIR from now.
      o <- fleet::run_simulation_ode(timesteps = s$years * 365, parameters = s$p,
                                     tuning = list(rtol = 1e-6, step_size_max = 10))
    )[["elapsed"]]
    r <- summarise_run(o[-1, ], s$years)          # drop the day-0 seed row
    r$timing <- data.frame(years = s$years, elapsed_s = el)
    tag_parts(r, nm, "fleet", 0L)
  })
  names(out) <- names(scen)
  log_msg("fleet done: %.1f s total", sum(vapply(out, function(o) o$timing$elapsed_s, numeric(1))))
  out
}
