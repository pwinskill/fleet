# Run every comparison scenario through BOTH models and save tidy CSVs.
#
#   Rscript comparison/run_replicates.R          # ~25 min on 10 workers (see cost note)
#   Rscript comparison/render_figures.R          # seconds: figures from the CSVs
#   CMP_ONLY=nets,smc Rscript comparison/run_replicates.R   # re-run a subset, merge into the CSVs
#
# The IBM (malariasimulation) is run N_REP times per scenario with different
# seeds, in parallel on a PSOCK cluster, and summarised per replicate; blink is
# run once (it is deterministic). Both get the SAME parameter list, with the
# rendering bands set once so their output columns line up exactly.
#
# Every scenario is burned in BURN_Y years in the IBM so it reaches its own
# stochastic steady state; blink is seeded at the malariaEquilibrium fixed point
# and integrated over the same horizon so both are on a common clock.
#
# Cost note: the IBM's per-band rendering dominates its run time at 10k people
# (~2.5 s per simulated year with the three default ranges, ~5 s with the 12-band
# age profile added to every family), so only the reference scenario carries the
# age-profile bands. Jobs are load-balanced, longest first.

## No absolute paths anywhere in here. BLINK_LIB prepends an R library, for
## installations that do not pick up R_LIBS_USER (the Windows-arm64 setup this was
## developed on); leave it unset and your normal library is used. ROOT is found by
## walking up to the DESCRIPTION, so these scripts run from any working directory
## and on anyone's checkout, whether via Rscript or source().
if (nzchar(.l <- Sys.getenv("BLINK_LIB"))) .libPaths(.l)
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(dirname(.f), "/") else getwd()
while (!file.exists(file.path(ROOT, "DESCRIPTION")) && dirname(ROOT) != ROOT)
  ROOT <- dirname(ROOT)
if (!file.exists(file.path(ROOT, "DESCRIPTION")))
  stop("run this from inside the blink checkout (no DESCRIPTION found above ", getwd(), ")")
suppressMessages(library(malariasimulation))
source(file.path(ROOT, "comparison", "theme.R"))     # scenario constants
DDIR <- file.path(ROOT, "comparison", "data"); dir.create(DDIR, showWarnings = FALSE)
N_WORKERS <- 10L
## CMP_SMOKE=1 -> a few-minute end-to-end check: 4-year horizon, one replicate,
## interventions at year 1 so every builder actually fires, output to data/smoke/.
SMOKE <- nzchar(Sys.getenv("CMP_SMOKE"))
if (SMOKE) { N_REP <- 1L; N_WORKERS <- 3L; BURN_Y <- 1L
  DDIR <- file.path(DDIR, "smoke"); dir.create(DDIR, showWarnings = FALSE) }
PROG <- file.path(DDIR, "progress.log"); if (file.exists(PROG)) invisible(file.remove(PROG))
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

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
if (length(ONLY)) { stopifnot(all(ONLY %in% names(scenarios))); scenarios <- scenarios[ONLY] }
## CMP_BLINK_ONLY=1 -> re-run only blink (the IBM rows are kept): for a blink-side
## model change, when the IBM results are unaffected
BLINK_ONLY <- nzchar(Sys.getenv("CMP_BLINK_ONLY"))

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
  ## realised EIR (IBM emits EIR_<species> as total bites; blink emits per-adult-per-year)
  eir <- if (any(grepl("^EIR_", names(df)))) {
    ec <- grep("^EIR_", names(df), value = TRUE)
    sum(rowSums(df[obs, ec, drop = FALSE])) / length(obs) / POP * 365
  } else mean(df$EIR[obs])
  eq$eir_realised <- eir
  list(eq = eq, age = age, monthly = monthly, doy = doy)
}
tag_parts <- function(r, nm, model, rep)
  lapply(r, function(x) if (is.null(x)) NULL else cbind(scenario = nm, model = model, rep = rep, x))

## ---- run blink (deterministic, seconds) --------------------------------------
suppressMessages(library(blink))
log_msg("blink: %d scenarios", length(scenarios))
ode <- lapply(names(scenarios), function(nm) {
  s <- scenarios[[nm]]
  el <- system.time(
    ## s$p already carries init_EIR: every scenario is built through
    ## set_equilibrium(), which is where blink reads the target EIR from now.
    out <- run_simulation_ode(timesteps = s$years * 365, parameters = s$p,
                              tuning = list(rtol = 1e-6, step_size_max = 10)))[["elapsed"]]
  r <- summarise_run(out[-1, ], s$years)          # drop the day-0 seed row
  r$timing <- data.frame(years = s$years, elapsed_s = el)
  tag_parts(r, nm, "blink", 0L)
})
names(ode) <- names(scenarios)
log_msg("blink done: %.1f s total", sum(vapply(ode, function(o) o$timing$elapsed_s, numeric(1))))

## ---- run the IBM replicates in parallel --------------------------------------
ibm <- list()
if (!BLINK_ONLY) {
jobs <- expand.grid(scenario = names(scenarios), rep = seq_len(N_REP),
                    stringsAsFactors = FALSE)
## rough cost in default-band sim-years, so the load balancer starts the longest jobs first
jobs$cost <- vapply(jobs$scenario, function(nm) scenarios[[nm]]$years *
  (if (length(scenarios[[nm]]$p$prevalence_rendering_min_ages) > 1) 2 else 1), numeric(1))
jobs <- jobs[order(-jobs$cost, jobs$rep), ]; rownames(jobs) <- NULL
log_msg("IBM: %d runs (%d scenarios x %d reps) on %d workers; ~%.0f min at 2.5 s per sim-year",
        nrow(jobs), length(scenarios), N_REP, N_WORKERS, sum(jobs$cost) * 2.5 / N_WORKERS / 60)
cl <- parallel::makeCluster(N_WORKERS)
## workers do not inherit .libPaths(), so hand them the parent's rather than
## hardcoding one: whatever library this session is using, they use too
.libs <- .libPaths()
parallel::clusterExport(cl, c("jobs", "scenarios", "summarise_run", "tag_parts", "band_lo",
                              "band_hi", "AGE_TAGS", "POP", "BURN_Y", "PROG", ".libs"))
invisible(parallel::clusterEvalQ(cl, {
  .libPaths(.libs)
  suppressMessages(library(malariasimulation))
}))
t0 <- Sys.time()
ibm <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(j) {
  nm <- jobs$scenario[j]; k <- jobs$rep[j]; s <- scenarios[[nm]]
  set.seed(1000L + 7L * k)
  el <- system.time(out <- run_simulation(timesteps = s$years * 365, parameters = s$p))[["elapsed"]]
  r <- summarise_run(out, s$years)
  r$timing <- data.frame(years = s$years, elapsed_s = el)
  cat(sprintf("[%s] %-10s rep %d  %5.1f min\n", format(Sys.time(), "%H:%M:%S"), nm, k, el / 60),
      file = PROG, append = TRUE)
  tag_parts(r, nm, "IBM", k)
}, chunk.size = 1L)
parallel::stopCluster(cl)
log_msg("IBM done in %.1f min", as.numeric(Sys.time() - t0, units = "mins"))
}

## ---- combine and write -------------------------------------------------------
bind <- function(part) do.call(rbind, c(lapply(ode, `[[`, part), lapply(ibm, `[[`, part)))
for (part in c("eq", "age", "monthly", "doy", "timing")) {
  d <- bind(part); f <- file.path(DDIR, paste0("rep_", part, ".csv"))
  if ((length(ONLY) || BLINK_ONLY) && file.exists(f)) {   # partial run: replace just those rows
    old <- read.csv(f, stringsAsFactors = FALSE)
    drop <- old$scenario %in% names(scenarios) & (if (BLINK_ONLY) old$model == "blink" else TRUE)
    old <- old[!drop, ]
    for (nm in setdiff(names(d), names(old))) old[[nm]] <- NA
    for (nm in setdiff(names(old), names(d))) d[[nm]] <- NA
    d <- rbind(old[names(d)], d)
  }
  write.csv(d, f, row.names = FALSE)
  log_msg("wrote rep_%s.csv (%d rows)", part, nrow(d))
}
log_msg("ALL DONE")
