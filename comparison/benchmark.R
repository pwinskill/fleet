# Indicative run times for fleet, for the table in README.md (issue #1).
#
#   Rscript comparison/benchmark.R          # ~6 min; writes comparison/data/timing.csv
#                                           # and prints the markdown tables
#
# Three questions a user actually has:
#   1. How long does a run take, by scenario and horizon?
#   2. Does population size matter?  (it must not -- the ODE is per-capita)
#   3. What do the solver settings cost?
#
# Every figure is the MINIMUM of N_REP repeats of a complete
# run_simulation_ode() call, including building the inputs and seeding the
# equilibrium, because that is what a user pays. The minimum, not the mean or
# median: contention from anything else on the machine can only ADD time, so the
# fastest repeat is the least contaminated estimate of the cost of the work. (An
# earlier median-of-3 pass, run while an install was going on, reported the same
# 30-year seasonal configuration as 15.7 s in one table and 5.8 s in another.)
# Run this on an otherwise idle machine. Timings are wall-clock on one core;
# nothing here is parallel.

## No absolute paths anywhere in here. FLEET_LIB is prepended to the library
## path, for installations that do not pick up R_LIBS_USER (the Windows-arm64
## setup this was developed on); the libraries already on the path are kept, so a
## FLEET_LIB holding only some of the dependencies still works. Leave it unset and
## your normal library is used. ROOT is found by walking up to the DESCRIPTION, so
## these scripts run from any working directory and on anyone's checkout, whether
## via Rscript or source().
if (nzchar(.l <- Sys.getenv("FLEET_LIB"))) .libPaths(c(.l, .libPaths()))
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(dirname(.f), "/") else getwd()
while (!file.exists(file.path(ROOT, "DESCRIPTION")) && dirname(ROOT) != ROOT)
  ROOT <- dirname(ROOT)
if (!file.exists(file.path(ROOT, "DESCRIPTION")))
  stop("run this from inside the fleet checkout (no DESCRIPTION found above ", getwd(), ")")
suppressMessages({library(malariasimulation); library(fleet)})
DDIR <- file.path(ROOT, "comparison", "data")
N_REP <- 5L
POP   <- 1e5                       # output scaling only; see the population table
SEASON <- list(g0 = 0.285, g = c(-0.33, -0.13, 0.052), h = c(-0.35, 0.020, 0.10))
## the fast preset documented in ?ode_tuning (atol stays at its default)
FAST <- list(rtol = 1e-6, step_size_max = 10)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

## ---- scenarios ---------------------------------------------------------------
## Each builds a parameter list at EIR 20 (SMC scenarios at 15, where SMC is
## actually deployed). YEARS is supplied at run time.
base_p <- function(seasonal = FALSE, pop = POP) {
  ov <- list(human_population = pop)
  if (seasonal) ov <- c(ov, list(model_seasonality = TRUE), SEASON)
  get_parameters(ov)
}
add_nets_irs <- function(p, years) {
  n <- max(1L, years %/% 3L)
  ts <- seq(365, by = 3 * 365, length.out = n)
  ## every time-varying net argument must have one entry per distribution
  p <- set_bednets(p, timesteps = ts, coverages = rep(0.7, n), retention = 5 * 365,
                   dn0 = matrix(0.387, n), rn = matrix(0.563, n),
                   rnm = matrix(0.24, n), gamman = rep(2.64 * 365, n))
  m <- function(v) matrix(v, nrow = n, ncol = 1)
  set_spraying(p, timesteps = ts, coverages = rep(0.8, n),
               ls_theta = m(2.025), ls_gamma = m(-0.009),
               ks_theta = m(-2.222), ks_gamma = m(0.008),
               ms_theta = m(-1.232), ms_gamma = m(-0.009))
}
add_smc <- function(p, years) {
  p <- set_drugs(p, list(SP_AQ_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.45)
  rounds <- as.vector(sapply(seq_len(years) - 1L,
                             function(y) 365 + y * 365 + c(0, 30, 60, 90) + 200))
  rounds <- rounds[rounds <= years * 365]
  set_smc(p, drug = 1, timesteps = rounds, coverages = rep(0.9, length(rounds)),
          min_ages = rep(round(0.25 * 365), length(rounds)),
          max_ages = rep(round(5 * 365), length(rounds)))
}
add_pev <- function(p) {
  set_pev_epi(p, profile = rtss_profile, timesteps = 365, coverages = 0.9,
              min_wait = 0, age = 5 * 30, booster_spacing = 12 * 30,
              booster_coverage = matrix(0.8), booster_profile = list(rtss_booster_profile))
}

SCEN <- list(
  `No interventions`                 = function(y) list(p = base_p(), eir = 20),
  `Seasonal`                         = function(y) list(p = base_p(TRUE), eir = 20),
  `Seasonal + treatment (AL)`        = function(y) {
    p <- set_drugs(base_p(TRUE), list(AL_params))
    list(p = set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4), eir = 20) },
  `Seasonal + nets + IRS`            = function(y) list(p = add_nets_irs(base_p(TRUE), y), eir = 20),
  `Seasonal + SMC (4 rounds/yr)`     = function(y) list(p = add_smc(base_p(TRUE), y), eir = 15),
  `All of the above + RTS,S`         = function(y) {
    p <- add_pev(add_nets_irs(add_smc(base_p(TRUE), y), y))
    list(p = p, eir = 15) }
)
YEARS <- c(5L, 10L, 30L)

## ---- timing helper -----------------------------------------------------------
## the whole user-visible call: build inputs, seed, integrate, render outputs
## p already carries init_EIR (every caller builds it through set_equilibrium),
## which is where run_simulation_ode reads the target EIR from.
time_run <- function(p, years, ctrl = FAST, reps = N_REP) {
  el <- numeric(reps); gc(verbose = FALSE)
  for (k in seq_len(reps)) {
    el[k] <- system.time(invisible(run_simulation_ode(
      timesteps = years * 365, parameters = p, tuning = ctrl)))[["elapsed"]]
  }
  min(el)
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

## ---- table 1: scenario x horizon ---------------------------------------------
log_msg("scenario x horizon: %d scenarios x %d horizons x %d reps",
        length(SCEN), length(YEARS), N_REP)
for (nm in names(SCEN)) {
  for (y in YEARS) {
    s <- SCEN[[nm]](y)
    sec <- time_run(set_equilibrium(s$p, init_EIR = s$eir), y)
    add_row(table = "scenario", scenario = nm, years = y, pop = POP,
            settings = "fast", seconds = sec)
    log_msg("  %-30s %2d y : %6.2f s", nm, y, sec)
  }
}

## ---- table 2: population size (must be flat) ---------------------------------
log_msg("population independence, 30-year seasonal run")
for (pop in c(1e3, 1e4, 1e5, 1e6, 1e7)) {
  p <- set_equilibrium(base_p(TRUE, pop = pop), init_EIR = 20)
  sec <- time_run(p, 30L)
  add_row(table = "population", scenario = "Seasonal", years = 30L, pop = pop,
          settings = "fast", seconds = sec)
  log_msg("  pop %-9s : %6.2f s", format(pop, scientific = TRUE), sec)
}

## ---- table 3: solver settings ------------------------------------------------
log_msg("solver settings, 30-year seasonal run")
PRESETS <- list(
  `default (atol = rtol = 1e-8, step <= 1 d)`  = list(atol = 1e-8, rtol = 1e-8, step_size_max = 1),
  `fast (atol = 1e-8, rtol = 1e-6, step <= 10 d)` = FAST,
  `loose (atol = rtol = 1e-6, step <= 10 d)`   = list(atol = 1e-6, rtol = 1e-6, step_size_max = 10)
)
for (nm in names(PRESETS)) {
  p <- set_equilibrium(base_p(TRUE), init_EIR = 20)
  sec <- time_run(p, 30L, ctrl = PRESETS[[nm]])
  add_row(table = "settings", scenario = nm, years = 30L, pop = POP,
          settings = nm, seconds = sec)
  log_msg("  %-46s : %6.2f s", nm, sec)
}

## ---- write + print -----------------------------------------------------------
timing <- do.call(rbind, rows)
timing$seconds <- round(timing$seconds, 2)
timing$s_per_year <- round(timing$seconds / timing$years, 3)
write.csv(timing, file.path(DDIR, "timing.csv"), row.names = FALSE)

md <- function(df) {
  df[] <- lapply(df, as.character)
  cat(paste0("| ", paste(names(df), collapse = " | "), " |\n"))
  cat(paste0("|", paste(rep(" --- ", ncol(df)), collapse = "|"), "|\n"))
  for (i in seq_len(nrow(df))) cat(paste0("| ", paste(df[i, ], collapse = " | "), " |\n"))
  cat("\n")
}
cat("\n## Scenario x horizon (seconds)\n\n")
sc <- timing[timing$table == "scenario", ]
wide <- reshape(sc[, c("scenario", "years", "seconds")], idvar = "scenario",
                timevar = "years", direction = "wide")
names(wide) <- c("Scenario", paste0(YEARS, " years"))
wide$Scenario <- factor(wide$Scenario, levels = names(SCEN))
md(wide[order(wide$Scenario), ])

cat("## Population size, 30-year seasonal run\n\n")
pp <- timing[timing$table == "population", c("pop", "seconds")]
md(data.frame(Population = format(pp$pop, big.mark = ",", scientific = FALSE, trim = TRUE),
              Seconds = pp$seconds))

cat("## Solver settings, 30-year seasonal run\n\n")
st <- timing[timing$table == "settings", c("settings", "seconds")]
md(data.frame(Settings = st$settings, Seconds = st$seconds))

cat("machine: ", R.version$version.string, " / ", R.version$platform, "\n", sep = "")
log_msg("wrote comparison/data/timing.csv (%d rows)", nrow(timing))
