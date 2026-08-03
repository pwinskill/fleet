# Clinical & severe incidence time series: ODE (malariaode) vs IBM (malariasimulation).
# Reuses the (already-validated) dynamic scenario configs — bed nets, seasonality,
# seasonal SMC — but renders clinical AND severe incidence and extracts monthly,
# annualised rates for the under-5 band. Both models get the same rendering bands
# so their n_inc_clinical_* / n_inc_severe_* / n_age_* columns line up.
#
# Outputs: comparison/data/inc_*.csv and comparison/plots/inc_*.png.

.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")
suppressMessages({library(pkgload); library(malariasimulation)})
ROOT <- "C:/Users/pwinskil/Documents/dev/blink2/malariaode/comparison"
source(file.path(ROOT, "plot_helpers.R"))
pkgload::load_all("C:/Users/pwinskil/Documents/dev/blink2/malariaode", quiet = TRUE)
DDIR <- file.path(ROOT, "data"); PDIR <- file.path(ROOT, "plots")
set.seed(1L)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# clinical + severe + prevalence rendered on the same bands (0-5y, 2-10y, all-ages)
set_inc_bands <- function(p) {
  lo <- c(0, 730, 0); hi <- c(1825, 3650, 36500)
  p$clinical_incidence_rendering_min_ages <- lo; p$clinical_incidence_rendering_max_ages <- hi
  p$severe_incidence_rendering_min_ages   <- lo; p$severe_incidence_rendering_max_ages   <- hi
  p$prevalence_rendering_min_ages <- lo; p$prevalence_rendering_max_ages <- hi
  p
}
run_ibm <- function(p, years) run_simulation(timesteps = years * 365, parameters = p)
run_ode <- function(p, years, E) run_simulation_ode(years * 365, p, init_EIR = E)

# monthly annualised incidence for a band: sum(cases)/sum(person-days)*365*per
inc_ts <- function(df, kind, tag, per = 1) {
  cases <- df[[paste0("n_inc_", kind, "_", tag)]]
  pop   <- df[[paste0("n_age_", tag)]]
  stopifnot(!is.null(cases), !is.null(pop))
  n <- nrow(df); mon <- (seq_len(n) - 1) %/% 30
  agg <- tapply(seq_len(n), mon, function(ix) {
    pd <- sum(pop[ix]); if (pd > 0) sum(cases[ix]) / pd * 365 * per else NA_real_
  })
  data.frame(year = as.numeric(names(agg)) * 30 / 365, inc = as.numeric(agg))
}
# long incidence frame (both metrics x both models) for one band. The ODE table
# carries a day-0 seed row the IBM lacks; drop it so both bin on days 1..N.
inc_long <- function(ibm, ode, tag) {
  ode <- ode[-1, , drop = FALSE]
  mk <- function(df, model) rbind(
    cbind(inc_ts(df, "clinical", tag, 1),    model = model, metric = INC_LEVELS[1]),
    cbind(inc_ts(df, "severe",   tag, 1000), model = model, metric = INC_LEVELS[2]))
  d <- rbind(mk(ibm, "IBM"), mk(ode, "ODE"))
  d$metric <- factor(d$metric, levels = INC_LEVELS)
  d
}

BAND <- "0_1825"   # under-5 (where clinical/severe burden concentrates)

# --- nets ---------------------------------------------------------------------
scenario_nets <- function() {
  log_msg("incidence: bed nets")
  yrs <- BURN_Y + 6
  p <- set_inc_bands(get_parameters(list(human_population = POP)))
  p <- set_bednets(p, timesteps = BURN_Y * 365, coverages = 0.8, retention = 5 * 365,
                   dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
                   rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)
  p <- set_equilibrium(p, init_EIR = 20)
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 20)
  d <- inc_long(ibm, ode, BAND)
  write.csv(d, file.path(DDIR, "inc_nets.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "inc_nets.png"), render_incidence(INC_SPECS$nets, d),
         width = 8, height = 6, dpi = 130)
}

# --- seasonality --------------------------------------------------------------
scenario_seasonal <- function() {
  log_msg("incidence: seasonality")
  yrs <- BURN_Y + 3
  p <- set_inc_bands(get_parameters(list(human_population = POP,
    model_seasonality = TRUE, g0 = 0.285, g = c(-0.33, -0.13, 0.052),
    h = c(-0.35, 0.020, 0.10))))
  p <- set_equilibrium(p, init_EIR = 20)
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 20)
  d <- inc_long(ibm, ode, BAND)
  write.csv(d, file.path(DDIR, "inc_seasonal.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "inc_seasonal.png"), render_incidence(INC_SPECS$seasonal, d),
         width = 8, height = 6, dpi = 130)
}

# --- seasonal SMC -------------------------------------------------------------
scenario_smc <- function() {
  log_msg("incidence: SMC")
  yrs <- BURN_Y + 3
  p <- set_inc_bands(get_parameters(list(human_population = POP,
    model_seasonality = TRUE, g0 = 0.285, g = c(-0.33, -0.13, 0.052),
    h = c(-0.35, 0.020, 0.10))))
  p <- set_drugs(p, list(SP_AQ_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.45)
  p <- set_smc(p, drug = 1, timesteps = SMC_ROUNDS, coverages = rep(0.9, length(SMC_ROUNDS)),
               min_ages = rep(round(0.25 * 365), length(SMC_ROUNDS)),
               max_ages = rep(round(5 * 365), length(SMC_ROUNDS)))
  p <- set_equilibrium(p, init_EIR = 15)
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 15)
  d <- inc_long(ibm, ode, BAND)
  write.csv(d, file.path(DDIR, "inc_smc.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "inc_smc.png"), render_incidence(INC_SPECS$smc, d),
         width = 8, height = 6, dpi = 130)
}

t0 <- Sys.time()
scenario_nets(); scenario_seasonal(); scenario_smc()
log_msg("incidence DONE in %.1f min", as.numeric(Sys.time() - t0, units = "mins"))
