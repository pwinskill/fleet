# ODE (blink) vs IBM (malariasimulation) comparison.
#
# For each scenario the SAME parameter list is passed to both models: the IBM is
# burned in ~30 years to settle onto its own dynamic steady state, the ODE (seeded
# at the malariaEquilibrium fixed point) holds/relaxes over the same horizon.
# Prevalence age bands are set once on the parameter list so both models emit the
# same n_detect_lm_<lo>_<hi> / n_age_<lo>_<hi> columns.
#
# Outputs: comparison/data/*.csv (tidy) and comparison/plots/*.png.

.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")
suppressMessages({
  library(pkgload); library(malariasimulation)
  library(ggplot2); library(dplyr); library(tidyr)
})
pkgload::load_all("C:/Users/pwinskil/Documents/dev/blink2/blink", quiet = TRUE)

ROOT   <- "C:/Users/pwinskil/Documents/dev/blink2/blink/comparison"
DDIR   <- file.path(ROOT, "data"); PDIR <- file.path(ROOT, "plots")
# shared figure builders + constants (POP, BURN_Y, AGE_EDGES, band_lo/hi, themes)
source(file.path(ROOT, "plot_helpers.R"))
set.seed(1L)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

add_render_bands <- function(p) {
  # prevalence bands = the age profile + the canonical 2-10y (730..3650) band
  p$prevalence_rendering_min_ages <- c(band_lo, 730)
  p$prevalence_rendering_max_ages <- c(band_hi, 3650)
  p$clinical_incidence_rendering_min_ages <- c(0.25 * 365, 730)
  p$clinical_incidence_rendering_max_ages <- c(5 * 365, 3650)
  p
}

# banded LM prevalence from an output frame, pooled over a row window
pfpr <- function(df, tag, rows = seq_len(nrow(df))) {
  nd <- df[[paste0("n_detect_lm_", tag)]][rows]
  na <- df[[paste0("n_age_", tag)]][rows]
  sum(nd) / sum(na)
}
tag_of <- function(lo, hi) paste0(round(lo), "_", round(hi))

# IBM realised per-person EIR (bites/person/year) from EIR_<species> columns
ibm_eir <- function(df, rows = seq_len(nrow(df))) {
  ec <- grep("^EIR_", names(df), value = TRUE)
  sum(rowSums(df[rows, ec, drop = FALSE])) / length(rows) / POP * 365
}

run_ibm <- function(p, years) run_simulation(timesteps = years * 365, parameters = p)
run_ode <- function(p, years, init_EIR)
  run_simulation_ode(timesteps = years * 365, parameters = p, init_EIR = init_EIR)

# ============================================================================
# Scenario A — equilibrium PfPR(2-10) vs EIR
# ============================================================================
scenario_A <- function() {
  log_msg("A: PfPR-EIR curve")
  EIRs <- c(1, 3, 10, 20, 50, 120)
  rows_obs <- function(df) (nrow(df) - 3 * 365 + 1):nrow(df)   # last 3 years
  res <- lapply(EIRs, function(E) {
    p <- add_render_bands(get_parameters(list(human_population = POP)))
    p <- set_equilibrium(p, init_EIR = E)
    ibm <- run_ibm(p, BURN_Y + 3)
    ode <- run_ode(p, 5, E)
    log_msg("   EIR %.0f done (IBM realised EIR=%.1f)", E, ibm_eir(ibm, rows_obs(ibm)))
    data.frame(init_EIR = E,
               IBM = pfpr(ibm, "730_3650", rows_obs(ibm)),
               ODE = pfpr(ode, "730_3650"),
               IBM_EIR = ibm_eir(ibm, rows_obs(ibm)))
  })
  d <- do.call(rbind, res)
  write.csv(d, file.path(DDIR, "A_pfpr_eir.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "A_pfpr_eir.png"), plot_A(d), width = 7, height = 5, dpi = 130)
  d
}

# ============================================================================
# Scenario B — age-prevalence profile at EIR 20
# ============================================================================
scenario_B <- function() {
  log_msg("B: age-prevalence profile @ EIR 20")
  p <- add_render_bands(get_parameters(list(human_population = POP)))
  p <- set_equilibrium(p, init_EIR = 20)
  ibm <- run_ibm(p, BURN_Y + 3)
  ode <- run_ode(p, 5, 20)
  rows_obs <- (nrow(ibm) - 3 * 365 + 1):nrow(ibm)
  tags <- tag_of(band_lo, band_hi)
  d <- data.frame(
    age_mid = (band_lo + band_hi) / 2 / 365,
    IBM = vapply(tags, function(tg) pfpr(ibm, tg, rows_obs), numeric(1)),
    ODE = vapply(tags, function(tg) pfpr(ode, tg), numeric(1)))
  write.csv(d, file.path(DDIR, "B_age_profile.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "B_age_profile.png"), plot_B(d), width = 7, height = 5, dpi = 130)
  d
}

# monthly-smoothed 2-10 prevalence time series (years on x-axis)
prev_ts <- function(df) {
  n <- nrow(df); mon <- ((seq_len(n) - 1) %/% 30)
  agg <- tapply(seq_len(n), mon, function(ix)
    sum(df$n_detect_lm_730_3650[ix]) / sum(df$n_age_730_3650[ix]))
  data.frame(year = (as.numeric(names(agg)) * 30) / 365, pfpr = as.numeric(agg))
}

# ============================================================================
# Scenario C — bed-net campaign deployed after burn-in
# ============================================================================
scenario_C <- function() {
  log_msg("C: bed nets at year %d", BURN_Y)
  yrs <- BURN_Y + 6
  p <- add_render_bands(get_parameters(list(human_population = POP)))
  p <- set_bednets(p, timesteps = BURN_Y * 365, coverages = 0.8, retention = 5 * 365,
                   dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
                   rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)
  p <- set_equilibrium(p, init_EIR = 20)
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 20)
  a <- prev_ts(ibm); a$model <- "IBM"; b <- prev_ts(ode); b$model <- "ODE"
  d <- rbind(a, b)
  write.csv(d, file.path(DDIR, "C_bednets.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "C_bednets.png"), plot_C(d), width = 8, height = 5, dpi = 130)
  d
}

# ============================================================================
# Scenario D — seasonal transmission, settled annual cycle
# ============================================================================
scenario_D <- function() {
  log_msg("D: seasonality")
  p <- add_render_bands(get_parameters(list(human_population = POP,
    model_seasonality = TRUE, g0 = 0.285, g = c(-0.33, -0.13, 0.052),
    h = c(-0.35, 0.020, 0.10))))
  p <- set_equilibrium(p, init_EIR = 20)
  yrs <- BURN_Y + 2
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 20)
  # final settled year, day-of-year prevalence
  doy <- function(df) {
    rows <- (nrow(df) - 365 + 1):nrow(df)
    data.frame(doy = seq_len(365),
               pfpr = df$n_detect_lm_730_3650[rows] / df$n_age_730_3650[rows])
  }
  a <- doy(ibm); a$model <- "IBM"; b <- doy(ode); b$model <- "ODE"
  # smooth the IBM daily noise with a 15-day rolling mean for display
  sm <- function(x, k = 15) stats::filter(c(x, x, x), rep(1 / k, k))[(length(x) + 1):(2 * length(x))]
  a$pfpr <- as.numeric(sm(a$pfpr))
  d <- rbind(a, b)
  write.csv(d, file.path(DDIR, "D_seasonal.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "D_seasonal.png"), plot_D(d), width = 8, height = 5, dpi = 130)
  d
}

# ============================================================================
# Scenario E — seasonal SMC in young children
# ============================================================================
scenario_E <- function() {
  log_msg("E: SMC")
  p <- add_render_bands(get_parameters(list(human_population = POP,
    model_seasonality = TRUE, g0 = 0.285, g = c(-0.33, -0.13, 0.052),
    h = c(-0.35, 0.020, 0.10))))
  p <- set_drugs(p, list(SP_AQ_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.45)
  # 4 monthly SMC rounds each year for 3 years after burn-in (peak-season)
  rounds <- as.vector(sapply(0:2, function(y) BURN_Y * 365 + y * 365 +
                               c(0, 30, 60, 90) + 200))
  p <- set_smc(p, drug = 1, timesteps = rounds, coverages = rep(0.9, length(rounds)),
               min_ages = rep(round(0.25 * 365), length(rounds)),
               max_ages = rep(round(5 * 365), length(rounds)))
  p <- set_equilibrium(p, init_EIR = 15)
  yrs <- BURN_Y + 3
  ibm <- run_ibm(p, yrs); ode <- run_ode(p, yrs, 15)
  # prevalence in the SMC target band (0.25-5y ~ nearest rendered band 0-5)
  ts_band <- function(df) {
    n <- nrow(df); mon <- ((seq_len(n) - 1) %/% 30)
    tg1 <- tag_of(0 * 365, 1 * 365); tg2 <- tag_of(1 * 365, 2 * 365)
    tg3 <- tag_of(2 * 365, 3 * 365); tg4 <- tag_of(3 * 365, 5 * 365)
    nd <- df[[paste0("n_detect_lm_", tg1)]] + df[[paste0("n_detect_lm_", tg2)]] +
          df[[paste0("n_detect_lm_", tg3)]] + df[[paste0("n_detect_lm_", tg4)]]
    na <- df[[paste0("n_age_", tg1)]] + df[[paste0("n_age_", tg2)]] +
          df[[paste0("n_age_", tg3)]] + df[[paste0("n_age_", tg4)]]
    agg <- tapply(seq_len(n), mon, function(ix) sum(nd[ix]) / sum(na[ix]))
    data.frame(year = as.numeric(names(agg)) * 30 / 365, pfpr = as.numeric(agg))
  }
  a <- ts_band(ibm); a$model <- "IBM"; b <- ts_band(ode); b$model <- "ODE"
  d <- rbind(a, b)
  write.csv(d, file.path(DDIR, "E_smc.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "E_smc.png"), plot_E(d), width = 8, height = 5, dpi = 130)
  d
}

# ============================================================================
# Scenario F — custom demography: age structure + age-prevalence
# ============================================================================
scenario_F <- function() {
  log_msg("F: custom demography")
  dr <- c(0.048, 0.007, 0.003, 0.004, 0.008, 0.020, 0.050, 0.120) / 365
  ag <- round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365)
  p <- add_render_bands(get_parameters(list(human_population = POP)))
  p <- set_demography(p, agegroups = ag, timesteps = 0, deathrates = matrix(dr, nrow = 1))
  p <- set_equilibrium(p, init_EIR = 20)
  ibm <- run_ibm(p, BURN_Y + 3); ode <- run_ode(p, 5, 20)
  rows_obs <- (nrow(ibm) - 3 * 365 + 1):nrow(ibm)
  tags <- tag_of(band_lo, band_hi); amid <- (band_lo + band_hi) / 2 / 365
  # age structure (fraction of population in each band)
  struct <- function(df, rows) {
    na <- vapply(tags, function(tg) mean(df[[paste0("n_age_", tg)]][rows]), numeric(1))
    na / sum(na)
  }
  ds <- data.frame(age_mid = amid, IBM = struct(ibm, rows_obs),
                   ODE = struct(ode, seq_len(nrow(ode))))
  dp <- data.frame(age_mid = amid,
                   IBM = vapply(tags, function(tg) pfpr(ibm, tg, rows_obs), numeric(1)),
                   ODE = vapply(tags, function(tg) pfpr(ode, tg), numeric(1)))
  write.csv(ds, file.path(DDIR, "F_demog_structure.csv"), row.names = FALSE)
  write.csv(dp, file.path(DDIR, "F_demog_prevalence.csv"), row.names = FALSE)
  ggsave(file.path(PDIR, "F_demog_structure.png"), plot_F_structure(ds), width = 7, height = 5, dpi = 130)
  ggsave(file.path(PDIR, "F_demog_prevalence.png"), plot_F_prev(dp), width = 7, height = 5, dpi = 130)
  list(structure = ds, prevalence = dp)
}

# ---- run all ----------------------------------------------------------------
t0 <- Sys.time()
A <- scenario_A(); B <- scenario_B(); C <- scenario_C()
D <- scenario_D(); E <- scenario_E(); F <- scenario_F()
log_msg("ALL DONE in %.1f min", as.numeric(Sys.time() - t0, units = "mins"))

cat("\n=== Scenario A (PfPR2-10 vs EIR) ===\n"); print(A, digits = 3)
cat("\n=== Scenario B (age profile @20) — abs diff summary ===\n")
cat(sprintf("max |IBM-ODE| = %.3f, mean = %.3f\n",
            max(abs(B$IBM - B$ODE)), mean(abs(B$IBM - B$ODE))))
cat("\n=== Scenario F structure — under-5 fraction ===\n")
cat(sprintf("IBM=%.3f ODE=%.3f\n", sum(F$structure$IBM[F$structure$age_mid < 5]),
            sum(F$structure$ODE[F$structure$age_mid < 5])))
