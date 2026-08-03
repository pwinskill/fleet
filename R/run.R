.malariaode_env <- new.env(parent = emptyenv())

#' Locate the odin model source (installed or in-development).
#' @noRd
odin_model_path <- function() {
  path <- system.file("odin", "malaria_ode.R", package = "malariaode")
  if (nzchar(path) && file.exists(path)) return(path)
  # development fallbacks
  for (cand in c(
    file.path("inst", "odin", "malaria_ode.R"),
    file.path(getwd(), "inst", "odin", "malaria_ode.R")
  )) {
    if (file.exists(cand)) return(cand)
  }
  stop("Could not locate inst/odin/malaria_ode.R")
}

#' Compile (and cache) the odin2 generator for the mean-field model.
#' @param odin_file optional explicit path to the odin source.
#' @noRd
get_generator <- function(odin_file = NULL) {
  if (is.null(odin_file)) odin_file <- odin_model_path()
  key <- normalizePath(odin_file, mustWork = TRUE)
  if (is.null(.malariaode_env$generators)) .malariaode_env$generators <- list()
  if (is.null(.malariaode_env$generators[[key]])) {
    .malariaode_env$generators[[key]] <- odin2::odin(odin_file, quiet = TRUE)
  }
  .malariaode_env$generators[[key]]
}

#' Run the mean-field (ODE) malaria model.
#'
#' @param timesteps number of days to simulate. Output has `timesteps + 1` rows
#'   (day 0 is the seeded equilibrium; days 1..timesteps are integrated).
#' @param parameters a malariasimulation::get_parameters() list (falciparum),
#'   optionally with interventions layered on via the `set_*` builders. The
#'   following are applied automatically from the list (no extra arguments):
#'   clinical treatment (time-varying), antimalarial resistance, bed nets, IRS,
#'   MDA/SMC/PMC, PEV (EPI + mass) and TBV vaccines, seasonality and flexible
#'   carrying capacity, and custom demography (`set_demography`: age-specific
#'   mortality and the resulting equilibrium age structure). P. vivax is rejected;
#'   the model is always compartmental (the individual-mosquito path does not
#'   apply). See the README for the mean-field approximations used by each module.
#' @details Seasonal runs oscillate around a limit cycle rather than holding flat;
#'   the state is seeded at the annual-mean (aseasonal) equilibrium, so the first
#'   ~10 years are a transient onto the cycle and the seasonal annual-mean EIR
#'   sits a few percent below the aseasonal `init_EIR` target (nonlinear
#'   averaging). Use a burned-in cycle for calibration/comparison.
#' @param init_EIR target adult EIR (bites/adult/year). If NULL, taken from
#'   `parameters$init_EIR` (set by malariasimulation::set_equilibrium()).
#' @param age_lower age-group lower edges in **years** (default graded grid).
#' @param n_eir,n_foim,n_eip Erlang-chain stage counts for the EIR lag, FOIM lag
#'   and mosquito EIP. Larger values sharpen the (otherwise gamma-shaped) lags
#'   toward the IBM's fixed delays; equilibrium is exact for any value.
#' @param atol,rtol,step_size_max dust2 ODE-solver controls. The defaults
#'   (`1e-8`, `1e-8`, `1`) preserve the flat equilibrium exactly; for long dynamic
#'   projections a looser tolerance and larger step cap (e.g. `1e-6`, `1e-6`, `10`)
#'   run several times faster with negligible effect on aggregate outputs.
#' @param odin_file optional path to the odin source (development use).
#' @return a wide, malariasimulation-style daily count table, one row per output
#'   day. Columns:
#'   \itemize{
#'     \item `timestep` — output day (0..timesteps), the row key.
#'     \item per age band (tags in **days**, `<lo>_<hi>`): `n_age_*` (population),
#'       `n_detect_lm_*` and `p_detect_lm_*` (LM-positive count and PfPR
#'       proportion), `n_detect_pcr_*`, `n_inc_clinical_*`, `n_inc_severe_*` and
#'       `n_inc_*` (all-infection incidence) — all per-day counts.
#'     \item `ft` (treated fraction), `EIR` (per adult per year), `FOIM`.
#'     \item population-total infection-state counts `S_count`, `D_count`,
#'       `A_count`, `U_count`, `Tr_count`, `Ph_count` (diagnostic).
#'   }
#'   Incidence columns are per-day counts; do not thin rows before postie. Can be
#'   passed directly to postie::get_rates() / postie::get_prevalence().
#' @examples
#' \dontrun{
#' p <- malariasimulation::get_parameters()
#'
#' # daily wide count table, init_EIR supplied directly
#' out <- run_simulation_ode(timesteps = 3650, parameters = p, init_EIR = 20)
#' head(out[, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650", "EIR")])
#'
#' # or let set_equilibrium() seed init_EIR into the parameter list
#' p2 <- malariasimulation::set_equilibrium(
#'   malariasimulation::get_parameters(), init_EIR = 5)
#' out2 <- run_simulation_ode(timesteps = 3650, parameters = p2)
#' }
#' @export
run_simulation_ode <- function(timesteps, parameters, init_EIR = NULL,
                               age_lower = default_age_lower(),
                               n_eir = 10L, n_foim = 10L, n_eip = 20L,
                               atol = 1e-8, rtol = 1e-8, step_size_max = 1,
                               odin_file = NULL) {
  if (is.null(init_EIR)) {
    init_EIR <- parameters$init_EIR
    if (is.null(init_EIR)) {
      stop("Provide `init_EIR` or run malariasimulation::set_equilibrium() first.")
    }
  }
  inp <- build_inputs(parameters, init_EIR, age_lower, n_eir, n_foim, n_eip,
                      timesteps = timesteps)
  generator <- get_generator(odin_file)
  # Robustness: the interpolate grids (pev/tbv "linear") are built to extend past
  # `timesteps` so the stepper never extrapolates them. step_size_max caps the step;
  # at the default 1 it binds only in the near-equilibrium regime (derivatives ~ 0,
  # where the adaptive stepper would otherwise propose a huge step past an interpolate
  # knot). For long dynamic projections a larger cap + looser atol/rtol trade a little
  # accuracy for speed; the equilibrium-preserving defaults (1, 1e-8) are unchanged.
  ctrl <- dust2::dust_ode_control(max_steps = max(1e5, 100 * timesteps),
                                  atol = atol, rtol = rtol, step_size_max = step_size_max)
  sys <- dust2::dust_system_create(generator, pars = inp$pars, n_particles = 1,
                                   ode_control = ctrl)
  dust2::dust_system_set_state_initial(sys)
  times <- seq(0, timesteps)
  events <- chemoprevention_events(inp$meta$parameters, timesteps)
  if (length(events) == 0) {
    y <- dust2::dust_system_simulate(sys, times)
  } else {
    y <- simulate_with_pulses(sys, times, events, inp$meta)
  }
  render_output(sys, y, times, inp)
}

#' Integrate with chemoprevention pulses applied between segments.
#' @noRd
simulate_with_pulses <- function(sys, times, events, meta) {
  uidx <- dust2::dust_unpack_index(sys)
  timesteps <- max(times)
  ev_days <- sort(unique(vapply(events, function(e) e$time, numeric(1))))
  ev_days <- ev_days[ev_days >= 1 & ev_days <= timesteps]
  bnds <- sort(unique(c(0, ev_days, timesteps)))
  cols <- vector("list", length(bnds) - 1)
  recorded_hi <- -1
  for (k in seq_len(length(bnds) - 1)) {
    lo <- bnds[k]; hi <- bnds[k + 1]
    seg_days <- times[times > recorded_hi & times <= hi]
    if (length(seg_days)) {
      cols[[k]] <- dust2::dust_system_simulate(sys, seg_days)
      recorded_hi <- max(seg_days)
    }
    if (hi %in% ev_days) {
      for (e in events[vapply(events, function(x) x$time == hi, logical(1))]) {
        apply_chemoprevention_pulse(sys, uidx, meta, e)
      }
    }
  }
  do.call(cbind, cols[!vapply(cols, is.null, logical(1))])
}

#' Default output age bands (days), taken from the parameter list's rendering
#' fields where present, always including 2-10y and all-ages.
#' @noRd
output_bands <- function(p) {
  gather <- function(lo, hi) {
    lo <- p[[lo]]; hi <- p[[hi]]
    if (is.null(lo) || length(lo) == 0) return(NULL)
    Map(c, lo, hi)
  }
  bands <- c(
    gather("prevalence_rendering_min_ages", "prevalence_rendering_max_ages"),
    gather("incidence_rendering_min_ages", "incidence_rendering_max_ages"),
    gather("clinical_incidence_rendering_min_ages", "clinical_incidence_rendering_max_ages"),
    gather("severe_incidence_rendering_min_ages", "severe_incidence_rendering_max_ages"),
    gather("age_group_rendering_min_ages", "age_group_rendering_max_ages"),
    list(c(730, 3650), c(0, 36500))
  )
  unique(bands)
}

#' Render dust2 output into a wide malariasimulation-style count table that can
#' be passed directly to postie::get_rates() / postie::get_prevalence().
#' @noRd
render_output <- function(sys, y, times, inp) {
  st <- dust2::dust_unpack_state(sys, y)
  meta <- inp$meta
  hp <- meta$human_population
  age_mid <- meta$age_mid
  nt <- length(times)
  # per-age-group time series are [n_age, nt] matrices
  mat <- function(a) if (is.null(dim(a))) matrix(a, nrow = 1) else a
  n_g <- mat(st$n_g); det_lm_g <- mat(st$det_lm_g); det_pcr_g <- mat(st$det_pcr_g)
  clin_g <- mat(st$clin_g); sev_g <- mat(st$sev_g); inc_g <- mat(st$inc_g)

  band_sum <- function(m, lo, hi) {
    idx <- which(age_mid >= lo & age_mid < hi)
    if (length(idx) == 0) return(rep(0, ncol(m)))
    colSums(m[idx, , drop = FALSE])
  }

  out <- data.frame(timestep = times)
  for (b in output_bands(meta$parameters)) {
    lo <- b[1]; hi <- b[2]
    tag <- paste0(round(lo), "_", round(hi))
    n_age  <- band_sum(n_g, lo, hi) * hp
    n_lm   <- band_sum(det_lm_g, lo, hi) * hp
    n_pcr  <- band_sum(det_pcr_g, lo, hi) * hp
    n_clin <- band_sum(clin_g, lo, hi) * hp   # per-day incidence count
    n_sev  <- band_sum(sev_g, lo, hi) * hp
    out[[paste0("n_age_", tag)]] <- n_age
    out[[paste0("n_detect_lm_", tag)]] <- n_lm
    out[[paste0("p_detect_lm_", tag)]] <- ifelse(n_age > 0, n_lm / n_age, 0)
    out[[paste0("n_detect_pcr_", tag)]] <- n_pcr
    out[[paste0("n_inc_clinical_", tag)]] <- n_clin
    out[[paste0("n_inc_severe_", tag)]] <- n_sev
    out[[paste0("n_inc_", tag)]] <- band_sum(inc_g, lo, hi) * hp   # all-infection incidence
  }
  out$ft <- as.numeric(st$ft_out)
  out$EIR <- as.numeric(st$EIR_yr)
  out$FOIM <- as.numeric(st$FOIM)
  # population-total infection-state fractions (diagnostic)
  out$S_count <- colSums(mat(st$S_g)) * hp
  out$D_count <- colSums(mat(st$D_g)) * hp
  out$A_count <- colSums(mat(st$A_g)) * hp
  out$U_count <- colSums(mat(st$U_g)) * hp
  out$Tr_count <- colSums(mat(st$Tr_g)) * hp
  out$Ph_count <- colSums(mat(st$Ph_g)) * hp
  out
}
