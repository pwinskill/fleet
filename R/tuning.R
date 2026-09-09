# Solver and discretisation settings for run_simulation_ode().
#
# These are gathered into one object rather than sitting on the run signature so
# that run_simulation_ode(timesteps, parameters, correlations) matches
# malariasimulation::run_simulation() argument for argument. Nothing
# epidemiological belongs here: every field is a numerical-approximation knob,
# and changing one should move the answer only by its own error term. Model
# parameters, the target EIR included, live on the parameter list.

#' ODE solver and discretisation settings.
#'
#' Numerical settings for [run_simulation_ode()]. Every default is the validated
#' choice; a run with `ode_tuning()` untouched is the reference configuration.
#' Pass either this object or a plain named list of the fields you want to
#' change; anything left out keeps its default.
#'
#' @param age_lower age-group lower edges in **years** (default graded grid:
#'   monthly to 1 year, quarterly to 5, yearly to 15, then 5-yearly to an
#'   absorbing top group). Must start at 0 and increase strictly.
#' @param n_eir,n_foim,n_eip Erlang-chain stage counts for the EIR lag, FOIM lag
#'   and mosquito EIP. Larger values sharpen the (otherwise gamma-shaped) lags
#'   toward the IBM's fixed delays; equilibrium is exact for any value.
#' @param n_ph,n_phc Erlang-chain stage counts for the post-treatment (`Ph`) and
#'   chemoprevention (`Ph_c`) prophylaxis compartments. `NULL` (default) matches the
#'   chain's variance to the drug's Weibull protection curve, capped at 20. For
#'   `Ph_c` that is `1/CV²` of the Weibull: 14 for SP-AQ, 15 for DHA-PQP. `Ph`
#'   follows the exponential treated stage `Tr`, so its count matches the variance
#'   of the whole `Tr + Ph` sojourn and its mean is the integrated protection left
#'   after `Tr`: 16 stages for SP-AQ, 20 for DHA-PQP, and 1 for AL, whose 10-day
#'   protection is already less variable than `Tr` itself. A drug mixture is
#'   moment-matched as a mixture. `1` is a single exponential stage, which for
#'   `Ph_c` leaks protection early between monthly SMC rounds. The count is fixed
#'   at the seed's drug mix: a first-line switch moves the chain's mean, not its
#'   shape.
#' @param atol,rtol,step_size_max dust2 ODE-solver controls. The defaults
#'   (`1e-8`, `1e-8`, `1`) preserve the flat equilibrium exactly; for long dynamic
#'   projections `rtol = 1e-6` with `step_size_max = 10` runs about 1.4–1.7x
#'   faster on seasonal projections, and barely faster on aseasonal ones (there
#'   the daily output grid, not the tolerance, sets the step count), with
#'   negligible effect on aggregate outputs. Keep `atol` at `1e-8`: the individual
#'   prophylaxis chain stages hold occupancies of order `1e-6`, which a looser
#'   absolute tolerance lets dip below zero.
#' @param odin_file optional path to an odin source to compile instead of the
#'   generator built into the package (development use; needs 'odin2' and a C++
#'   toolchain).
#' @return a `blink_ode_tuning` list.
#' @seealso [run_simulation_ode()]
#' @examples
#' ode_tuning()$atol
#' # a faster long projection
#' ode_tuning(rtol = 1e-6, step_size_max = 10)$step_size_max
#' @export
ode_tuning <- function(age_lower = default_age_lower(),
                       n_eir = 10L, n_foim = 10L, n_eip = 20L,
                       n_ph = NULL, n_phc = NULL,
                       atol = 1e-8, rtol = 1e-8, step_size_max = 1,
                       odin_file = NULL) {
  pos_int <- function(x, nm, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(invisible())
    if (length(x) != 1L || is.na(x) || !is.numeric(x) || x < 1 || x != round(x)) {
      stop(sprintf("`%s` must be a single positive whole number.", nm), call. = FALSE)
    }
  }
  pos_num <- function(x, nm) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) || x <= 0) {
      stop(sprintf("`%s` must be a single positive number.", nm), call. = FALSE)
    }
  }
  ## age_lower is validated here rather than in build_inputs() because a bad grid
  ## does not fail loudly downstream -- it silently re-bins every output column
  if (!is.numeric(age_lower) || length(age_lower) < 2L || anyNA(age_lower) ||
      !all(is.finite(age_lower)) || age_lower[1] != 0 || is.unsorted(age_lower, strictly = TRUE)) {
    stop("`age_lower` must be a finite numeric vector of at least two age-group ",
         "lower edges in years, starting at 0 and strictly increasing.", call. = FALSE)
  }
  for (nm in c("n_eir", "n_foim", "n_eip")) pos_int(get(nm), nm)
  for (nm in c("n_ph", "n_phc")) pos_int(get(nm), nm, allow_null = TRUE)
  for (nm in c("atol", "rtol", "step_size_max")) pos_num(get(nm), nm)
  if (!is.null(odin_file) && (length(odin_file) != 1L || !is.character(odin_file))) {
    stop("`odin_file` must be NULL or a single path.", call. = FALSE)
  }
  structure(
    list(age_lower = age_lower, n_eir = n_eir, n_foim = n_foim, n_eip = n_eip,
         n_ph = n_ph, n_phc = n_phc, atol = atol, rtol = rtol,
         step_size_max = step_size_max, odin_file = odin_file),
    class = "blink_ode_tuning")
}

#' Coerce a partial named list to a full tuning object.
#'
#' Anything the caller left out keeps its default, and an unrecognised name is an
#' error rather than a silently ignored field, because a mistyped `rtol` that quietly
#' did nothing would look like the tolerance simply not mattering.
#' @noRd
as_ode_tuning <- function(x) {
  if (inherits(x, "blink_ode_tuning")) return(x)
  if (is.null(x)) return(ode_tuning())
  if (!is.list(x)) {
    stop("`tuning` must be an ode_tuning() object or a named list of its fields.",
         call. = FALSE)
  }
  known <- names(formals(ode_tuning))
  nm <- names(x)
  if (length(x) && (is.null(nm) || !all(nzchar(nm)))) {
    stop("every element of `tuning` must be named. Fields: ",
         paste(known, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(nm, known)
  if (length(bad)) {
    stop("unknown `tuning` field(s): ", paste(bad, collapse = ", "),
         ".\nFields are: ", paste(known, collapse = ", "), call. = FALSE)
  }
  do.call(ode_tuning, x)
}

## Arguments that were on run_simulation_ode()'s signature before it was cut back
## to malariasimulation::run_simulation()'s. Caught by name and reported with the
## call that replaces them: R's own "unused argument" error names the argument
## but not where it went, and every one of these has moved somewhere specific.
MOVED_TO_TUNING <- c("age_lower", "n_eir", "n_foim", "n_eip", "n_ph", "n_phc",
                     "atol", "rtol", "step_size_max", "odin_file")

#' @noRd
.reject_moved_args <- function(...) {
  if (...length() == 0L) return(invisible(NULL))
  nm <- names(list(...))
  if (is.null(nm)) nm <- rep("", ...length())
  if ("init_EIR" %in% nm) {
    stop("`init_EIR` is no longer an argument to run_simulation_ode(); it is read ",
         "from the parameter list, which is where malariasimulation keeps it:\n",
         "  parameters <- malariasimulation::set_equilibrium(parameters, init_EIR = 20)\n",
         "  run_simulation_ode(timesteps, parameters)", call. = FALSE)
  }
  moved <- intersect(nm, MOVED_TO_TUNING)
  if (length(moved)) {
    stop(paste0("`", moved, "`", collapse = ", "),
         if (length(moved) > 1) " are now fields " else " is now a field ",
         "of `tuning`:\n  run_simulation_ode(timesteps, parameters, ",
         "tuning = ode_tuning(", paste0(moved, " = ...", collapse = ", "), "))",
         call. = FALSE)
  }
  stop("unused argument(s): ", paste(nm[nzchar(nm)], collapse = ", "), call. = FALSE)
}
