# Discretisation settings for run_simulation_ode().
#
# These are gathered into one object rather than sitting on the run signature so
# that run_simulation_ode(timesteps, parameters, correlations) matches
# malariasimulation::run_simulation() argument for argument. Nothing
# epidemiological belongs here: every field is a numerical-approximation knob,
# and changing one should move the answer only by its own error term. Model
# parameters, the target EIR included, live on the parameter list.

#' Discretisation settings
#'
#' Numerical settings for [run_simulation_ode()], passed as its `tuning`
#' argument. Every default is the validated choice, and a run with the defaults
#' is the reference configuration. `run_simulation_ode()` also accepts a plain
#' named list of the fields to change: fields left out keep their defaults, and
#' an unknown field is an error.
#'
#' @param age_lower age-group lower edges in years. The default is
#'   [default_age_lower()], a 118-group grid log-spaced between pinned reporting
#'   ages and weighted towards the under-fives; its documentation gives the
#'   grid's accuracy and cost. Must be finite, start at 0 and increase strictly,
#'   with at least two edges, and must be in years: a top edge above 1000 is
#'   rejected as a grid given in days.
#'
#'   Each day every group loses its ageing and mortality fraction out of the
#'   same stock as the day's other exits, so the narrowest group must be wide
#'   enough to hold them together. Progression and clearance do not depend on
#'   transmission, and a grid with a group too narrow for them, a little over a
#'   day, is an error before the run. Infection is checked as the run goes: a
#'   run in which any compartment falls below zero stops with an error naming
#'   the day and the age group. That takes several consecutive groups a day or
#'   two wide at an EIR in the hundreds; the default grid's narrowest groups are
#'   16.6 days wide.
#' @param n_ph,n_phc stage counts of the post-treatment (`Ph`) and
#'   chemoprevention (`Ph_c`) prophylaxis chains. `NULL` (default) matches each
#'   chain to the drug's Weibull protection curve. A stage is left with a fixed
#'   probability each day, so it lasts a geometric number of days, and a chain
#'   matches the protection's mean `M` and variance `V` with `M^2 / (V + M)`
#'   stages: 10 for SP-AQ as chemoprevention, 9 for DHA-PQP. `Ph` follows the
#'   treated stage `Tr`, so its count matches the variance of the whole
#'   `Tr + Ph` sojourn and its mean is the protection left after `Tr`: 9 stages
#'   for SP-AQ, 10 for DHA-PQP, and 1 for AL, whose 10-day protection is already
#'   less variable than `Tr` itself. A drug mixture is moment-matched as a
#'   mixture.
#'
#'   A stage lasts at least a day and also loses its ageing and mortality
#'   fraction each day, so the default is limited to what fits the
#'   shortest-protecting drug on the schedule and capped at 20, and an explicit
#'   count that does not fit is an error. The count is fixed for the run, sized
#'   at the seed's drug mix, so a first-line switch moves the chain's mean, not
#'   its shape.
#' @param n_sub sub-steps per day of the mosquito model, integrated by an
#'   exponential midpoint rule (second order, and stable however stiff the
#'   larval equations become). Once a run is on its cycle the error falls
#'   fourfold with each doubling. At the default of 32, a seasonal run's daily
#'   EIR stays within 2.4e-4 of the converged solution at the seasonal trough
#'   (5.5e-6 of its peak), comparable to the per-step tolerance to which
#'   `malariasimulation` integrates its own mosquito model
#'   (`r_tol = a_tol = 1e-4`); 8 leaves it 3.6e-3 off, a slight shift in the
#'   phase of the seasonal cycle. Annual means agree to about 4e-6 at 32 and
#'   5e-5 at 8, and a run with no seasonality and no interventions is exact to
#'   1e-11 at any count. The mosquito model is a small part of the day, so 32
#'   costs 2 to 13% more than 8.
#' @param batch_bins under *P. vivax*, whether hypnozoite batch counts share
#'   levels. `malariasimulation` gives each person a count of hypnozoite
#'   batches, 0 to `kmax`, and `fleet` holds every vivax compartment and
#'   immunity stock at each level of that ladder, so its length accounts for
#'   most of a vivax run's cost. With `TRUE` (default), counts 0 to 3 and `kmax`
#'   have a level each and the counts between share levels in pairs, the people
#'   within a pair spread as a Poisson fitted to the levels either side. That
#'   gives eight levels for the default `kmax = 10`, which take a quarter less
#'   time than eleven and move no output of the `fleetcheck` vivax suite by more
#'   than 0.4%. `FALSE` gives every count its own level. Ignored under
#'   *P. falciparum*.
#' @param odin_file optional path to an odin source to compile instead of the
#'   generator built into the package: a discrete-time model with the built-in
#'   one's parameters and states. For development; needs 'odin2' and a C++
#'   toolchain. The compiled model is cached by path for the session, so an
#'   edit to the file needs a new session.
#' @param ... the ODE solver's controls (`atol`, `rtol`, `step_size_max`) and the
#'   Erlang lag stage counts (`n_eir`, `n_foim`, `n_eip`) are accepted with a
#'   warning and ignored: the model advances one day at a time, with the IBM's
#'   own delay lines, so there is no solver to control and no lag to
#'   discretise. Any other name is an error.
#' @return A list of class `fleet_ode_tuning`.
#' @seealso [run_simulation_ode()], [default_age_lower()]
#' @examples
#' ode_tuning()$n_sub
#' # a finer age grid, to check how much of a severe-incidence level is the grid's
#' length(ode_tuning(age_lower = default_age_lower(n_group = 417))$age_lower)
#' @export
ode_tuning <- function(age_lower = default_age_lower(), n_ph = NULL, n_phc = NULL,
                       n_sub = 32L, odin_file = NULL, ..., batch_bins = TRUE) {
  .retire_tuning_fields(...)
  pos_int <- function(x, nm, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(invisible())
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 1 || x != round(x) ||
        x > .Machine$integer.max) {
      stop(sprintf("`%s` must be a single positive whole number.", nm), call. = FALSE)
    }
  }
  ## age_lower is validated here rather than in build_inputs() because a bad grid
  ## does not fail loudly downstream -- it silently re-bins every output column
  if (!is.numeric(age_lower) || length(age_lower) < 2L || anyNA(age_lower) ||
      !all(is.finite(age_lower)) || age_lower[1] != 0 || is.unsorted(age_lower, strictly = TRUE)) {
    stop("`age_lower` must be a finite numeric vector of at least two age-group ",
         "lower edges in years, starting at 0 and strictly increasing.", call. = FALSE)
  }
  ## ...and in YEARS. A grid handed over in days (default_age_lower() * 365 is the
  ## easy slip, since every rendering band and every output tag is in days) passes
  ## every check above, then quietly puts the whole population in the first age
  ## group: n_age_730_3650 comes back 0 and the run looks merely uninteresting.
  ##
  ## The threshold only has to separate years from days, and the two scales differ
  ## by 365x, so it belongs an order of magnitude above any plausible lifespan
  ## rather than just above it. At 200 it rejected grids the package itself builds
  ## -- default_age_lower(max_age = 250) is accepted there (its own guard asks only
  ## for >= 20) and returns a perfectly usable absorbing top group -- while the
  ## smallest days-valued grid it needs to catch, default_age_lower(max_age = 20) *
  ## 365, tops out at 7300. 1000 years rejects every days grid with 7x to spare and
  ## accepts every years grid a caller could mean.
  if (max(age_lower) > 1000) {
    stop("`age_lower` must be in YEARS, but its top edge is ", signif(max(age_lower), 6),
         " -- not a plausible human age. A grid supplied in days (for example ",
         "default_age_lower() * 365) is the usual cause; divide by 365.", call. = FALSE)
  }
  for (nm in c("n_ph", "n_phc")) pos_int(get(nm), nm, allow_null = TRUE)
  pos_int(n_sub, "n_sub")
  if (!(is.logical(batch_bins) && length(batch_bins) == 1L && !is.na(batch_bins))) {
    stop("`batch_bins` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(odin_file) && (length(odin_file) != 1L || !is.character(odin_file))) {
    stop("`odin_file` must be NULL or a single path.", call. = FALSE)
  }
  structure(
    list(age_lower = age_lower, n_ph = n_ph, n_phc = n_phc, n_sub = as.integer(n_sub),
         batch_bins = batch_bins,
         odin_file = odin_file),
    class = "fleet_ode_tuning")
}

## Fields of the ODE solver fleet no longer has. A script written for it keeps
## running -- each is ignored, and says so -- rather than failing on a setting
## that no longer means anything.
RETIRED_TUNING <- c("atol", "rtol", "step_size_max", "n_eir", "n_foim", "n_eip")

#' @noRd
.retire_tuning_fields <- function(...) {
  if (...length() == 0L) return(invisible(NULL))
  nm <- names(list(...))
  if (is.null(nm) || !all(nzchar(nm))) {
    stop("every `tuning` field must be named. Fields: ",
         paste(setdiff(names(formals(ode_tuning)), "..."), collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(nm, RETIRED_TUNING)
  if (length(bad)) {
    stop("unknown `tuning` field(s): ", paste(bad, collapse = ", "),
         ".\nFields are: ", paste(setdiff(names(formals(ode_tuning)), "..."), collapse = ", "),
         call. = FALSE)
  }
  warning("`tuning` field(s) ", paste0("`", nm, "`", collapse = ", "), " ignored: fleet ",
          "advances one day at a time, as malariasimulation does, with the IBM's own ",
          "delay lines, so there is no ODE solver to control and no lag chain to ",
          "size. Drop ", if (length(nm) > 1) "them" else "it", ".", call. = FALSE)
  invisible(NULL)
}

#' Coerce a partial named list to a full tuning object.
#'
#' Anything the caller left out keeps its default, and an unrecognised name is an
#' error rather than a silently ignored field, because a mistyped `n_sub` that
#' quietly did nothing would look like the setting simply not mattering.
#' @noRd
as_ode_tuning <- function(x) {
  if (is.null(x)) return(ode_tuning())
  # A tuning object is re-validated rather than trusted: one saved by an older
  # fleet, or edited by hand, goes through the same checks as a fresh one.
  if (inherits(x, "fleet_ode_tuning")) x <- unclass(x)
  if (!is.list(x)) {
    stop("`tuning` must be an ode_tuning() object or a named list of its fields.",
         call. = FALSE)
  }
  known <- setdiff(names(formals(ode_tuning)), "...")
  nm <- names(x)
  if (length(x) && (is.null(nm) || !all(nzchar(nm)))) {
    stop("every element of `tuning` must be named. Fields: ",
         paste(known, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(nm, c(known, RETIRED_TUNING))
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
MOVED_TO_TUNING <- c("age_lower", "n_ph", "n_phc", "n_sub", "odin_file")

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
  retired <- intersect(nm, RETIRED_TUNING)
  if (length(retired)) {
    stop(paste0("`", retired, "`", collapse = ", "), " belonged to the ODE model, its ",
         "solver and its Erlang lag chains, and fleet now advances one day at a time ",
         "with the IBM's own delay lines: drop ", if (length(retired) > 1) "them." else "it.",
         call. = FALSE)
  }
  stop("unused argument(s): ", paste(nm[nzchar(nm)], collapse = ", "), call. = FALSE)
}
