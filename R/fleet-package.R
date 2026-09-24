#' fleet: a mean-field twin of malariasimulation
#'
#' A fast, deterministic mean-field counterpart to the \pkg{malariasimulation}
#' individual-based model of *Plasmodium falciparum* malaria, built on the
#' \pkg{odin2} / \pkg{dust2} stack and advanced one day at a time on the
#' individual-based model's own clock. It reproduces the age- and
#' biting-heterogeneity-structured human model (states S/D/A/U/Tr plus treatment
#' and chemoprevention prophylaxis, and the six immunity functions) coupled to
#' the compartmental mosquito model, and honours the same intervention modules.
#' It accepts an unmodified `malariasimulation::get_parameters()` list, is
#' seeded at the \pkg{malariaEquilibrium} solution, and returns
#' malariasimulation-style outputs that feed straight into \pkg{postie}. A
#' multi-decade run completes in under a second, independent of population
#' size.
#'
#' @section Getting started:
#' \itemize{
#'   \item [run_simulation_ode()]: run the model; returns a wide, daily count
#'     table (P. falciparum only).
#'   \item [default_age_lower()]: the default graded age grid.
#'   \item [ode_tuning()]: discretisation settings, passed as `tuning =`.
#' }
#' The count table is malariasimulation-shaped, so post-process it with
#' \pkg{postie} exactly as you would an IBM run: `postie::get_rates(out)` and
#' `postie::get_prevalence(out, diagnostic = "lm")`.
#' See `vignette("fleet")` for a worked tour and `vignette("using")` for where
#' and how much the mean field departs from the IBM and what to do about it. The
#' measured agreement behind those claims is a separate project,
#' [fleetcheck](https://pwinskill.github.io/fleetcheck/).
#'
#' @keywords internal
#' @useDynLib fleet, .registration = TRUE
"_PACKAGE"
