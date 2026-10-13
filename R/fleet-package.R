#' fleet: mean-field counterpart to malariasimulation
#'
#' A deterministic mean-field counterpart to the \pkg{malariasimulation}
#' individual-based model (IBM) of *Plasmodium falciparum* and *P. vivax*
#' malaria, advanced one day at a time on the IBM's daily clock and built on
#' \pkg{odin2} and \pkg{dust2}. It reproduces the human model, structured by age
#' and biting heterogeneity and, for *P. vivax*, by hypnozoite batch, coupled to
#' the compartmental mosquito model, with the same intervention modules. It
#' takes an unmodified `malariasimulation::get_parameters()` list, is seeded at
#' the \pkg{malariaEquilibrium} (or \pkg{malariaEquilibriumVivax}) solution, and
#' returns `malariasimulation`'s output table, which \pkg{postie} post-processes
#' as it would an IBM run. Run time does not depend on population size.
#'
#' @section Functions:
#' \itemize{
#'   \item [run_simulation_ode()] runs the model.
#'   \item [ode_tuning()] holds the discretisation settings, passed as `tuning`.
#'   \item [default_age_lower()] gives the default age grid.
#' }
#' The `ode` in two of the names is kept for compatibility: the model is a daily
#' update, not an ODE.
#'
#' `vignette("fleet")` is a worked introduction, and `vignette("using")` covers
#' where the mean field departs from the IBM and what to do about it. Agreement
#' with the IBM is measured in a separate project,
#' [fleetcheck](https://pwinskill.github.io/fleetcheck/).
#'
#' @keywords internal
#' @useDynLib fleet, .registration = TRUE
"_PACKAGE"
