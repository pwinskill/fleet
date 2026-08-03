#' blink: a mean-field (ODE) twin of malariasimulation
#'
#' A fast, deterministic mean-field (ODE) counterpart to the
#' \pkg{malariasimulation} individual-based model of *Plasmodium falciparum*
#' malaria, built on the \pkg{odin2} / \pkg{dust2} stack. It reproduces the
#' age- and biting-heterogeneity-structured human model (states S/D/A/U/Tr plus
#' treatment and chemoprevention prophylaxis, and the six immunity functions)
#' coupled to the compartmental mosquito model, and honours the same intervention
#' modules. It accepts an unmodified `malariasimulation::get_parameters()` list,
#' is seeded at the \pkg{malariaEquilibrium} fixed point, and returns
#' malariasimulation-style outputs that feed straight into \pkg{postie}. A
#' multi-decade run completes in seconds, independent of population size.
#'
#' @section Getting started:
#' \itemize{
#'   \item [run_simulation_ode()] — run the model; returns a wide, daily count
#'     table (P. falciparum only).
#'   \item [get_epi_outputs()] — post-process a run into postie `rates` (long)
#'     and `prevalence` (wide).
#'   \item [default_age_lower()] — the default graded age grid.
#' }
#' See `vignette("blink")` for a worked tour, and the README's *Mean-field
#' approximations* section for where and how much the ODE departs from the IBM.
#'
#' @keywords internal
"_PACKAGE"
