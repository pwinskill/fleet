# Post-processing to postie-format outputs. The data.frame returned by
# run_simulation_ode() is a wide, malariasimulation-style count table, so it can
# be passed straight to postie::get_rates() / postie::get_prevalence().

#' Post-process an ODE run into postie-format rates and prevalence.
#'
#' Returns a list with a `rates` table (long: one row per timestep x age band,
#' with `clinical`, `severe`, `mortality`, `yld`, `yll`, `dalys`, `person_days`)
#' and a `prevalence` table (wide: one `<diagnostic>_prevalence_<lo>_<hi>` column
#' per age band), matching `postie::get_rates()` / `postie::get_prevalence()`.
#'
#' @param x output of [run_simulation_ode()].
#' @param diagnostic prevalence diagnostic, "lm" or "pcr".
#' @param rates_args named list of extra arguments for [postie::get_rates()]
#'   (e.g. `scaler`, `treatment_scaler`, `life_expectancy`, `infer_ft`).
#' @param prevalence_args named list of extra arguments for
#'   [postie::get_prevalence()].
#' @param ... shared arguments accepted by both postie functions (e.g.
#'   `baseline_year`, `ages_as_years`).
#' @return list(rates, prevalence).
#' @examples
#' \dontrun{
#' p <- malariasimulation::get_parameters()
#' out <- run_simulation_ode(timesteps = 3650, parameters = p, init_EIR = 20)
#' epi <- get_epi_outputs(out, diagnostic = "lm")
#' head(epi$rates)
#' epi$prevalence$lm_prevalence_2_10
#' }
#' @export
get_epi_outputs <- function(x, diagnostic = "lm",
                            rates_args = list(), prevalence_args = list(), ...) {
  if (!requireNamespace("postie", quietly = TRUE)) {
    stop("Package 'postie' is required for get_epi_outputs(); install mrc-ide/postie.")
  }
  shared <- list(...)
  list(
    rates = do.call(postie::get_rates, c(list(x), rates_args, shared)),
    prevalence = do.call(
      postie::get_prevalence,
      c(list(x, diagnostic = diagnostic), prevalence_args, shared)
    )
  )
}
