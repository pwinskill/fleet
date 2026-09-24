# Parameter translation: malariasimulation parameter list -> malariaEquilibrium
# ("eq") parameter names. Replicated from malariasimulation R/compatibility.R
# (back_translations) so this package is self-contained and version-robust.
# The resulting eq_params supply BOTH the human equilibrium seed
# (via malariaEquilibrium) AND the human rate constants used by the odin model.

.inverse_param <- function(name, new_name) {
  function(params) list(new_name, 1 / params[[name]])
}
.mean_param <- function(new_name, name, weights) {
  function(params) list(new_name, stats::weighted.mean(params[[name]], params[[weights]]))
}

# IBM name -> eq name (character) or a function(params) -> list(eq_name, value)
.back_translations <- list(
  average_age = .inverse_param("average_age", "eta"),
  rho = "rho",
  a0  = "a0",
  da  = .inverse_param("da", "rA"),
  dd  = .inverse_param("dd", "rD"),
  du  = .inverse_param("du", "rU"),
  dt  = .inverse_param("dt", "rT"),
  de  = "dE",
  cd  = "cD",
  cu  = "cU",
  ct  = "cT",
  d1  = "d1",
  rid = "dd",
  id0 = "ID0",
  kd  = "kd",
  ud  = "ud",
  ad  = "ad0",
  gammad = "gd",
  b0  = "b0",
  b1  = "b1",
  rb  = "db",
  ib0 = "IB0",
  kb  = "kb",
  ub  = "ub",
  phi0 = "phi0",
  phi1 = "phi1",
  rc  = "dc",
  ic0 = "IC0",
  kc  = "kc",
  uc  = "uc",
  rm  = "dm",
  mum = .mean_param("mu", "mum", "species_proportions"),
  sigma_squared = "s2",
  fd0 = "fd0",
  gamma1 = "g_inf",
  pcm = "PM",
  dem = "tau",
  Q0  = .mean_param("Q0", "Q0", "species_proportions"),
  blood_meal_rates = .mean_param("f", "blood_meal_rates", "species_proportions"),
  delay_gam = "tl",
  uv  = "uv",
  rva = "dv",
  pvm = "PVM",
  rvm = "dvm",
  fv0 = "fv0",
  av  = "av",
  gammav = "gammav",
  theta0 = "theta0",
  theta1 = "theta1",
  iv0 = "IV0",
  kv  = "kv"
)

#' Translate a malariasimulation parameter list to malariaEquilibrium eq params.
#'
#' @param params a list produced by malariasimulation::get_parameters() (+ set_*).
#' @return a named list of malariaEquilibrium parameters (a plain list; starts
#'   from malariaEquilibrium::load_parameter_set() defaults and overrides).
#' @noRd
translate_parameters <- function(params) {
  translated <- as.list(malariaEquilibrium::load_parameter_set())
  for (name in names(params)) {
    tr <- .back_translations[[name]]
    if (is.null(tr)) next
    if (is.character(tr)) {
      translated[[tr]] <- params[[name]]
    } else if (is.function(tr)) {
      kv <- tr(params)
      translated[[kv[[1]]]] <- kv[[2]]
    }
  }
  translated
}
