# Parameter translation: malariasimulation parameter list -> malariaEquilibrium
# ("eq") parameter names. Replicated from malariasimulation R/compatibility.R
# (back_translations) so this package is self-contained and version-robust.
# The resulting eq_params supply BOTH the human equilibrium seed
# (via malariaEquilibrium) AND the human ODE rate constants used by the odin model.

.inverse_param <- function(name, new_name) {
  function(params) list(new_name, 1 / params[[name]])
}
.mean_param <- function(new_name, name, weights) {
  function(params) list(new_name, stats::weighted.mean(params[[name]], params[[weights]]))
}
.product_param <- function(new_name, name_1, name_2) {
  function(params) list(new_name, params[[name_1]] * params[[name_2]])
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

# ---- P. vivax ---------------------------------------------------------------
# Replicated from malariasimulation R/compatibility.R (vivax_translations), and
# deliberately keyed the SAME way round as that table -- White/equilibrium name
# on the left, malariasimulation name on the right -- so the two can be diffed
# line by line. That is the opposite direction from .back_translations above,
# and the reason is that the vivax map is ADDITIVE rather than a rename: the
# vivax equilibrium reads some fields under their malariasimulation names
# (sigma_squared, n_heterogeneity_groups) and others under White's, from the
# same list. Renaming would break the first kind.
#
# The one entry that is not a plain alias: `phi1` in a malariasimulation vivax
# list is the RATIO phi_D_min/phi_D_max, not the level, so White's phi_D_min is
# the product of the two. fleet's odin already uses the ratio form directly --
# phi0 * (phi1 + (1 - phi1) / (...)) -- so this conversion is needed only for
# the equilibrium call, not for the model.
.vivax_translations <- list(
  mean_age = "average_age",
  rho_age  = "rho",
  age_0    = "a0",
  N_het    = "n_heterogeneity_groups",

  bb = "b",          # mosquito -> human transmission probability (constant; no IB)

  c_PCR = "cu",      # human -> mosquito infectivity, PCR-detectable
  c_LM  = "ca",      # human -> mosquito infectivity, LM-detectable (constant for vivax)
  c_D   = "cd",
  c_T   = "ct",

  d_E = "de",
  r_D = .inverse_param("dd", "r_D"),
  r_T = .inverse_param("dt", "r_T"),

  r_par  = .inverse_param("ra", "r_par"),   # anti-parasite immunity decay
  r_clin = .inverse_param("rc", "r_clin"),  # clinical immunity decay

  mu_M = .mean_param("mu_M", "mum", "species_proportions"),
  Q0   = .mean_param("Q0", "Q0", "species_proportions"),
  blood_meal_rates = .mean_param("blood_meal_rates", "blood_meal_rates",
                                 "species_proportions"),
  tau_M = "dem",

  ff      = "f",       # relapse rate per batch
  gamma_L = "gammal",  # hypnozoite batch clearance rate
  K_max   = "kmax",    # maximum hypnozoite batches

  u_par      = "ua",
  phi_LM_max = "philm_max",
  phi_LM_min = "philm_min",
  A_LM_50pc  = "alm50",
  K_LM       = "klm",
  u_clin     = "uc",
  phi_D_max  = "phi0",
  phi_D_min  = .product_param("phi_D_min", "phi0", "phi1"),
  A_D_50pc   = "ic0",
  K_D        = "kc",

  A_d_PCR_50pc = "apcr50",
  K_d_PCR      = "kpcr",
  d_PCR_max    = "dpcr_max",
  d_PCR_min    = "dpcr_min",
  d_LM         = "da",

  P_MI = "pcm",
  d_MI = "rm"
)

#' Add White-format aliases to a malariasimulation vivax parameter list.
#'
#' Additive, not a rename: the returned list is the input plus the equilibrium's
#' expected names, because `malariaEquilibriumVivax::vivax_equilibrium()` reads
#' both conventions from the same object.
#'
#' @param params a list from malariasimulation::get_parameters(parasite = "vivax").
#' @return the same list with the White-format aliases added.
#' @noRd
translate_vivax_parameters <- function(params) {
  translated <- params
  for (nm in names(.vivax_translations)) {
    tr <- .vivax_translations[[nm]]
    if (is.character(tr)) {
      translated[[nm]] <- params[[tr]]
    } else if (is.function(tr)) {
      translated[[nm]] <- tr(params)[[2]]
    }
  }
  translated
}

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
