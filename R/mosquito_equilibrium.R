# Mosquito equilibrium / biology helpers, replicated from malariasimulation
# R/mosquito_biology.R and src/mosquito_biology.cpp (compartmental path only).
# Griffin et al., "Modelling the impact of vector control interventions on
# Anopheles gambiae population dynamics".

# Per-capita daily oviposition rate (eggs_laid in src/mosquito_biology.cpp).
eggs_laid <- function(beta, mu, f) {
  eov <- beta / mu * (exp(mu / f) - 1)
  eov * mu * exp(-mu / f) / (1 - exp(-mu / f))
}

# 'omega' auxiliary for aquatic equilibrium (R/mosquito_biology.R:calculate_omega).
calculate_omega <- function(p, species) {
  sub_omega <- p$gamma * p$ml / p$me - (p$del / p$dl) +
    ((p$gamma - 1) * p$ml * p$del)
  mum <- p$mum[[species]]
  beta_eff <- eggs_laid(p$beta, mum, p$blood_meal_rates[[species]])
  -0.5 * sub_omega + sqrt(
    0.25 * sub_omega^2 +
      0.5 * p$gamma * beta_eff * p$ml * p$del /
        (p$me * mum * p$dl * (1 + p$dpl * p$mup))
  )
}

# Baseline larval carrying capacity K0 for m adults of a species.
calculate_carrying_capacity <- function(p, m, species) {
  omega <- calculate_omega(p, species)
  m * 2 * p$dl * p$mum[[species]] * (1 + p$dpl * p$mup) * p$gamma * (omega + 1) /
    (omega / (p$ml * p$del) - (1 / (p$ml * p$dl)) - 1)
}

# Equilibrium compartment counts c(E, L, P, Sm, Pm, Im) for m adults of a species.
initial_mosquito_counts <- function(p, species, foim, m) {
  omega <- calculate_omega(p, species)
  mum <- p$mum[[species]]
  n_E <- 2 * omega * mum * p$dl * (1 + p$dpl * p$mup) * m
  n_L <- 2 * mum * p$dl * (1 + p$dpl * p$mup) * m
  n_P <- 2 * p$dpl * mum * m
  n_Sm <- m * mum / (foim + mum)
  inc_surv <- exp(-mum * p$dem)
  n_Pm <- m * foim / (foim + mum) * (1 - inc_surv)
  n_Im <- m * foim / (foim + mum) * inc_surv
  c(E = n_E, L = n_L, P = n_P, Sm = n_Sm, Pm = n_Pm, Im = n_Im)
}

# (The former equilibrium_total_M() helper was removed: total_M is solved inline
# in build_inputs() to be consistent with the incubation survival it applies.)
