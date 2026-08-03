# Build the odin2 input list (constants + equilibrium initial conditions) from
# a malariasimulation parameter list and a target adult init_EIR.

#' Default (graded) age grid
#'
#' Lower edges (in years) of the age groups. Fine in infancy, where immunity and
#' maternal dynamics move fast, and coarse in adulthood. Band aggregation in the
#' outputs assigns each age group to a band by its midpoint, so bands whose edges
#' fall inside a group (on a coarse custom grid) may be mis-binned; the default
#' grid places edges at 2, 5, 10 and 15 years.
#' @param max_age oldest age-group lower edge, in years (absorbing top group).
#' @return numeric vector of age-group lower edges in years.
#' @examples
#' default_age_lower()
#' # coarser top of the grid
#' default_age_lower(max_age = 60)
#' @export
default_age_lower <- function(max_age = 80) {
  c(
    seq(0, 1 - 1 / 12, by = 1 / 12),   # monthly 0-1y   (12)
    seq(1, 5 - 0.25, by = 0.25),       # quarterly 1-5y (16)
    seq(5, 15 - 1, by = 1),            # yearly 5-15y    (10)
    seq(15, max_age - 5, by = 5),      # 5-yearly        (13)
    max_age                            # absorbing top    (1)
  )
}

#' Total clinical treatment coverage active at timestep t.
#' @noRd
get_ft <- function(p, t = 1) {
  cov <- p$clinical_treatment_coverages
  if (is.null(cov) || length(cov) == 0) return(0)
  total <- 0
  for (d in seq_along(cov)) {
    ts <- p$clinical_treatment_timesteps[[d]]
    idx <- which(ts <= t)
    if (length(idx)) total <- total + cov[[d]][max(idx)]
  }
  min(total, 1)                          # cap as treatment_series() does (parity)
}

#' Re-solve the S/T/D/A/U/P equilibrium per age with the CORRECTED prophylaxis
#' aging recursion.
#'
#' `malariaEquilibrium::human_equilibrium_no_het` computes the prophylaxis inflow
#' as `bP <- rT*bT + r[i-1]*P[i-1]/betaP` — the `rT*bT` term is not divided by
#' `betaP`. That makes its P/S columns not the fixed point of the clean-flux
#' prophylaxis ODE this model implements, so seeding directly from it leaves a
#' small (ft>0) transient. We reuse its FOI/phi/prop/r columns but re-solve the
#' disease block with the corrected `bP = (rT*bT + r[i-1]*P[i-1]) / betaP`.
#' @param FOI per-age equilibrium force of infection (length n_age).
#' @param phi per-age probability of clinical disease on infection.
#' @param prop per-age equilibrium population fraction (reference age structure).
#' @param r per-age aging rate (1/band width; 0 in the absorbing top group).
#' @param eta birth/death rate (1/average_age).
#' @param rA,rD,rU,rT recovery rates out of the A, D, U and Tr compartments.
#' @param rP prophylaxis exit rate (drug-linked mean duration).
#' @param ft effective treated fraction (coverage x drug efficacy).
#' @return list(S, T, D, A, U, P) of per-age equilibrium population fractions.
#' @noRd
solve_disease_block <- function(FOI, phi, prop, r, eta, rA, rD, rU, rT, rP, ft) {
  n <- length(FOI)
  S <- Tc <- D <- A <- U <- P <- numeric(n)
  for (i in seq_len(n)) {
    re <- r[i] + eta
    betaT <- rT + re; betaD <- rD + re
    betaA <- FOI[i] * phi[i] + rA + re
    betaU <- FOI[i] + rU + re; betaP <- rP + re
    aT <- ft * phi[i] * FOI[i] / betaT
    aP <- rT * aT / betaP
    aD <- (1 - ft) * phi[i] * FOI[i] / betaD
    if (i == 1) {
      bT <- 0; bD <- 0; bP <- 0
    } else {
      bT <- r[i - 1] * Tc[i - 1] / betaT
      bD <- r[i - 1] * D[i - 1] / betaD
      bP <- (rT * bT + r[i - 1] * P[i - 1]) / betaP   # corrected (see @description)
    }
    Y <- (prop[i] - (bT + bD + bP)) / (1 + aT + aD + aP)
    Tc[i] <- aT * Y + bT
    D[i]  <- aD * Y + bD
    P[i]  <- aP * Y + bP
    rA_in <- if (i == 1) 0 else r[i - 1] * A[i - 1]
    rU_in <- if (i == 1) 0 else r[i - 1] * U[i - 1]
    A[i] <- (rA_in + (1 - phi[i]) * Y * FOI[i] + rD * D[i]) /
      (betaA + (1 - phi[i]) * FOI[i])
    U[i] <- (rU_in + rA * A[i]) / betaU
    S[i] <- Y - A[i] - U[i]
  }
  list(S = S, T = Tc, D = D, A = A, U = U, P = P)
}

#' Guard for malariasimulation features not represented in the mean-field model.
#'
#' Every intervention module is now modelled, including custom demography,
#' time-varying EPI PEV coverage and coverage-weighted multi-drug resistance. The
#' resistance arms that malariasimulation actually implements (early-treatment-
#' failure, slow-parasite-clearance) are both modelled; its partner-drug and
#' late-failure arms are rejected upstream by set_antimalarial_resistance, so
#' cannot reach us. Individual-level effects with no mean-field analogue
#' (intervention correlation, per-person antibody variation) are documented
#' approximations rather than errors. Kept as an extension point.
#' @noRd
check_unsupported <- function(p) {
  # Surface the remaining mean-field approximation so it is never silent: seasonal
  # PEV boosters are scheduled at a fixed delay from the primary series, not aligned
  # to the transmission season. (Clinical-drug mixes, incl. first-line switches, are
  # now modelled time-varyingly via drug_mix_series().)
  if (isTRUE(p$pev) && isTRUE(p$pev_epi_seasonal_boosters)) warning(
    "seasonal_boosters is approximated as a fixed days-since-primary schedule ",
    "(not aligned to the transmission season).", call. = FALSE)
  invisible(NULL)
}

#' Build odin2 inputs.
#' @param parameters a malariasimulation::get_parameters() list (falciparum).
#' @param init_EIR target adult EIR (infectious bites per adult per year).
#' @param age_lower age-group lower edges in years.
#' @param n_eir,n_foim,n_eip Erlang-chain stage counts for the EIR lag, FOIM lag
#'   and mosquito EIP respectively.
#' @param timesteps simulation horizon in days. Used only to size the
#'   intervention time-series grids (vector control, PEV, TBV, carrying capacity),
#'   which are built to extend past `timesteps` so the ODE stepper never
#'   extrapolates the interpolation tables.
#' @return list with `pars` (odin parameter list) and `meta` (grid metadata).
#' @noRd
build_inputs <- function(parameters, init_EIR, age_lower = default_age_lower(),
                         n_eir = 10L, n_foim = 10L, n_eip = 20L, timesteps = 3650) {
  p <- parameters
  if (!is.null(p$parasite) && p$parasite == "vivax") {
    stop("blink supports P. falciparum only (parasite = 'vivax' not supported)")
  }
  if (length(init_EIR) != 1 || is.na(init_EIR) || init_EIR <= 0) {
    stop("init_EIR must be a single positive number (bites per adult per year).")
  }
  check_unsupported(p)
  eqp <- translate_parameters(p)
  # honour a user-supplied custom equilibrium parameter set (set_equilibrium's
  # eq_params) explicitly, rather than relying on it round-tripping through the
  # IBM parameter names.
  if (!is.null(p$eq_params)) eqp <- utils::modifyList(eqp, as.list(p$eq_params))
  ft <- get_ft(p, 1)                       # baseline treatment coverage at t=1
  dser <- drug_mix_series(p, eqp)          # time-varying drug_eff/cT/rP; seed = t=0 blend
  eqp[["cT"]] <- dser$seed$cT              # drug-linked treated infectivity for the seed
  rP <- dser$seed$rP                       # drug-linked prophylaxis rate for the seed
  rP_c <- chemoprevention_prophylaxis_rate(p, eqp)  # chemoprevention prophylaxis rate
  ft_seed <- ft * dser$seed$drug_eff       # effective treated fraction = coverage * efficacy

  ## heterogeneity nodes
  n_het <- p$n_heterogeneity_groups
  if (is.null(n_het)) n_het <- 5
  if (isFALSE(p$enable_heterogeneity)) {
    zeta <- 1
    het_wt <- 1
    n_het <- 1
    gq <- list(nodes = 0, weights = 1)
  } else {
    gq <- malariaEquilibrium::gq_normal(n_het)
    s2 <- eqp[["s2"]]
    zeta <- exp(gq$nodes * sqrt(s2) - s2 / 2)
    het_wt <- gq$weights
  }

  ## age grid + demography (replicates malariaEquilibrium::human_equilibrium_no_het)
  n_age <- length(age_lower)
  age_days <- age_lower * 365
  width <- c(diff(age_days), Inf)
  r_age <- 1 / width
  r_age[n_age] <- 0
  eta <- eqp[["eta"]]
  age_mid <- c((age_days[-n_age] + age_days[-1]) / 2, age_days[n_age])
  psi <- 1 - eqp[["rho"]] * exp(-age_mid / eqp[["a0"]])
  # constant-hazard age structure = the malariaEquilibrium seed reference
  prop_ce <- numeric(n_age)
  for (i in seq_len(n_age)) {
    if (i == 1) prop_ce[i] <- eta / (r_age[1] + eta)
    else prop_ce[i] <- prop_ce[i - 1] * r_age[i - 1] / (r_age[i] + eta)
  }
  # age-specific mortality. Custom demography (set_demography) is time-varying:
  # mu_age is an interpolated [n_age, n_mut] series over deathrate_timesteps. The
  # baseline (t = 0) row seeds the equilibrium age structure below. Constant
  # demography is a flat two-knot series at 1/average_age.
  if (isTRUE(p$custom_demography)) {
    if (!all(is.finite(p$deathrates)) || any(p$deathrates <= 0)) {
      stop("set_demography deathrates must be finite and > 0.", call. = FALSE)
    }
    if (max(age_mid) > max(p$deathrate_agegroups)) {
      warning("The oldest model age group (", round(max(age_mid) / 365), "y) exceeds ",
              "the top set_demography age group (", round(max(p$deathrate_agegroups) / 365),
              "y); ages above it use the top death rate here, whereas the IBM removes ",
              "them. Raise default_age_lower(max_age=) or extend deathrate_agegroups ",
              "to match.", call. = FALSE)
    }
    # right-closed bins (edge_{g-1}, edge_g], matching ms .bincode(age, c(0, agegroups));
    # ages above the top edge cap to the last rate here (the warned divergence from the
    # IBM, which removes them).
    grp <- pmin(pmax(findInterval(age_mid, c(0, p$deathrate_agegroups),
                                  left.open = TRUE), 1L), ncol(p$deathrates))
    mu_age <- p$deathrates[1, grp]                    # baseline (seed) hazard by model age
    mu_age_z <- t(p$deathrates[, grp, drop = FALSE])  # [n_age, n_mut] (time last)
    mu_age_t <- p$deathrate_timesteps
  } else {
    mu_age <- rep(eta, n_age)
    mu_age_z <- matrix(eta, n_age, 2)
    mu_age_t <- c(0, timesteps + 365)
  }
  if (ncol(mu_age_z) < 2) {                  # interpolate needs >= 2 knots
    mu_age_z <- cbind(mu_age_z, mu_age_z[, 1])
    mu_age_t <- c(mu_age_t[1], timesteps + 365)
  }
  # equilibrium age structure under this mortality (McKendrick aging chain)
  prop <- numeric(n_age)
  prop[1] <- 1 / (r_age[1] + mu_age[1])
  if (n_age >= 2) for (i in 2:n_age) prop[i] <- prop[i - 1] * r_age[i - 1] / (r_age[i] + mu_age[i])
  prop <- prop / sum(prop)
  mean_psi <- sum(prop * psi)
  resc <- prop / prop_ce           # rescale the constant-eta seed to this age structure
  age20 <- which.min(abs(age_mid - 20 * 365))
  mask20 <- numeric(n_age); mask20[age20] <- 1

  dm <- eqp[["dm"]]; dvm <- eqp[["dvm"]]
  icm_factor <- numeric(n_age); ivm_factor <- numeric(n_age)
  for (i in seq_len(n_age - 1)) {
    w <- age_days[i + 1] - age_days[i]
    icm_factor[i] <- dm / w * (exp(-age_days[i] / dm) - exp(-age_days[i + 1] / dm))
    ivm_factor[i] <- dvm / w * (exp(-age_days[i] / dvm) - exp(-age_days[i + 1] / dvm))
  }

  ## human equilibrium seed (per het node), on the model age grid
  S0 <- D0 <- A0 <- U0 <- Tr0 <- Ph0 <- Phc0 <- matrix(0, n_age, n_het)
  IB_init <- ICA_init <- ID_init <- IVA_init <- matrix(0, n_age, n_het)
  cA0 <- matrix(0, n_age, n_het)
  for (j in seq_len(n_het)) {
    eqn <- malariaEquilibrium::human_equilibrium_no_het(
      EIR = init_EIR * zeta[j], ft = ft_seed, p = eqp, age = age_lower
    )
    w <- het_wt[j]
    db <- solve_disease_block(
      eqn[, "FOI"], eqn[, "phi"], eqn[, "prop"], eqn[, "r"],
      eqp[["eta"]], eqp[["rA"]], eqp[["rD"]], eqp[["rU"]], eqp[["rT"]], rP, ft_seed
    )
    S0[, j]  <- w * resc * db$S
    Tr0[, j] <- w * resc * db$T
    D0[, j]  <- w * resc * db$D
    A0[, j]  <- w * resc * db$A
    U0[, j]  <- w * resc * db$U
    Ph0[, j] <- w * resc * db$P
    IB_init[, j]  <- eqn[, "IB"]
    ICA_init[, j] <- eqn[, "ICA"]
    ID_init[, j]  <- eqn[, "ID"]
    IVA_init[, j] <- eqn[, "IVA"]
    cA0[, j] <- eqn[, "cA"]
  }
  ## equilibrium infectivity sum (Xf0) matching the odin inf_sum formula
  inf0 <- eqp[["cD"]] * D0 + cA0 * A0 + eqp[["cU"]] * U0 + eqp[["cT"]] * Tr0
  Xf0 <- sum(outer(psi, zeta) * inf0) / mean_psi
  Xe0 <- init_EIR / 365   # equilibrium adult EIR per day

  ## equilibrium FOIM (het-integrated) for mosquito seeding
  eq_full <- malariaEquilibrium::human_equilibrium(
    EIR = init_EIR, ft = ft_seed, p = eqp, age = age_lower, h = gq
  )
  init_foim <- eq_full$FOIM

  ## mosquito parameters + equilibrium seed (per-human densities)
  ## EIP is an Erlang chain of n_eip stages; per-stage rate reip = n_eip/dem.
  ## Chain survival g_s replaces exp(-mum*dem); total_M is matched to g so that
  ## sum_s a_s*Im_s = init_EIR/365 exactly.
  species_prop <- p$species_proportions
  n_spp <- length(species_prop)
  hp <- p$human_population
  dem <- p$dem
  reip <- n_eip / dem
  mum_v <- a_spp <- beta_eff <- g_s <- numeric(n_spp)
  for (s in seq_len(n_spp)) {
    fmr <- p$blood_meal_rates[[s]]
    mum_v[s] <- p$mum[[s]]
    beta_eff[s] <- eggs_laid(p$beta, mum_v[s], fmr)
    a_spp[s] <- p$Q0[[s]] * fmr
    g_s[s] <- (reip / (reip + mum_v[s]))^n_eip
  }
  ## per-species FOIM consistent with THIS model's seeded human infectivity
  foim0 <- a_spp * Xf0
  denom <- sum(a_spp * species_prop * foim0 * g_s / (foim0 + mum_v))
  total_M <- (init_EIR / 365) * hp / denom

  ME0 <- ML0 <- MP0 <- Sm0 <- Im0 <- K0 <- numeric(n_spp)
  Em0 <- matrix(0, n_spp, n_eip)
  for (s in seq_len(n_spp)) {
    m_s <- species_prop[s] * total_M
    counts <- initial_mosquito_counts(p, s, foim0[s], m_s) / hp
    ME0[s] <- counts[["E"]]; ML0[s] <- counts[["L"]]; MP0[s] <- counts[["P"]]
    Sm0[s] <- counts[["Sm"]]                       # = m_s*mum/(foim+mum)/hp
    lam <- reip; mm <- mum_v[s]
    Em0[s, 1] <- Sm0[s] * foim0[s] / (lam + mm)
    if (n_eip >= 2) for (k in 2:n_eip) Em0[s, k] <- Em0[s, k - 1] * lam / (lam + mm)
    Im0[s] <- lam * Em0[s, n_eip] / mm
    K0[s] <- calculate_carrying_capacity(p, m_s, s) / hp
  }
  # A species with proportion 0 (common in site files: one dominant vector, others
  # at 0) gives K0 = 0, so its larval carrying capacity Kcap = 0. Its larval count
  # is also 0, so the aquatic ODE evaluates 0 / Kcap = 0/0 = NaN and the stepper
  # crashes. Floor K0 to a negligible positive value: the species stays inert (0
  # mosquitoes, 0 biting) but the RHS is finite.
  if (max(K0) > 0) K0 <- pmax(K0, max(K0) * 1e-9)

  ## intervention time series (Phase 2). Each starts at its baseline at t = 0,
  ## so the equilibrium seed above is preserved; interventions act from their
  ## scheduled timesteps. (VC t=0 row equals the baseline a/mu/beta used above.)
  trt <- treatment_series(p)
  res <- resistance_series(p, eqp)
  vc <- vector_control_series(p, timesteps)
  pevs <- pev_series(p, age_mid, timesteps)
  tbvs <- tbv_series(p, age_mid, timesteps)
  ccs <- carrying_capacity_series(p, K0, timesteps)

  pars <- list(
    n_age = n_age, n_het = n_het, n_spp = n_spp,
    n_eir = n_eir, n_foim = n_foim, n_eip = n_eip,
    r_age = r_age, psi = psi,
    n_mut = ncol(mu_age_z), mu_age_t = mu_age_t, mu_age_z = mu_age_z,
    age_mid = age_mid, mask20 = mask20,
    icm_factor = icm_factor, ivm_factor = ivm_factor, mean_psi = mean_psi,
    zeta = zeta, het_wt = het_wt,
    rA = eqp[["rA"]], rD = eqp[["rD"]], rU = eqp[["rU"]],
    rT = eqp[["rT"]], rP_c = rP_c, rT_slow = res$rT_slow,
    d_ib = eqp[["db"]], d_ica = eqp[["dc"]], d_id = eqp[["dd"]], d_iva = eqp[["dv"]],
    ub = eqp[["ub"]], uc = eqp[["uc"]], ud = eqp[["ud"]], uv = eqp[["uv"]],
    # Acquired-immunity offset in the b/phi/theta Hill calls. The IBM adds +0.5 per
    # individual; empirically (A/B vs the ms 3.0.0 IBM ensemble mean) offset 0 is the
    # better mean-field match, so default 0. Set parameters$acquired_immunity_offset =
    # 0.5 to reproduce the IBM's literal per-individual Hill functions.
    acq_offset = if (!is.null(p$acquired_immunity_offset)) p$acquired_immunity_offset else 0,
    b0 = eqp[["b0"]], b1 = eqp[["b1"]], ib0 = eqp[["IB0"]], kb = eqp[["kb"]],
    phi0 = eqp[["phi0"]], phi1 = eqp[["phi1"]], ic0 = eqp[["IC0"]], kc = eqp[["kc"]],
    d1 = eqp[["d1"]], id0 = eqp[["ID0"]], kd = eqp[["kd"]],
    fd0 = eqp[["fd0"]], ad0 = eqp[["ad0"]], gd = eqp[["gd"]],
    theta0 = eqp[["theta0"]], theta1 = eqp[["theta1"]], iv0 = eqp[["IV0"]], kv = eqp[["kv"]],
    fv0 = eqp[["fv0"]], av = eqp[["av"]], gammav = eqp[["gammav"]],
    cD = eqp[["cD"]], cU = eqp[["cU"]], g_inf = eqp[["g_inf"]],
    PM = eqp[["PM"]], PVM = eqp[["PVM"]],
    n_ftt = length(trt$times), ft_times = trt$times, ft_vals = trt$vals,
    n_dmix = length(dser$times), dmix_times = dser$times,
    cT_vals = dser$cT, drug_eff_vals = dser$drug_eff, rP_vals = dser$rP,
    n_rest = length(res$times), res_times = res$times,
    etf_vals = res$etf, spc_vals = res$spc,
    del = p$del, dl = p$dl, dpl = p$dpl, me = p$me, ml = p$ml, mup = p$mup,
    mosq_gamma = p$gamma, beta_eff = beta_eff,
    n_cct = length(ccs$times), cc_times = ccs$times, K_vals = ccs$vals,
    n_vct = length(vc$times), vc_times = vc$times,
    a_vals = vc$a, mum_vals = vc$mum,
    n_pevt = length(pevs$times), pev_times = pevs$times, pev_vals = pevs$vals,
    n_tbvt = length(tbvs$times), tbv_times = tbvs$times,
    tbv_fU_vals = tbvs$fU, tbv_fA_vals = tbvs$fA,
    tbv_fD_vals = tbvs$fD, tbv_fT_vals = tbvs$fT,
    de = p$de, tl = p$delay_gam, dem = p$dem,
    S0 = S0, D0 = D0, A0 = A0, U0 = U0, Tr0 = Tr0, Ph0 = Ph0, Phc0 = Phc0,
    IB_init = IB_init, ICA_init = ICA_init, ID_init = ID_init, IVA_init = IVA_init,
    ME0 = ME0, ML0 = ML0, MP0 = MP0, Sm0 = Sm0, Em0 = Em0, Im0 = Im0,
    Xe0 = Xe0, Xf0 = Xf0
  )

  list(
    pars = pars,
    meta = list(
      n_age = n_age, n_het = n_het, n_spp = n_spp, age_lower = age_lower,
      age_mid = age_mid, age_lo = age_days, age_hi = age_days + width,
      prop = prop, init_EIR = init_EIR, init_foim = init_foim,
      total_M = total_M, ft = ft, eq_full = eq_full,
      human_population = hp, parameters = p, eqp = eqp, n_eip = n_eip
    )
  )
}
