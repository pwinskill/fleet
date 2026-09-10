# Build the odin2 input list (constants + equilibrium initial conditions) from
# a malariasimulation parameter list and a target adult init_EIR.

#' Default (graded) age grid
#'
#' Lower edges (in years) of the age groups. Fine in infancy, where immunity and
#' maternal dynamics move fast, and coarse in adulthood. Band aggregation in the
#' outputs weights each age group by the exact fraction of its own width that
#' falls inside the band, so a band edge landing inside a group (as it can on a
#' coarse custom grid) apportions that group between the two bands instead of
#' handing it whole to one of them; the default grid places edges at 2, 5, 10 and
#' 15 years, so the usual rendering bands fall on group boundaries exactly and
#' every weight is 0 or 1.
#' @param max_age oldest age-group lower edge, in **years** (absorbing top
#'   group). Must be a single finite number `>= 20`: the grid is graded up to a
#'   5-yearly section starting at 15, so there is no room for an absorbing top
#'   group below 20.
#' @return numeric vector of age-group lower edges in years.
#' @examples
#' default_age_lower()
#' # coarser top of the grid
#' default_age_lower(max_age = 60)
#' @export
default_age_lower <- function(max_age = 80) {
  # Below 20 the 5-yearly section seq(15, max_age - 5, by = 5) runs backwards and
  # dies with base R's opaque "wrong sign in 'by' argument", which says nothing
  # about age grids. Reject it here, in the units the caller is thinking in.
  if (length(max_age) != 1L || !is.numeric(max_age) || !is.finite(max_age) ||
      max_age < 20) {
    stop("`max_age` must be a single finite number >= 20 (years). The graded grid ",
         "runs 5-yearly from 15, so an absorbing top group below 20 years leaves ",
         "no room for it.", call. = FALSE)
  }
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
#' aging recursion, prophylaxis being an Erlang chain of `n_ph` stages.
#'
#' `malariaEquilibrium::human_equilibrium_no_het` computes the prophylaxis inflow
#' as `bP <- rT*bT + r[i-1]*P[i-1]/betaP`, in which the `rT*bT` term is not divided by
#' `betaP`. That makes its P/S columns not the fixed point of the clean-flux
#' prophylaxis ODE this model implements, so seeding directly from it leaves a
#' small (ft>0) transient. We reuse its FOI/phi/prop/r columns but re-solve the
#' disease block with the corrected `bP = (rT*bT + r[i-1]*P[i-1]) / betaP`,
#' generalised to a chain: stage 1 is fed by `rT*T`, stage m by `n_ph*rP*P[m-1]`,
#' and every stage ages at `r` and is left at `n_ph*rP` (chain mean `1/rP`).
#' @param FOI per-age equilibrium force of infection (length n_age).
#' @param phi per-age probability of clinical disease on infection.
#' @param prop per-age equilibrium population fraction (reference age structure).
#' @param r per-age aging rate (1/band width; 0 in the absorbing top group).
#' @param eta birth/death rate (1/average_age).
#' @param rA,rD,rU,rT recovery rates out of the A, D, U and Tr compartments.
#' @param rP prophylaxis exit rate: 1 / the chain's mean (drug-linked) duration.
#' @param ft effective treated fraction (coverage x drug efficacy).
#' @param n_ph number of prophylaxis chain stages.
#' @return list(S, T, D, A, U, P) of per-age equilibrium population fractions; `P`
#'   is an n_age x n_ph matrix of stage occupancies.
#' @noRd
solve_disease_block <- function(FOI, phi, prop, r, eta, rA, rD, rU, rT, rP, ft, n_ph = 1L) {
  n <- length(FOI); k <- as.integer(n_ph); rPk <- rP * k
  S <- Tc <- D <- A <- U <- numeric(n)
  P <- matrix(0, n, k)
  for (i in seq_len(n)) {
    re <- r[i] + eta
    betaT <- rT + re; betaD <- rD + re
    betaA <- FOI[i] * phi[i] + rA + re
    betaU <- FOI[i] + rU + re; betaP <- rPk + re
    aT <- ft * phi[i] * FOI[i] / betaT
    aD <- (1 - ft) * phi[i] * FOI[i] / betaD
    if (i == 1) {
      bT <- 0; bD <- 0
    } else {
      bT <- r[i - 1] * Tc[i - 1] / betaT
      bD <- r[i - 1] * D[i - 1] / betaD
    }
    # prophylaxis chain: each stage's occupancy is its inflow / betaP, split into
    # the part proportional to the at-risk pool Y (aP) and the rest (bP): stage 1
    # is fed by rT*T (corrected form, see @description) plus aging-in, stage m by
    # rPk*P[m-1] plus aging-in.
    aP <- bP <- numeric(k)
    aging_in <- function(m) if (i == 1) 0 else r[i - 1] * P[i - 1, m]
    aP[1] <- rT * aT / betaP
    bP[1] <- (rT * bT + aging_in(1)) / betaP
    for (m in seq_len(k)[-1]) {
      aP[m] <- rPk * aP[m - 1] / betaP
      bP[m] <- (rPk * bP[m - 1] + aging_in(m)) / betaP
    }
    Y <- (prop[i] - (bT + bD + sum(bP))) / (1 + aT + aD + sum(aP))
    Tc[i] <- aT * Y + bT
    D[i]  <- aD * Y + bD
    P[i, ] <- aP * Y + bP
    rA_in <- if (i == 1) 0 else r[i - 1] * A[i - 1]
    rU_in <- if (i == 1) 0 else r[i - 1] * U[i - 1]
    A[i] <- (rA_in + (1 - phi[i]) * Y * FOI[i] + rD * D[i]) /
      (betaA + (1 - phi[i]) * FOI[i])
    U[i] <- (rU_in + rA * A[i]) / betaU
    S[i] <- Y - A[i] - U[i]
  }
  list(S = S, T = Tc, D = D, A = A, U = U, P = P)
}

#' malariasimulation's mosquito sizing for a target EIR, replicated exactly.
#'
#' `set_equilibrium()` (compatibility.R) solves the human equilibrium under the
#' DEFAULT exponential age structure on its own 0.1-year grid (`EQUILIBRIUM_AGES`),
#' stores the resulting FOIM as `init_foim`, and `equilibrium_total_M()`
#' (mosquito_biology.R) converts the requested EIR into an adult-mosquito density
#' with a species-weighted `mum`. Custom demography enters nowhere. Reproduced here
#' rather than read from `parameters$total_M` so it holds for any `init_EIR` passed
#' at run time.
#' @param eqp_ibm translated equilibrium parameters with the IBM's raw (unconverted)
#'   sojourn rates.
#' @param ft_raw total clinical-treatment coverage at t = 1 (the IBM passes raw
#'   coverage here, not coverage x efficacy).
#' @noRd
ibm_total_M <- function(p, init_EIR, eqp_ibm, ft_raw) {
  n_het <- if (is.null(p$n_heterogeneity_groups)) 5L else p$n_heterogeneity_groups
  eq <- malariaEquilibrium::human_equilibrium(
    EIR = init_EIR, ft = ft_raw, p = eqp_ibm, age = 0:999 / 10,
    h = malariaEquilibrium::gq_normal(n_het)
  )
  mum <- stats::weighted.mean(p$mum, p$species_proportions)
  lifetime <- eq$FOIM * exp(-mum * p$dem) / (eq$FOIM + mum)
  (init_EIR * p$human_population / 365) /
    sum(p$species_proportions * p$blood_meal_rates * p$Q0 * lifetime)
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
  if (length(p$pev_epi_timesteps) > 0 && isTRUE(p$pev_epi_seasonal_boosters)) warning(
    "seasonal_boosters is approximated as a fixed days-since-primary schedule ",
    "(not aligned to the transmission season).", call. = FALSE)
  invisible(NULL)
}

#' Which frozen `eq_params` entries no longer agree with the live parameter list.
#'
#' `set_equilibrium()` stores its back-translation as `parameters$eq_params` (even
#' when passed `eq_params = NULL`), and `build_inputs()` merges that OVER
#' `translate_parameters()`, so the stored copy wins for every shared constant. On
#' the ordinary path the two are bit-identical, so this returns nothing; it returns
#' names only when the live list has moved underneath the frozen copy.
#'
#' Compared with a relative tolerance so that a value which merely took a different
#' arithmetic route to the same number (1/x re-inverted, a weighted mean recomputed)
#' cannot warn. Non-scalar and non-finite entries are skipped rather than guessed at.
#' @param live `translate_parameters(p)`, the current translation.
#' @param stored `as.list(p$eq_params)`, the frozen copy.
#' @param tol relative tolerance for "the same number".
#' @return character vector of offenders, labelled `eq_name (ibm_name)` where the
#'   malariasimulation name that feeds the constant can be identified; empty when
#'   the two agree.
#' @noRd
.stale_eq_params <- function(live, stored, tol = 1e-8) {
  nms <- intersect(names(live), names(stored))
  differs <- vapply(nms, function(k) {
    a <- live[[k]]; b <- stored[[k]]
    if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L) return(FALSE)
    if (!is.finite(a) || !is.finite(b)) return(FALSE)
    abs(a - b) > tol * max(1, abs(a), abs(b))
  }, logical(1))
  bad <- nms[differs]
  if (!length(bad)) return(character(0))
  # Only now (never on the clean path) pay for the reverse name map: users edit
  # `du`, not `rU`, so the warning has to name the field they actually typed.
  eq2ibm <- list()
  for (ibm in names(.back_translations)) {
    tr <- .back_translations[[ibm]]
    eq <- if (is.character(tr)) tr else environment(tr)[["new_name"]]
    if (!is.character(eq) || length(eq) != 1L) next
    eq2ibm[[eq]] <- c(eq2ibm[[eq]], ibm)
  }
  vapply(bad, function(k) {
    ibm <- eq2ibm[[k]]
    if (is.null(ibm)) k else sprintf("%s (eq `%s`)", paste(ibm, collapse = "/"), k)
  }, character(1), USE.NAMES = FALSE)
}

#' Build odin2 inputs.
#' @param parameters a malariasimulation::get_parameters() list (falciparum).
#' @param init_EIR target adult EIR (infectious bites per adult per year).
#' @param age_lower age-group lower edges in years.
#' @param n_eir,n_foim,n_eip Erlang-chain stage counts for the EIR lag, FOIM lag
#'   and mosquito EIP respectively.
#' @param n_ph,n_phc Erlang-chain stage counts for the post-treatment and the
#'   chemoprevention prophylaxis compartments; NULL picks the count that matches
#'   the variance of the drug's Weibull protection (capped at 20; see
#'   `erlang_stages()`), fixed at the seed's drug mix.
#' @param timesteps simulation horizon in days. Used only to size the
#'   intervention time-series grids (vector control, PEV, TBV, carrying capacity),
#'   which are built to extend past `timesteps` so the ODE stepper never
#'   extrapolates the interpolation tables.
#' @return list with `pars` (odin parameter list) and `meta` (grid metadata).
#' @noRd
build_inputs <- function(parameters, init_EIR, age_lower = default_age_lower(),
                         n_eir = 10L, n_foim = 10L, n_eip = 20L, n_ph = NULL, n_phc = NULL,
                         timesteps = 3650) {
  p <- parameters
  if (!is.null(p$parasite) && p$parasite == "vivax") {
    stop("blink supports P. falciparum only (parasite = 'vivax' not supported)")
  }
  # A single finite positive NUMBER. `is.na() || <= 0` alone let two values through
  # to die far downstream in messages that never mention init_EIR: Inf (reaches
  # `if (max(K0) > 0)` as NaN) and the string "20" (compared as a string, "20" <= 0
  # is FALSE, then dies at `EIR * zeta[j]`).
  if (length(init_EIR) != 1 || !is.numeric(init_EIR) || !is.finite(init_EIR) ||
      init_EIR <= 0) {
    stop("init_EIR must be a single positive number (bites per adult per year).")
  }
  check_unsupported(p)
  eqp <- translate_parameters(p)
  # honour a user-supplied custom equilibrium parameter set (set_equilibrium's
  # eq_params) explicitly, rather than relying on it round-tripping through the
  # IBM parameter names.
  if (!is.null(p$eq_params)) {
    eqp_stored <- as.list(p$eq_params)
    # eq_params wins, which means set_equilibrium() FREEZES every translated
    # biological constant at the value it had when it was called: an edit to the
    # live list afterwards (p$du <- 10) is a silent no-op. That contract stands --
    # set_equilibrium() is meant to be the last call on the list -- but a violation
    # of it must not be silent, so name the parameters whose live translation no
    # longer agrees with the frozen copy.
    stale <- .stale_eq_params(eqp, eqp_stored)
    eqp <- utils::modifyList(eqp, eqp_stored)
    if (length(stale)) {
      warning("parameters$eq_params (stored by set_equilibrium()) overrides the live ",
              "parameter list, so the following no longer has any effect: ",
              paste(stale, collapse = ", "), ". set_equilibrium() freezes the ",
              "translated biological constants at the values they had when it was ",
              "called, so call it LAST, after every other set_*() and manual edit. ",
              "(If you passed a custom `eq_params` to set_equilibrium(), this is ",
              "expected and can be ignored.)", call. = FALSE)
    }
  }
  eqp_ibm <- eqp                           # the IBM's own (unconverted) rates, for ibm_total_M()
  # Disease-progression rates: malariasimulation advances states once per whole day
  # with exit probability rate_to_prob(1/d) = 1 - exp(-1/d) (competing_hazards.R:78,
  # utils.R:118), so the REALISED mean dwell is 1/(1 - exp(-1/d)) -- e.g. 5.517 d for
  # dd = dt = 5, not 5. Use that exit rate so blink's dwell matches the IBM's. Must be
  # applied AFTER the eq_params merge above, or set_equilibrium's stored rates silently
  # overwrite it. Not applied to rP: ms models prophylaxis as a hazard multiplier, not
  # a compartment, so there is no per-day census of it.
  for (nm in c("rA", "rD", "rU", "rT")) {
    if (!is.null(eqp[[nm]]) && is.finite(eqp[[nm]]) && eqp[[nm]] > 0) {
      eqp[[nm]] <- 1 - exp(-eqp[[nm]])
    }
  }
  ft <- get_ft(p, 1)                       # baseline treatment coverage at t=1
  dser <- drug_mix_series(p, eqp)          # time-varying drug_eff/cT/rP; seed = t=0 blend
  eqp[["cT"]] <- dser$seed$cT              # drug-linked treated infectivity for the seed
  rP <- dser$seed$rP                       # drug-linked prophylaxis rate for the seed
  chp <- chemoprevention_prophylaxis(p, eqp)   # chemoprevention prophylaxis: chain rate + stages
  rP_c <- chp$rate
  # Erlang stage counts of the two prophylaxis chains: matched to the Weibull shape
  # of the drug (mix) unless overridden. 1 = the old single exponential compartment.
  chk_stages <- function(k, nm) {
    if (length(k) != 1 || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k))
      stop("`", nm, "` must be a single whole number >= 1", if (nm %in% c("n_ph", "n_phc"))
        " (or NULL for the default)", ".", call. = FALSE)
    as.integer(k)
  }
  n_eir <- chk_stages(n_eir, "n_eir"); n_foim <- chk_stages(n_foim, "n_foim")
  n_eip <- chk_stages(n_eip, "n_eip")
  n_ph  <- if (is.null(n_ph))  dser$seed$n_ph else chk_stages(n_ph, "n_ph")
  n_phc <- if (is.null(n_phc)) chp$n_stages   else chk_stages(n_phc, "n_phc")
  if (!is.null(p$hold_init_EIR) && !(is.logical(p$hold_init_EIR) &&
                                     length(p$hold_init_EIR) == 1 && !is.na(p$hold_init_EIR)))
    stop("parameters$hold_init_EIR must be TRUE or FALSE.", call. = FALSE)
  # The two mean-field fidelity knobs, validated here rather than passed straight
  # through to dust2. Both went unchecked, and both silently accept nonsense:
  # `bite_dedup` is a 0/1 switch that the odin model uses as a blend weight, so 2
  # and -1 both run and both extrapolate outside either documented endpoint; and
  # `acquired_immunity_offset = 50` runs, changing clinical incidence 12-fold,
  # while -5 makes (IB + offset)^kb non-finite and dies in the stepper.
  bite_dedup <- if (is.null(p$bite_dedup)) 1 else p$bite_dedup
  if (is.logical(bite_dedup) && length(bite_dedup) == 1L && !is.na(bite_dedup)) {
    bite_dedup <- as.numeric(bite_dedup)   # FALSE is the natural spelling of "0"
  }
  if (length(bite_dedup) != 1 || !is.numeric(bite_dedup) || !is.finite(bite_dedup) ||
      !(bite_dedup == 0 || bite_dedup == 1))
    stop("parameters$bite_dedup must be 0 or 1 (TRUE/FALSE also accepted): 1 ",
         "reproduces the IBM's per-timestep bite deduplication, 0 the linear ",
         "b*EPS form. It is a switch, not a dial.", call. = FALSE)
  acq_offset <- if (is.null(p$acquired_immunity_offset)) 0 else p$acquired_immunity_offset
  if (length(acq_offset) != 1 || !is.numeric(acq_offset) || !is.finite(acq_offset) ||
      acq_offset < 0 || acq_offset > 1)
    stop("parameters$acquired_immunity_offset must be a single number in [0, 1]: 0 ",
         "(the default) is the better mean-field match, 0.5 reproduces the IBM's ",
         "literal per-individual +0.5 in the b/phi/theta Hill functions.", call. = FALSE)
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
    # A Gauss-Hermite rule needs enough nodes to be a quadrature at all. Below 3 it
    # is not one: at n = 1 the single node sits at zeta = exp(-s2/2) = 0.434, so
    # every human in the model experiences 43% of the EIR the output column reports,
    # and nothing said so. Unvalidated, n = 0 and NA both died in gq_normal() or in
    # `if` with messages that never mention heterogeneity.
    if (length(n_het) != 1 || !is.numeric(n_het) || is.na(n_het) ||
        !is.finite(n_het) || n_het != round(n_het) || n_het < 3) {
      stop("parameters$n_heterogeneity_groups must be a single whole number >= 3: ",
           "fewer Gauss-Hermite nodes than that do not integrate the log-normal ",
           "biting distribution, they just move every human off it. For a genuine ",
           "no-heterogeneity run set parameters$enable_heterogeneity = FALSE, which ",
           "collapses to one stratum at zeta = 1 exactly.", call. = FALSE)
    }
    n_het <- as.integer(n_het)
    gq <- malariaEquilibrium::gq_normal(n_het)
    s2 <- eqp[["s2"]]
    zeta <- exp(gq$nodes * sqrt(s2) - s2 / 2)
    het_wt <- gq$weights
    # The nodes are used exactly as malariaEquilibrium emits them. The continuous
    # log-normal has E[zeta] = 1 by construction; the finite quadrature of it does
    # not (sum(w*zeta) is 0.99972 at the default 5 nodes, 0.97506 at 3), and
    # neither does a finite IBM sample. Do NOT renormalise to force sum(w*zeta) = 1:
    # malariasimulation::set_equilibrium() seeds through
    # malariaEquilibrium::human_equilibrium(h = gq_normal(n)) with these nodes
    # untouched, and the IBM's individuals draw zeta from the continuous
    # distribution, so that invariant is one neither library imposes. Rescaling
    # here would be a mean-matched fudge that moves blink off the mechanism it
    # replicates.
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
    # Age used to look up each model group's death rate. The open-ended top group is
    # represented by its lower bound in age_mid, which a right-closed bin would assign
    # to the band BELOW it (80 y -> the 60-80 rate on the default grid, so the over-80s
    # died at the wrong rate and the group held ~2.4x too many people). Everyone in that
    # group is older than its lower bound, so look it up one day above.
    age_lookup <- age_mid
    age_lookup[n_age] <- age_days[n_age] + 1
    if (max(age_lookup) > max(p$deathrate_agegroups)) {
      warning("The oldest model age group (", round(max(age_mid) / 365), "y) exceeds ",
              "the top set_demography age group (", round(max(p$deathrate_agegroups) / 365),
              "y); ages above it use the top death rate here, whereas the IBM removes ",
              "them. Raise default_age_lower(max_age=) or extend deathrate_agegroups ",
              "to match.", call. = FALSE)
    }
    # right-closed bins (edge_{g-1}, edge_g], matching ms .bincode(age, c(0, agegroups));
    # ages above the top edge cap to the last rate here (the warned divergence from the
    # IBM, which removes them).
    grp <- pmin(pmax(findInterval(age_lookup, c(0, p$deathrate_agegroups),
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
  # Maternal immunity is inherited from mothers aged [20, 21) years:
  # malariasimulation selects them with trunc(age / 365) == 20
  # (mortality_processes.R:42). Pick the model group CONTAINING age 20, not the
  # nearest midpoint: on the default grid 20 y is a band EDGE, so the midpoints
  # 17.5 y and 22.5 y are exactly equidistant and which.min() silently took the
  # FIRST, sourcing ICM/IVM from 15-20 year-olds and running them ~9-15% low.
  age20 <- max(1L, findInterval(20 * 365, age_days))
  mask20 <- numeric(n_age); mask20[age20] <- 1

  dm <- eqp[["dm"]]; dvm <- eqp[["dvm"]]
  icm_factor <- numeric(n_age); ivm_factor <- numeric(n_age)
  for (i in seq_len(n_age - 1)) {
    w <- age_days[i + 1] - age_days[i]
    icm_factor[i] <- dm / w * (exp(-age_days[i] / dm) - exp(-age_days[i + 1] / dm))
    ivm_factor[i] <- dvm / w * (exp(-age_days[i] / dvm) - exp(-age_days[i + 1] / dvm))
  }

  ## mosquito parameters (per species)
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
    # incubation survival exactly as malariasimulation (adult_mosquito_eqs.cpp:
    # incubation_survival = exp(-mu * tau)); the EIP chain is loss-free, so no
    # chain-survival compensation is needed and this stays correct as mum varies.
    g_s[s] <- exp(-mum_v[s] * dem)
  }

  ## human equilibrium seed (per het node) on the model age grid at adult EIR `EIR`,
  ## plus the adult-mosquito density at which sum_s a_s*Im_s = EIR/365 exactly under
  ## THIS age structure (the EIP chain is loss-free; survival g_s applies at exit).
  seed_human <- function(EIR) {
    S0 <- D0 <- A0 <- U0 <- Tr0 <- matrix(0, n_age, n_het)
    Ph0 <- array(0, c(n_age, n_het, n_ph)); Phc0 <- array(0, c(n_age, n_het, n_phc))
    IB_init <- ICA_init <- ID_init <- IVA_init <- matrix(0, n_age, n_het)
    cA0 <- matrix(0, n_age, n_het)
    for (j in seq_len(n_het)) {
      eqn <- malariaEquilibrium::human_equilibrium_no_het(
        EIR = EIR * zeta[j], ft = ft_seed, p = eqp, age = age_lower
      )
      w <- het_wt[j]
      db <- solve_disease_block(
        eqn[, "FOI"], eqn[, "phi"], eqn[, "prop"], eqn[, "r"],
        eqp[["eta"]], eqp[["rA"]], eqp[["rD"]], eqp[["rU"]], eqp[["rT"]], rP, ft_seed, n_ph
      )
      S0[, j]  <- w * resc * db$S
      Tr0[, j] <- w * resc * db$T
      D0[, j]  <- w * resc * db$D
      A0[, j]  <- w * resc * db$A
      U0[, j]  <- w * resc * db$U
      Ph0[, j, ] <- w * resc * db$P                 # [n_age, n_ph] stage occupancies
      IB_init[, j]  <- eqn[, "IB"]
      ICA_init[, j] <- eqn[, "ICA"]
      ID_init[, j]  <- eqn[, "ID"]
      IVA_init[, j] <- eqn[, "IVA"]
      cA0[, j] <- eqn[, "cA"]
    }
    ## equilibrium infectivity sum (Xf0) matching the odin inf_sum formula
    inf0 <- eqp[["cD"]] * D0 + cA0 * A0 + eqp[["cU"]] * U0 + eqp[["cT"]] * Tr0
    Xf0 <- sum(outer(psi, zeta) * inf0) / mean_psi
    ## per-species FOIM consistent with THIS model's seeded human infectivity
    foim0 <- a_spp * Xf0
    denom <- sum(a_spp * species_prop * foim0 * g_s / (foim0 + mum_v))
    list(S0 = S0, D0 = D0, A0 = A0, U0 = U0, Tr0 = Tr0, Ph0 = Ph0, Phc0 = Phc0,
         IB_init = IB_init, ICA_init = ICA_init, ID_init = ID_init, IVA_init = IVA_init,
         Xf0 = Xf0, foim0 = foim0, total_M = (EIR / 365) * hp / denom)
  }

  ## Which EIR to seed at. Under the default demography blink's own sizing is
  ## malariasimulation's formula evaluated at blink's own equilibrium (its age
  ## grid, whole-day rates, drug-linked cT and ft*eff), which lands within ~1-2% of
  ## the IBM's total_M, so the seed is at init_EIR. Under a custom demography the two conventions part
  ## ways: set_equilibrium() sizes total_M from the equilibrium under the DEFAULT
  ## exponential age structure (compatibility.R set_equilibrium ->
  ## equilibrium_total_M; custom mortality never enters), and the IBM then drifts to
  ## whatever transmission that mosquito density supports under the custom
  ## mortality. Replicate that: take the IBM's total_M and find the EIR at which
  ## blink's equilibrium under the custom age structure has exactly that density --
  ## a fixed point, so no burn-in is needed and the same parameter list realises the
  ## same transmission in both models. parameters$hold_init_EIR = TRUE keeps the
  ## previous behaviour (init_EIR is the EIR blink realises).
  eir_seed <- init_EIR
  total_M_ibm <- NA_real_
  sub_threshold <- FALSE
  if (isTRUE(p$custom_demography) && !isTRUE(p$hold_init_EIR)) {
    total_M_ibm <- ibm_total_M(p, init_EIR, eqp_ibm, ft)
    # blink's total_M(EIR) is increasing in EIR and tends to the transmission
    # threshold M_crit as EIR -> 0, so a root exists iff M_crit < total_M_ibm:
    # test at a near-zero EIR rather than at an arbitrary fraction of init_EIR.
    f <- function(lE) seed_human(exp(lE))$total_M - total_M_ibm
    lo <- log(init_EIR * 1e-6); f_lo <- f(lo)
    if (f_lo > 0) {
      # The IBM's mosquito density is below blink's transmission threshold under
      # this demography: transmission is not sustainable, and the IBM decays to
      # elimination from its seed. Keep the IBM's density (so blink decays too) and
      # start the humans at a near-zero EIR.
      warning("Under this custom demography the mosquito population that ",
              "set_equilibrium() implies is below blink's transmission threshold (",
              signif(total_M_ibm, 4), " vs ", signif(f_lo + total_M_ibm, 4), " adults); ",
              "seeding the mosquitoes at the IBM's density and the humans at a near-zero ",
              "EIR, so transmission decays as it would in the IBM. Set ",
              "parameters$hold_init_EIR = TRUE to seed at init_EIR instead.", call. = FALSE)
      eir_seed <- exp(lo); sub_threshold <- TRUE
    } else {
      eir_seed <- exp(stats::uniroot(f, c(lo, log(init_EIR * 50)), f.lower = f_lo,
                                     extendInt = "upX", tol = 1e-10)$root)
    }
  }
  hs <- seed_human(eir_seed)
  S0 <- hs$S0; D0 <- hs$D0; A0 <- hs$A0; U0 <- hs$U0; Tr0 <- hs$Tr0
  Ph0 <- hs$Ph0; Phc0 <- hs$Phc0
  IB_init <- hs$IB_init; ICA_init <- hs$ICA_init; ID_init <- hs$ID_init; IVA_init <- hs$IVA_init
  Xf0 <- hs$Xf0; foim0 <- hs$foim0
  total_M <- if (sub_threshold) total_M_ibm else hs$total_M
  Xe0 <- eir_seed / 365   # equilibrium adult EIR per day

  ## equilibrium FOIM (het-integrated) at the seed, for reference
  eq_full <- malariaEquilibrium::human_equilibrium(
    EIR = eir_seed, ft = ft_seed, p = eqp, age = age_lower, h = gq
  )
  init_foim <- eq_full$FOIM

  ME0 <- ML0 <- MP0 <- Sm0 <- Im0 <- Em_inc0 <- K0 <- numeric(n_spp)
  Xi0 <- matrix(0, n_spp, n_eip)
  for (s in seq_len(n_spp)) {
    m_s <- species_prop[s] * total_M
    counts <- initial_mosquito_counts(p, s, foim0[s], m_s) / hp
    ME0[s] <- counts[["E"]]; ML0[s] <- counts[["L"]]; MP0[s] <- counts[["P"]]
    Sm0[s] <- counts[["Sm"]]                       # = m_s*mum/(foim+mum)/hp
    # the EIP chain is a loss-free delay of Sm*foim, so at equilibrium every stage
    # holds the (constant) inflow; survival exp(-mum*dem) is applied at the exit.
    Xi0[s, ] <- Sm0[s] * foim0[s]
    Em_inc0[s] <- counts[["Pm"]]                   # incubating stock (ms E state)
    Im0[s] <- counts[["Im"]]                       # = m_s*foim/(foim+mum)*exp(-mum*dem)/hp
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
    n_eir = n_eir, n_foim = n_foim, n_eip = n_eip, n_ph = n_ph, n_phc = n_phc,
    r_age = r_age, psi = psi,
    n_mut = ncol(mu_age_z), mu_age_t = mu_age_t, mu_age_z = mu_age_z,
    age_mid = age_mid, mask20 = mask20,
    icm_factor = icm_factor, ivm_factor = ivm_factor, mean_psi = mean_psi,
    zeta = zeta, het_wt = het_wt,
    rA = eqp[["rA"]], rD = eqp[["rD"]], rU = eqp[["rU"]],
    rT = eqp[["rT"]], rP_c = rP_c, rT_slow = res$rT_slow,
    # slow-clearance fraction at t = 0, used to split the Tr seed between the fast
    # and slow treated compartments (the IBM assigns each treated individual to one
    # or the other by a Bernoulli draw, so Tr is a two-component mixture).
    spc0 = if (length(res$spc)) res$spc[1] else 0,
    d_ib = eqp[["db"]], d_ica = eqp[["dc"]], d_id = eqp[["dd"]], d_iva = eqp[["dv"]],
    # integer refractory windows: ms tests (timestep - last_boosted) >= u with an
    # integer timestep, so the realised wait is ceil(u) days -> u_eff = ceil(u) - 1.
    ub_eff = ceiling(eqp[["ub"]]) - 1, uc_eff = ceiling(eqp[["uc"]]) - 1,
    ud_eff = ceiling(eqp[["ud"]]) - 1, uv_eff = ceiling(eqp[["uv"]]) - 1,
    # Acquired-immunity offset in the b/phi/theta Hill calls. The IBM adds +0.5 per
    # individual; empirically (A/B vs the ms 3.0.0 IBM ensemble mean) offset 0 is the
    # better mean-field match, so default 0. Set parameters$acquired_immunity_offset =
    # 0.5 to reproduce the IBM's literal per-individual Hill functions.
    acq_offset = acq_offset,
    # Reproduce the IBM's per-timestep bite deduplication (saturating hazard) by default;
    # parameters$bite_dedup = 0 restores the linear b*EPS form, for which the
    # malariaEquilibrium seed is an exact fixed point (used by the flat-equilibrium tests).
    bite_dedup = bite_dedup,
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
    ME0 = ME0, ML0 = ML0, MP0 = MP0, Sm0 = Sm0, Xi0 = Xi0, Em_inc0 = Em_inc0, Im0 = Im0,
    Xe0 = Xe0, Xf0 = Xf0
  )

  list(
    pars = pars,
    meta = list(
      n_age = n_age, n_het = n_het, n_spp = n_spp, age_lower = age_lower,
      age_mid = age_mid, age_lo = age_days, age_hi = age_days + width,
      prop = prop, init_EIR = init_EIR, eir_seed = eir_seed, init_foim = init_foim,
      total_M = total_M, total_M_ibm = total_M_ibm, ft = ft, eq_full = eq_full,
      human_population = hp, parameters = p, eqp = eqp, n_eip = n_eip,
      n_ph = n_ph, n_phc = n_phc
    )
  )
}
