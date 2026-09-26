# Build the odin2 input list (constants + equilibrium initial conditions) from
# a malariasimulation parameter list and a target adult init_EIR.

#' Anchor ages, in years, that a default grid always puts a group edge on.
#'
#' The boundaries malaria burden is conventionally reported at. Pinning them
#' means the usual rendering bands fall on group boundaries exactly and every
#' band weight is 0 or 1, so an age profile is a sum of whole groups rather than
#' a re-apportionment of partial ones. 21 is there for maternal immunity: the IBM
#' draws a newborn's from a mother aged [20, 21) years (trunc(age / 365) == 20),
#' and edges at 20 and 21 give the mean field whole groups making up that year.
#' @noRd
AGE_ANCHORS <- c(0, 1, 2, 3, 5, 7, 10, 15, 20, 21, 30, 40, 60)

#' Default (graded) age grid
#'
#' Lower edges (in years) of the age groups: fine in childhood, where immunity and
#' maternal dynamics move fast, and coarse in adulthood. Group edges are pinned
#' at the conventional reporting boundaries -- 0, 1, 2, 3, 5, 7, 10, 15, 20, 30,
#' 40, 60 years -- and at 21, so that whole groups make up exactly the year the
#' IBM draws maternal immunity from, [20, 21). The groups between two anchors are
#' of equal width. The budget of groups is spread across the anchor intervals in
#' proportion to their width in `log(1 + age)`, weighted 1.25 below 5 years and
#' 0.5 from 21 years; each interval gets at least one group, and the remainder is
#' shared out largest-remainder.
#'
#' Log rather than linear width, because what the grid has to resolve is the
#' rise of immunity with age, and that is far closer to a function of log age
#' than of age. (Of `log(1 + age)` rather than `log(age)`, which is infinite
#' over the first year. The offset is in years, so the allocation is not
#' scale-free.) The grading matters more than the count: at 53 groups an
#' equal-width grid is four times as far from the converged solution as a
#' log-graded one (rms over the age profile), and a grid of fixed monthly /
#' quarterly / yearly / 5-yearly sections, which over-resolves infancy and
#' under-resolves everything above 15, 1.7 times as far. The weights move
#' groups to where the grid's error comes from. It is largest in clinical and
#' severe incidence among young children at high transmission, each age band's
#' error comes from the groups inside it, and the groups above 21 years move no
#' child's incidence measurably. Weighted, the grid matches the unweighted one's
#' accuracy with about a quarter fewer groups.
#'
#' The default, 118 groups, is the smallest at which every falciparum claim in
#' `fleetcheck` passes: the one that decides it is the clinical age profile's
#' 3-5 year band at EIR 120, which needs 14 groups between 3 and 5 years. It is
#' the smallest grid that passes, not the most accurate. The departure from
#' `fleet`'s own converged profile is fleet-low and halves with each doubling of
#' `n_group`, first order. On the default the largest departure of the EIR 20
#' clinical age profile from the converged one is 3.2%, in the oldest band,
#' which the weights coarsen, and the rms 1.6%; at 53 groups they are 6.1% and
#' 3.2%. Measured against the IBM median, severe incidence at EIR 120 is 2% low.
#' What refinement does **not** remove is the mean-field approximation itself --
#' one immunity value per stratum, where the IBM holds a spread of infection
#' histories at the same age -- so a persistent difference from the IBM is not
#' evidence that the grid is too coarse. `validations/age-grid/run.R` in
#' `fleetcheck` measures both.
#'
#' Band aggregation weights each age group by the exact fraction of its own width
#' that falls inside the band, so a band edge landing inside a group (as one
#' can, away from the anchors) apportions that group between the two bands
#' instead of handing it whole to one of them. The same is true of intervention
#' age targeting. Neither is snapped to the grid.
#'
#' @param max_age oldest age-group lower edge, in **years** (the absorbing top
#'   group). Must be a single finite number `>= 20`: below that the grid has too
#'   few anchors left to resolve the ages over which immunity develops.
#' @param n_group total number of age groups, the absorbing top group included.
#'   The default is 118. Changing it refines or coarsens the whole grid while
#'   keeping its shape, which is what a grid-convergence check wants; run time
#'   is roughly in proportion.
#' @return numeric vector of age-group lower edges in years.
#' @examples
#' default_age_lower()
#' # coarser top of the grid
#' default_age_lower(max_age = 60)
#' # under half the resolution everywhere, same shape, over twice as fast
#' length(default_age_lower(n_group = 53))
#' @export
default_age_lower <- function(max_age = 80, n_group = 118L) {
  if (length(max_age) != 1L || !is.numeric(max_age) || !is.finite(max_age) ||
      max_age < 20) {
    stop("`max_age` must be a single finite number >= 20 (years). Below that ",
         "the grid has too few anchor ages left to resolve the ages over which ",
         "immunity develops.", call. = FALSE)
  }
  # An anchor a hair below `max_age` leaves a final band of microseconds, whose
  # ageing rate is then thousands per day against a grid whose next-narrowest is
  # 52 days -- a stiff system from a typo, accepted by every check downstream.
  # Absorb an anchor that close instead.
  anchors <- c(AGE_ANCHORS[AGE_ANCHORS < max_age - 1 / 365], max_age)
  n_iv <- length(anchors) - 1L
  if (length(n_group) != 1L || !is.numeric(n_group) || !is.finite(n_group) ||
      n_group != round(n_group) || n_group > 1e6 || n_group < n_iv + 1L) {
    stop("`n_group` must be a single whole number between ", n_iv + 1L,
         " and 1e6 for max_age = ", max_age, ": one group per interval between ",
         "the anchor ages, plus the absorbing top group. A fractional value was ",
         "silently truncated before this check existed, and one above the ",
         "integer range died inside as.integer() saying nothing about age grids.",
         call. = FALSE)
  }
  # one group per interval, then the remainder shared out by weighted log width:
  # 1.25 below 5 years, 0.5 from 21 (see above). Largest remainder rather than
  # rounding, so the groups allocated always sum to the budget exactly --
  # rounding each share independently loses or gains one.
  lo <- anchors[-length(anchors)]; hi <- anchors[-1]
  lw <- (log1p(hi) - log1p(lo)) * ifelse(hi <= 5, 1.25, 1) * ifelse(lo >= 21, 0.5, 1)
  spare <- as.integer(n_group) - 1L - n_iv
  cnt <- rep(1L, n_iv)
  if (spare > 0L) {
    share <- lw / sum(lw) * spare
    cnt <- cnt + as.integer(floor(share))
    short <- spare - sum(as.integer(floor(share)))
    if (short > 0L) {
      # Rounded before ordering, because some intervals have mathematically
      # IDENTICAL weighted log widths -- log1p(2)-log1p(1) and log1p(5)-log1p(3)
      # are the same number, under the same weight -- and their computed
      # remainders differ only in the last bits. Ordering on that noise let a
      # different libm return a different grid, silently, for some n_group.
      # Rounding makes true ties exact so order()'s stable index tie-break
      # decides them, which favours the younger interval and is reproducible
      # everywhere.
      top <- order(round(share - floor(share), 12), seq_len(n_iv),
                   decreasing = c(TRUE, FALSE), method = "radix")[seq_len(short)]
      cnt[top] <- cnt[top] + 1L
    }
  }
  # equal width within each interval; `max_age` closes the last one and is the
  # lower edge of the absorbing group
  c(unlist(lapply(seq_len(n_iv), function(i)
    lo[i] + (hi[i] - lo[i]) * (seq_len(cnt[i]) - 1L) / cnt[i])), max_age)
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
#' aging recursion, prophylaxis being a chain of `n_ph` stages.
#'
#' `malariaEquilibrium::human_equilibrium_no_het` computes the prophylaxis inflow
#' as `bP <- rT*bT + r[i-1]*P[i-1]/betaP`, in which the `rT*bT` term is not divided by
#' `betaP`. That makes its P/S columns not the fixed point of the clean-flux
#' prophylaxis this model implements, so seeding directly from it leaves a
#' small (ft>0) transient. We reuse its FOI/phi/prop/r columns but re-solve the
#' disease block with the corrected `bP = (rT*bT + r[i-1]*P[i-1]) / betaP`,
#' generalised to a chain: stage 1 is fed by `rT*T`, stage m by `n_ph*rP*P[m-1]`,
#' and every stage ages at `r` and is left at `n_ph*rP` (chain mean `1/rP`). A
#' chain stage's fixed point is its inflow / (exit + ageing) on the daily clock
#' as in continuous time, so the chain is seeded at its own equilibrium.
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

#' The same for P. vivax: set_equilibrium() takes its FOIM from
#' malariaEquilibriumVivax, on the same 0.1-year grid, and equilibrium_total_M()
#' is otherwise parasite-agnostic.
#' @noRd
ibm_total_M_vivax <- function(p, init_EIR, ft_raw) {
  eq <- malariaEquilibriumVivax::vivax_equilibrium(
    EIR = init_EIR, ft = ft_raw, p = translate_vivax_parameters(p), age = 0:999 / 10)
  mum <- stats::weighted.mean(p$mum, p$species_proportions)
  lifetime <- eq$FOIM * exp(-mum * p$dem) / (eq$FOIM + mum)
  (init_EIR * p$human_population / 365) /
    sum(p$species_proportions * p$blood_meal_rates * p$Q0 * lifetime)
}

# The vivax block's inputs for a falciparum run: switched off, collapsed to one
# cell, and empty. One cell rather than zero because odin2 cannot have a
# zero-length dimension -- the same reason enable_heterogeneity collapses n_het
# to 1 rather than removing it. The rate constants only have to be finite and
# keep every denominator away from zero; an empty block never uses them.
inert_vivax_block <- function() {
  z <- array(0, c(1L, 1L, 1L))
  c(list(pf_on = 1, pv_on = 0, n_age_v = 1L, n_het_v = 1L, n_bat = 1L, n_hyp = 1L,
         n_phv = 1L, rc_on = 0, n_ra = 1L, n_rc = 1L,
         Sv0 = z, Dv0 = z, Av0 = z, Uv0 = z, Trv0 = z, JAv0 = z, JCv0 = z,
         KAv0 = z, KCv0 = z, spread_on = 0, n_q = 1L, qz = 0, qw = 1,
         Phv0 = array(0, c(1L, 1L, 1L, 1L))),
    as.list(stats::setNames(rep(1, length(VIVAX_SCALARS)), VIVAX_SCALARS)))
}

# Upper edge given to malariaEquilibriumVivax for fleet's open-ended top age
# group. Nobody is alive at 200, so the bin holds everyone over fleet's last edge.
VIVAX_TOP_AGE <- 200

# Within-cell immunity spread under vivax (see the note at mA in the odin
# model): the Gauss-Hermite rule the gamma is read at, and the switch.
# Five nodes hold the clinical curve (kc = 5.4) within ~5% of the exact gamma
# expectation over the spreads the IBM shows (CV up to ~0.6 in school-age
# children, less above), and the gamma closure is itself no better than that in
# the far lower tail. parameters$immunity_spread = FALSE evaluates every curve
# at the cell mean.
VIVAX_GH_Z <- c(-2.856970013872806, -1.355626179974266, 0,
                1.355626179974266, 2.856970013872806)
VIVAX_GH_W <- c(0.01125741132772071, 0.2220759220056126, 0.5333333333333329,
                0.2220759220056126, 0.01125741132772071)
immunity_spread_rule <- function(p) {
  on <- p$immunity_spread
  if (is.null(on)) on <- TRUE
  if (!(is.logical(on) && length(on) == 1 && !is.na(on)))
    stop("parameters$immunity_spread must be TRUE or FALSE.", call. = FALSE)
  if (on) list(spread_on = 1, n_q = length(VIVAX_GH_Z), qz = VIVAX_GH_Z, qw = VIVAX_GH_W)
  else list(spread_on = 0, n_q = 1L, qz = 0, qw = 1)
}

# Append n empty levels to the hypnozoite dimension (the third) of a seed array:
# the liver-stage-protected levels, which the equilibrium has no notion of --
# malariaEquilibriumVivax has no radical cure, and the IBM starts with nobody
# protected too.
pad_levels <- function(x, n) {
  if (n == 0) return(x)
  d <- dim(x)
  out <- array(0, replace(d, 3, d[3] + n))
  idx <- lapply(d, seq_len)
  do.call(`[<-`, c(list(out), idx, list(value = x)))
}

# Stages in each chain the vivax refractory windows are carried as (see the
# refractory stocks in the odin model). ICA's 4-day window then runs exactly,
# a stage a day; IAA's 44-day one as four 11-day stages. Holding IAA's to exactly
# 44 days instead moves infants' IAA by 2.7% of itself and no other output by
# more than 0.7%, at 2.4 times the cost of a vivax run.
VIVAX_WINDOW_STAGES <- 4L

# The vivax constants the model has no falciparum counterpart for, read from a
# malariasimulation vivax list. Everything a vivax list carries under a
# falciparum name -- da, dd, dt, cd, cu, ct, rc, uc, rm, pcm and the clinical
# curve phi0, phi1, ic0, kc -- reaches the model through the shared parameters.
# The refractory window uses the same integer rule as the others: ms tests
# (timestep - last_boosted) >= u, so the realised wait is ceil(u) days and
# u_eff = ceil(u) - 1.
vivax_scalars <- function(p) {
  # checked under malariasimulation's own names, before anything is computed
  # from them
  src <- c("gammal", "f", "b", "philm_min", "philm_max", "alm50", "klm", "dpcr_min",
           "dpcr_max", "apcr50", "kpcr", "ca", "ra", "ua")
  bad <- src[!vapply(src, function(nm) {
    x <- p[[nm]]
    length(x) == 1 && is.numeric(x) && is.finite(x)
  }, logical(1))]
  if (length(bad))
    stop("vivax parameter(s) missing or non-finite: ", paste(bad, collapse = ", "),
         ". Build the list with malariasimulation::get_parameters(parasite = ",
         "\"vivax\").", call. = FALSE)
  v <- list(gammal = p$gammal, ff = p$f, bv = p$b,
            philm_min = p$philm_min, philm_max = p$philm_max,
            alm50 = p$alm50, klm = p$klm,
            dpcr_min = p$dpcr_min, dpcr_max = p$dpcr_max,
            apcr50 = p$apcr50, kpcr = p$kpcr,
            cA_v = p$ca, d_iaa = p$ra, ua_eff = max(0, ceiling(p$ua) - 1))
  stopifnot(setequal(names(v), VIVAX_SCALARS))
  v
}

# The vivax scalar constants the odin model declares, in one place so the inert
# and live builds cannot fall out of step.
VIVAX_SCALARS <- c("gammal", "ff", "bv", "philm_min", "philm_max", "alm50", "klm",
                   "dpcr_min", "dpcr_max", "apcr50", "kpcr", "cA_v", "d_iaa", "ua_eff")

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
#' @param parameters a malariasimulation::get_parameters() list, for either
#'   parasite.
#' @param init_EIR target adult EIR (infectious bites per adult per year).
#' @param age_lower age-group lower edges in years.
#' @param n_ph,n_phc stage counts for the post-treatment and the
#'   chemoprevention prophylaxis chains; NULL picks the count that matches
#'   the variance of the drug's Weibull protection on the daily clock (see
#'   `erlang_stages()`), fixed at the seed's drug mix.
#' @param n_sub mosquito sub-steps per day.
#' @param timesteps simulation horizon in days. Used only to size the
#'   intervention time-series grids (vector control, PEV, TBV, carrying
#'   capacity), which are built to extend past `timesteps`.
#' @return list with `pars` (odin parameter list) and `meta` (grid metadata).
#' @noRd
build_inputs <- function(parameters, init_EIR, age_lower = default_age_lower(),
                         n_ph = NULL, n_phc = NULL, n_sub = 32L, timesteps = 3650) {
  p <- parameters
  # One parasite per run, as in malariasimulation. Much of what follows is shared
  # unchanged, because a vivax list carries the falciparum names wherever the two
  # coincide: da, dd, dt arrive as rA, rD, rT, and rm, pcm as dm, PM, so the
  # whole-day rate conversion, the drug series and the maternal-immunity factor
  # below are already the vivax ones under vivax. Falciparum-only constants the
  # vivax list lacks (b0, ib0, ...) keep malariaEquilibrium's defaults: finite
  # values for a falciparum block that is held empty.
  is_vivax <- identical(p$parasite, "vivax")
  if (!is.null(p$parasite) && !(identical(p$parasite, "falciparum") || is_vivax)) {
    stop("parameters$parasite must be 'falciparum' or 'vivax'.", call. = FALSE)
  }
  # malariasimulation cannot run these under vivax -- update_mass_drug_admin()
  # passes variables$id, which does not exist for vivax, and fails with "attempt
  # to apply non-function" (present in 3.0.0 and 3.0.1). A twin cannot offer a
  # result the original cannot produce, so refuse rather than run a pulse into
  # the empty falciparum chain and report it as a vivax answer.
  if (is_vivax) {
    chemo <- c(mda = isTRUE(p$mda), smc = isTRUE(p$smc), pmc = isTRUE(p$pmc))
    if (any(chemo))
      stop("parameters$parasite = 'vivax' cannot be combined with ",
           paste(toupper(names(chemo)[chemo]), collapse = ", "), ": malariasimulation ",
           "itself fails on this combination (update_mass_drug_admin() reads the ",
           "falciparum-only `id` variable), so there is no IBM result to reproduce.",
           call. = FALSE)
    # The vivax fields, checked before the seed and the guards below use them,
    # so a list that is not a vivax list fails here, saying so.
    invisible(vivax_scalars(p))
    for (nm in c("kmax", "ua")) {
      x <- p[[nm]]
      if (length(x) != 1 || !is.numeric(x) || !is.finite(x) ||
          (nm == "kmax" && (x < 1 || x != round(x))) || x < 0)
        stop("parameters$", nm, " must be a single ",
             if (nm == "kmax") "whole number >= 1" else "finite number >= 0",
             ". Build the list with malariasimulation::get_parameters(parasite = \"vivax\").",
             call. = FALSE)
    }
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
  # dd = dt = 5, not 5. Use that exit rate so fleet's dwell matches the IBM's. Must be
  # applied AFTER the eq_params merge above, or set_equilibrium's stored rates silently
  # overwrite it. Not applied to rP: ms models prophylaxis as a hazard multiplier, not
  # a compartment, so there is no per-day census of it.
  for (nm in c("rA", "rD", "rU", "rT")) {
    if (!is.null(eqp[[nm]]) && is.finite(eqp[[nm]]) && eqp[[nm]] > 0) {
      eqp[[nm]] <- 1 - exp(-eqp[[nm]])
    }
  }
  ## Treatment coverage for the seed, on malariasimulation's own two clocks. Its
  ## set_equilibrium() sizes the mosquitoes (init_foim, total_M) from the human
  ## equilibrium under the coverage in force at timestep 1, but the IBM draws its
  ## initial humans (calculate_eq(), variables.R) from the equilibrium under the
  ## coverage at timestep 0. With the usual set_clinical_treatment(timesteps = 1)
  ## the humans therefore start untreated, among mosquitoes sized for a treated
  ## population, and relax from there. fleet starts where the IBM starts.
  ft <- get_ft(p, 1)                       # mosquito sizing: coverage at timestep 1
  ft0 <- get_ft(p, 0)                      # initial humans: coverage at timestep 0
  dser <- drug_mix_series(p, eqp)          # time-varying drug_eff/cT/rP; seed = t=0 blend
  eqp[["cT"]] <- dser$seed$cT              # drug-linked treated infectivity for the seed
  rP <- dser$seed$rP                       # drug-linked prophylaxis rate for the seed
  chp <- chemoprevention_prophylaxis(p, eqp)   # chemoprevention prophylaxis: chain rate + stages
  rP_c <- chp$rate
  # Stage counts of the two prophylaxis chains: matched to the drug (mix)'s
  # Weibull protection unless overridden. 1 = a single stage.
  chk_stages <- function(k, nm) {
    if (length(k) != 1 || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k))
      stop("`", nm, "` must be a single whole number >= 1", if (nm %in% c("n_ph", "n_phc"))
        " (or NULL for the default)", ".", call. = FALSE)
    as.integer(k)
  }
  if (!is.null(n_ph)) n_ph <- chk_stages(n_ph, "n_ph")
  if (!is.null(n_phc)) n_phc <- chk_stages(n_phc, "n_phc")
  n_sub <- chk_stages(n_sub, "n_sub")
  # The three delay lines, in whole days plus a fraction read by interpolation,
  # as malariasimulation's LaggedValue reads a fractional lag.
  lag_ok <- function(x, nm, pos = FALSE) {
    if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x < 0 || (pos && x <= 0))
      stop("parameters$", nm, " must be a single finite number ", if (pos) "> 0" else ">= 0",
           " (days).", call. = FALSE)
  }
  lag_ok(p$de, "de"); lag_ok(p$delay_gam, "delay_gam"); lag_ok(p$dem, "dem", pos = TRUE)
  de_floor <- as.integer(floor(p$de))
  fl_floor <- as.integer(floor(p$delay_gam))
  # the incubation queue holds ceiling(dem) days, and the flux leaving it today
  # entered ceiling(dem) - 1 days ago
  n_tau <- as.integer(ceiling(p$dem) - 1L)
  if (!is.null(p$hold_init_EIR) && !(is.logical(p$hold_init_EIR) &&
                                     length(p$hold_init_EIR) == 1 && !is.na(p$hold_init_EIR)))
    stop("parameters$hold_init_EIR must be TRUE or FALSE.", call. = FALSE)
  # The two mean-field fidelity knobs, validated here rather than passed straight
  # through to dust2. Both went unchecked, and both silently accept nonsense:
  # `bite_dedup` is a 0/1 switch that the odin model uses as a blend weight, so 2
  # and -1 both run and both extrapolate outside either documented endpoint; and
  # `acquired_immunity_offset = 50` runs, changing clinical incidence 12-fold,
  # while -5 makes (IB + offset)^kb non-finite and poisons the run.
  bite_dedup <- if (is.null(p$bite_dedup)) 1 else p$bite_dedup
  if (is.logical(bite_dedup) && length(bite_dedup) == 1L && !is.na(bite_dedup)) {
    bite_dedup <- as.numeric(bite_dedup)   # FALSE is the natural spelling of "0"
  }
  if (length(bite_dedup) != 1 || !is.numeric(bite_dedup) || !is.finite(bite_dedup) ||
      !(bite_dedup == 0 || bite_dedup == 1))
    stop("parameters$bite_dedup must be 0 or 1 (TRUE/FALSE also accepted): 1 ",
         "reproduces the IBM's per-timestep bite deduplication, 0 lets every bite ",
         "infect independently. It is a switch, not a dial.", call. = FALSE)
  # The IBM adds +0.5 to each person's positive acquired immunity inside the Hill
  # curves. Falciparum evaluates its curves at the cell mean, where a per-person
  # +0.5 has no exact counterpart and 0 is the better match; vivax, with its
  # immunity spread, evaluates them at quadrature nodes that stand for people, so
  # it takes the IBM's own 0.5 (at a single node, the mean, it takes 0 too).
  acq_offset <- if (!is.null(p$acquired_immunity_offset)) p$acquired_immunity_offset else
    if (is_vivax && !isFALSE(p$immunity_spread)) 0.5 else 0
  if (length(acq_offset) != 1 || !is.numeric(acq_offset) || !is.finite(acq_offset) ||
      acq_offset < 0 || acq_offset > 1)
    stop("parameters$acquired_immunity_offset must be a single number in [0, 1]: 0 ",
         "(the falciparum default) is the better match where the Hill curves are read ",
         "at a cell mean, 0.5 (the vivax default) is the IBM's literal per-person +0.5.",
         call. = FALSE)
  # The seed takes raw coverage, as set_equilibrium() and create_variables()
  # hand it to malariaEquilibrium and malariaEquilibriumVivax, which send that
  # share of clinical cases to treatment. Efficacy acts in the run, so a treated
  # seed relaxes onto the model's own treated state, as the IBM's does, among
  # mosquitoes sized as the IBM sizes them.
  ft_seed <- ft
  ft0_seed <- ft0

  ## heterogeneity nodes
  n_het <- p$n_heterogeneity_groups
  if (is.null(n_het)) n_het <- 5
  if (isFALSE(p$enable_heterogeneity)) {
    zeta <- 1
    het_wt <- 1
    n_het <- 1
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
    # here would be a mean-matched fudge that moves fleet off the mechanism it
    # replicates.
  }

  ## age grid + demography (replicates malariaEquilibrium::human_equilibrium_no_het)
  n_age <- length(age_lower)
  age_days <- age_lower * 365
  width <- c(diff(age_days), Inf)
  r_ce <- 1 / width; r_ce[n_age] <- 0   # malariaEquilibrium's own convention
  eta <- eqp[["eta"]]
  age_mid <- c((age_days[-n_age] + age_days[-1]) / 2, age_days[n_age])
  psi <- 1 - eqp[["rho"]] * exp(-age_mid / eqp[["a0"]])
  # constant-hazard age structure = the malariaEquilibrium seed reference. This
  # one keeps r = 1/width: it is not fleet's age structure, it is the structure
  # the seed arrives on, and it has to match what malariaEquilibrium assumed.
  prop_ce <- numeric(n_age)
  for (i in seq_len(n_age)) {
    if (i == 1) prop_ce[i] <- eta / (r_ce[1] + eta)
    else prop_ce[i] <- prop_ce[i - 1] * r_ce[i - 1] / (r_ce[i] + eta)
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
  # AGEING RATE: exponentially fitted, not 1/width.
  #
  # A linear chain empties a band at r, so the chain's stationary ratio between
  # consecutive bands is r / (r + mu). With the obvious r = 1/h that is
  # 1 / (1 + mu*h), while the continuous truth is exp(-mu*h) -- and
  # 1/(1+x) > exp(-x) for every x > 0, so the obvious rate ALWAYS decays too
  # slowly and always leaves too many people alive at old ages. It is a
  # first-order (donor-cell) discretisation of the McKendrick advection.
  #
  # Choosing instead
  #     r = mu / (expm1(mu * h))
  # makes that ratio exp(-mu*h) exactly, so the stationary age structure is the
  # analytic band-integrated survival curve rather than an approximation to it.
  # This is exponential fitting, the standard cure for a first-order advection
  # scheme. It tends to 1/h as mu*h -> 0, so it is a no-op where bands are
  # narrow -- 0.3% in the narrowest bands -- and does its work in the wide ones,
  # 11% in the 5-year bands above 60. Against 20 IBM replicates it took the worst
  # band error from 7.4% to 1.3% and the population age structure from 9 of 11
  # bands inside the replicate band to 11 of 11.
  #
  # r depends on mu, which reads oddly -- ageing should not depend on dying.
  # It is not a biological rate: it is the coefficient that makes the DISCRETE
  # scheme reproduce the CONTINUOUS solution, and that solution involves mu.
  #
  # Fitted at the baseline (t = 0) mortality. Under set_demography() with
  # time-varying rates mu moves and this rate does not follow it, so the fit is
  # exact at the seed and an improvement, not an identity, later. r_age is a
  # constant in the odin model; making it time-varying would be the next step if
  # a transient ever needs it.
  mu0 <- mu_age
  r_age <- c(mu0[-n_age] / expm1(mu0[-n_age] * width[-n_age]), 0)
  # a zero hazard gives 0/0; the limit is 1/h
  flat <- !is.finite(r_age) | mu0 * width <= 0
  r_age[flat] <- r_ce[flat]
  r_age[n_age] <- 0
  # Each day an age group loses r_age + mu_age of its people. That is a fraction
  # of the group, so it cannot reach 1: a group narrower than about a day would
  # empty more than completely and every compartment in it would go negative.
  leave <- r_age + apply(mu_age_z, 1, max)
  # ...and it loses that out of the same stock as its own exit probability, so the
  # two together cannot exceed 1 either. The disease exits are bounded by the
  # largest infection probability, b0 (an infection moves someone out of S, A or
  # U), together with the progression and clearance probabilities; the prophylaxis
  # chains, below, exit at n x rate. Both are checked against the narrowest group.
  res <- resistance_series(p, eqp)
  room <- 1 - max(leave)
  too_narrow <- function(disease_exit, what = "infection and progression",
                         g = which.max(leave)) {
    stop("age group ", g, " ([", signif(age_days[g], 6), ", ", signif(age_days[g] + width[g], 6),
         ") days) is too narrow for the daily clock: it would lose ", signif(leave[g], 3),
         " of its people a day to ageing and death on top of up to ", signif(disease_exit, 3),
         " to ", what, ", more than it holds. Widen it in `age_lower`.",
         call. = FALSE)
  }
  # (The vivax counterpart depends on the EIR, and is checked once the seed EIR
  # is known, below.)
  if (!is_vivax) {
    disease_exit <- max(1 - (1 - eqp[["b0"]]) * (1 - max(eqp[["rA"]], eqp[["rU"]])),
                        eqp[["rD"]], eqp[["rT"]], max(res$rT_slow))
    if (disease_exit > room) too_narrow(disease_exit)
  }
  # Stage counts of the two prophylaxis chains. The defaults are sized at the
  # seed's drug mix and then held to what the shortest-protecting mix on the
  # schedule leaves room for; an explicit count that asks for more is refused
  # rather than run with negative occupancies. A chain with no drug behind it
  # carries nobody, so its count is not checked.
  cap <- function(rate) max(1L, as.integer(floor(room / max(rate))))
  if (is.null(n_ph)) n_ph <- min(dser$seed$n_ph, cap(dser$rP))
  if (is.null(n_phc)) n_phc <- min(chp$n_stages, cap(rP_c))
  n_phct <- min(chp$n_stages_t, cap(chp$rate_t))
  # Integer refractory windows: the IBM boosts when (timestep - last_boosted) >= u
  # with an integer timestep, so the realised wait is ceil(u) days and
  # u_eff = ceil(u) - 1; u = 0 is no window at all rather than a negative one.
  u_eff <- vapply(c(ub = "ub", uc = "uc", ud = "ud", uv = "uv"),
                  function(nm) max(0, ceiling(eqp[[nm]]) - 1), numeric(1))
  in_use <- c(n_ph = sum(unlist(p$clinical_treatment_coverages)) > 0,
              n_phc = isTRUE(p$smc) || isTRUE(p$mda) || isTRUE(p$pmc))
  for (ch in list(list("n_ph", n_ph, dser$rP), list("n_phc", n_phc, rP_c))) {
    if (in_use[[ch[[1]]]] && ch[[2]] * max(ch[[3]]) > room)
      stop("`", ch[[1]], " = ", ch[[2]], "` would have each stage of the prophylaxis chain ",
           "pass on ", signif(ch[[2]] * max(ch[[3]]), 3), " of its occupants a day, on top ",
           "of the ", signif(max(leave), 3), " the narrowest age group loses to ageing and ",
           "death, which is more than it holds: at most ", cap(ch[[3]]), " stages fit this ",
           "drug's protection on this age grid.", call. = FALSE)
  }
  # Under vivax the post-treatment chain is the vivax block's, sized and checked
  # as falciparum's just was, and carrying the batch dimension; the falciparum
  # chains collapse to one inert stage. Liver-stage protection after radical
  # cure adds n_ls levels to the hypnozoite dimension, a chain matched to the
  # drug's protection like the others. The model applies it after the day's
  # other transitions, so it needs only that a stage pass on at most everyone:
  # n_ls / mean_ls <= 1 for the shortest-protecting radical-cure mix on the
  # schedule.
  n_phv <- 1L
  n_ls <- 0L
  if (is_vivax) {
    n_phv <- n_ph
    n_ph <- 1L; n_phc <- 1L; n_phct <- 1L
    if (dser$n_ls > 0)
      n_ls <- min(dser$n_ls, max(1L, as.integer(floor(min(dser$mean_ls)))))
  }

  # equilibrium age structure under this mortality (McKendrick aging chain)
  prop <- numeric(n_age)
  prop[1] <- 1 / (r_age[1] + mu_age[1])
  if (n_age >= 2) for (i in 2:n_age) prop[i] <- prop[i - 1] * r_age[i - 1] / (r_age[i] + mu_age[i])
  prop <- prop / sum(prop)
  mean_psi <- sum(prop * psi)
  resc <- prop / prop_ce           # rescale the constant-eta seed to this age structure
  # Maternal immunity is inherited from a mother aged [20, 21) years:
  # sample_maternal_immunity() draws each newborn's mother uniformly among those
  # with trunc(age / 365) == 20 in its heterogeneity group (mortality_processes.R),
  # so a newborn inherits, in expectation, their population-weighted mean
  # immunity. mask20 is each age group's share of that year -- the fraction of
  # its width inside it, as band aggregation weighs a group -- and the model
  # weights by the live population. The default grid pins edges at 20 and 21, so
  # there it is 1 for the groups making up [20, 21) and 0 elsewhere; the open top
  # group counts whole if the year reaches into it, as a rendering band's does.
  lo20 <- age_days; hi20 <- c(age_days[-1], Inf)
  ov20 <- pmax(0, pmin(hi20, 21 * 365) - pmax(lo20, 20 * 365))
  mask20 <- ifelse(is.finite(hi20), ov20 / (hi20 - lo20), as.numeric(ov20 > 0))

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
  mum_v <- a_spp <- beta_eff <- g_s <- numeric(n_spp)
  for (s in seq_len(n_spp)) {
    fmr <- p$blood_meal_rates[[s]]
    mum_v[s] <- p$mum[[s]]
    beta_eff[s] <- eggs_laid(p$beta, mum_v[s], fmr)
    a_spp[s] <- p$Q0[[s]] * fmr
    # incubation survival exactly as malariasimulation (adult_mosquito_eqs.cpp:
    # incubation_survival = exp(-mu * tau)), applied once as the flux leaves the
    # incubation delay line
    g_s[s] <- exp(-mum_v[s] * dem)
  }

  ## human equilibrium seed (per het node) on the model age grid at adult EIR `EIR`,
  ## plus the adult-mosquito density at which sum_s a_s*Im_s = EIR/365 exactly under
  ## THIS age structure (survival g_s applies as the flux leaves incubation).
  seed_human <- function(EIR, ft_s = ft_seed) {
    S0 <- D0 <- A0 <- U0 <- Tr0 <- matrix(0, n_age, n_het)
    Ph0 <- array(0, c(n_age, n_het, n_ph)); Phc0 <- array(0, c(n_age, n_het, n_phc))
    IB_init <- ICA_init <- ID_init <- IVA_init <- matrix(0, n_age, n_het)
    cA0 <- matrix(0, n_age, n_het)
    for (j in seq_len(n_het)) {
      eqn <- malariaEquilibrium::human_equilibrium_no_het(
        EIR = EIR * zeta[j], ft = ft_s, p = eqp, age = age_lower
      )
      w <- het_wt[j]
      db <- solve_disease_block(
        eqn[, "FOI"], eqn[, "phi"], eqn[, "prop"], eqn[, "r"],
        eqp[["eta"]], eqp[["rA"]], eqp[["rD"]], eqp[["rU"]], eqp[["rT"]], rP, ft_s, n_ph
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
    ## equilibrium infectivity reaching mosquitoes (Xf0), normalised as the odin
    ## model normalises it, by the population mean of zeta x psi
    inf0 <- eqp[["cD"]] * D0 + cA0 * A0 + eqp[["cU"]] * U0 + eqp[["cT"]] * Tr0
    Xf0 <- sum(outer(psi, zeta) * inf0) / (mean_psi * sum(het_wt * zeta))
    ## per-species FOIM consistent with THIS model's seeded human infectivity
    foim0 <- a_spp * Xf0
    denom <- sum(a_spp * species_prop * foim0 * g_s / (foim0 + mum_v))
    list(S0 = S0, D0 = D0, A0 = A0, U0 = U0, Tr0 = Tr0, Ph0 = Ph0, Phc0 = Phc0,
         IB_init = IB_init, ICA_init = ICA_init, ID_init = ID_init, IVA_init = IVA_init,
         Xf0 = Xf0, foim0 = foim0, total_M = (EIR / 365) * hp / denom)
  }

  ## The vivax seed, the IBM's own: malariaEquilibriumVivax, exactly as
  ## set_equilibrium() and create_variables() call it. It returns every state
  ## as [age, het, batch], with the heterogeneity strata already built on the
  ## same Gauss-Hermite nodes fleet and the IBM use. Three steps take it to
  ## fleet's grid:
  ##   * it treats its ages as bin EDGES, so fleet's lower edges would give one
  ##     group fewer and no open top; one extra upper edge gives fleet's
  ##     absorbing top group a bin to take its within-age fractions from.
  ##   * the IBM draws each person's state and batch count from the
  ##     equilibrium's S, D, A, U and T (initial_state_vivax()), leaving out its
  ##     prophylaxis state P: nobody starts drug-protected, and the draw
  ##     renormalises within each age and heterogeneity group. So does this.
  ##   * its own age structure is not fleet's. The distribution WITHIN each age
  ##     and heterogeneity group over (batch, state) is correct conditional on
  ##     them, so it is re-weighted onto fleet's stationary structure `prop` and
  ##     the quadrature weights, as resc does for falciparum.
  ## With heterogeneity off, the IBM still solves the heterogeneous equilibrium:
  ## set_equilibrium() sizes the mosquitoes from the whole of it, and
  ## create_variables() draws everyone's state and batch count from its first
  ## group (groups <- rep(1, size)), the least bitten. Both are replicated: the
  ## humans come from group 1, renormalised within each age, and the mosquitoes
  ## from the heterogeneous seed.
  seed_human_vivax <- function(EIR, ft_s) {
    e <- malariaEquilibriumVivax::vivax_equilibrium(
      EIR = EIR, ft = ft_s, p = translate_vivax_parameters(p), age = c(age_lower, VIVAX_TOP_AGE))
    st <- e$states
    n_het_eq <- dim(st$S)[2]
    if (n_het_eq == n_het) {
      zeta_eq <- zeta; wt_eq <- het_wt
    } else {
      gq_eq <- malariaEquilibrium::gq_normal(n_het_eq)
      zeta_eq <- exp(gq_eq$nodes * sqrt(eqp[["s2"]]) - eqp[["s2"]] / 2)
      wt_eq <- gq_eq$weights
    }
    tot <- apply(st$S + st$D + st$A + st$U + st$T, c(1, 2), sum)
    rw_eq <- function(x) sweep(x, c(1, 2), outer(prop, wt_eq) / tot, `*`)
    ## infectivity reaching mosquitoes and the mosquito sizing exactly as
    ## seed_human(), summed over the batch dimension; A carries the constant ca
    xf <- function(D, A, U, Tr, z, w) {
      inf0 <- apply(eqp[["cD"]] * D + p$ca * A + eqp[["cU"]] * U + eqp[["cT"]] * Tr, c(1, 2), sum)
      sum(outer(psi, z) * inf0) / (mean_psi * sum(w * z))
    }
    foim0 <- a_spp * xf(rw_eq(st$D), rw_eq(st$A), rw_eq(st$U), rw_eq(st$T), zeta_eq, wt_eq)
    denom <- sum(a_spp * species_prop * foim0 * g_s / (foim0 + mum_v))
    if (n_het_eq != n_het) {
      grp1 <- function(x) x[, 1, , drop = FALSE]
      st <- lapply(st, function(x) if (length(dim(x)) == 3L) grp1(x) else x)
      tot <- tot[, 1, drop = FALSE]
    }
    rw <- function(x) sweep(x, c(1, 2), outer(prop, het_wt) / tot, `*`)
    Sv0 <- rw(st$S); Uv0 <- rw(st$U); Av0 <- rw(st$A); Dv0 <- rw(st$D); Trv0 <- rw(st$T)
    ## the humans' own infectivity, where the human-infectivity line starts
    Xf0 <- xf(Dv0, Av0, Uv0, Trv0, zeta, het_wt)
    ## immunity as stocks (the model holds J = N * I): the equilibrium's mean per
    ## cell times its population. Everyone in a cell starts at its mean, as
    ## initial_immunity() starts the IBM, so the within-cell spread builds up from
    ## zero in both.
    Ncell <- Sv0 + Uv0 + Av0 + Dv0 + Trv0
    list(Sv0 = Sv0, Dv0 = Dv0, Av0 = Av0, Uv0 = Uv0, Trv0 = Trv0,
         JAv0 = st$IAA * Ncell, JCv0 = st$ICA * Ncell,
         KAv0 = st$IAA^2 * Ncell, KCv0 = st$ICA^2 * Ncell,
         Xf0 = Xf0, foim0 = foim0, total_M = (EIR / 365) * hp / denom)
  }
  seed <- if (is_vivax) seed_human_vivax else seed_human

  ## Which EIR to seed at. Under the default demography fleet's own sizing is
  ## malariasimulation's formula evaluated at fleet's own equilibrium (its age
  ## grid, whole-day rates, drug-linked cT and ft*eff), which lands within ~1-2% of
  ## the IBM's total_M, so the seed is at init_EIR. Under a custom demography the two conventions part
  ## ways: set_equilibrium() sizes total_M from the equilibrium under the DEFAULT
  ## exponential age structure (compatibility.R set_equilibrium ->
  ## equilibrium_total_M; custom mortality never enters), and the IBM then drifts to
  ## whatever transmission that mosquito density supports under the custom
  ## mortality. Replicate that: take the IBM's total_M and find the EIR at which
  ## fleet's equilibrium under the custom age structure has exactly that density --
  ## so the same parameter list realises the same transmission in both models.
  ## parameters$hold_init_EIR = TRUE seeds at init_EIR instead (init_EIR is then
  ## the EIR fleet realises).
  eir_seed <- init_EIR
  total_M_ibm <- NA_real_
  sub_threshold <- FALSE
  if (isTRUE(p$custom_demography) && !isTRUE(p$hold_init_EIR)) {
    total_M_ibm <- if (is_vivax) ibm_total_M_vivax(p, init_EIR, ft) else
      ibm_total_M(p, init_EIR, eqp_ibm, ft)
    # fleet's total_M(EIR) is increasing in EIR and tends to the transmission
    # threshold M_crit as EIR -> 0, so a root exists iff M_crit < total_M_ibm:
    # test at a near-zero EIR rather than at an arbitrary fraction of init_EIR.
    f <- function(lE) seed(exp(lE), ft_seed)$total_M - total_M_ibm
    lo <- log(init_EIR * 1e-6); f_lo <- f(lo)
    if (f_lo > 0) {
      # The IBM's mosquito density is below fleet's transmission threshold under
      # this demography: transmission is not sustainable, and the IBM decays to
      # elimination from its seed. Keep the IBM's density (so fleet decays too) and
      # start the humans at a near-zero EIR.
      warning("Under this custom demography the mosquito population that ",
              "set_equilibrium() implies is below fleet's transmission threshold (",
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
  ## The humans come from the timestep-0 equilibrium and the mosquitoes from the
  ## timestep-1 one (see ft0 above); the two coincide unless treatment starts at
  ## timestep 1 or later. The human-infectivity delay line starts at the humans'
  ## own infectivity, as the IBM's LaggedValue defaults to it (processes.R), while
  ## the mosquito states, the EIR line and the incubation queue start from init_foim.
  hs <- seed(eir_seed, ft0_seed)
  hm <- if (ft0_seed == ft_seed) hs else seed(eir_seed, ft_seed)
  if (is_vivax) {
    # the falciparum block is held empty under vivax: zero state, one inert
    # prophylaxis stage each (the chains were collapsed to 1 above)
    zf <- matrix(0, n_age, n_het)
    S0 <- D0 <- A0 <- U0 <- Tr0 <- zf
    IB_init <- ICA_init <- ID_init <- IVA_init <- zf
    Ph0 <- array(0, c(n_age, n_het, n_ph)); Phc0 <- array(0, c(n_age, n_het, n_phc))
  } else {
    S0 <- hs$S0; D0 <- hs$D0; A0 <- hs$A0; U0 <- hs$U0; Tr0 <- hs$Tr0
    Ph0 <- hs$Ph0; Phc0 <- hs$Phc0
    IB_init <- hs$IB_init; ICA_init <- hs$ICA_init; ID_init <- hs$ID_init; IVA_init <- hs$IVA_init
  }
  Xf0 <- hs$Xf0; foim0 <- hm$foim0
  total_M <- if (sub_threshold) total_M_ibm else hm$total_M

  ## The vivax counterpart of the too-narrow guard above. Each day a vivax cell
  ## loses its ageing and death fraction and the day's competing event --
  ## infection or progression, or a prophylaxis stage's exit -- out of the same
  ## stock; the batch clearance and the liver-stage clock act on what is left.
  ## Unlike falciparum's, whose deduplicated bites cap infection at b0, a vivax
  ## person's infection probability rises towards 1 with the EIR, so the bound
  ## is checked at the seed's EIR, for each age group at its own exposure in its
  ## most exposed stratum: the narrowest groups are the youngest, whom mosquitoes
  ## bite least. A seasonal peak can go above it.
  if (is_vivax) {
    # the seed's EPS (mpsi / mzp = 1 / sum(w zeta) on the stationary population)
    eps_g <- eir_seed / 365 * max(zeta) * psi / sum(het_wt * zeta)
    i_g <- 1 - exp(-(p$b * eps_g + p$kmax * p$f))
    ev_g <- pmax(1 - (1 - i_g) * (1 - max(eqp[["rA"]], eqp[["rD"]], 1 - exp(-1 / p$dpcr_min))),
                 eqp[["rT"]], max(res$rT_slow), n_phv * max(dser$rP))
    over <- leave + ev_g - 1
    if (max(over) > 0) {
      g <- which.max(over)
      too_narrow(ev_g[g], paste0("infection and progression at EIR ", signif(eir_seed, 3)),
                 g = g)
    }
  }

  ## There is deliberately no het-integrated malariaEquilibrium solve here: the
  ## FOIM the mosquitoes are seeded at is `foim0`, computed above from this model's
  ## own timestep-1 human infectivity, fleet's counterpart of the IBM's init_foim.
  ## A reference solve for a diagnostic belongs in the diagnostic, not on the hot
  ## path.

  ME0 <- ML0 <- MP0 <- Sm0 <- Im0 <- Em0 <- K0 <- numeric(n_spp)
  for (s in seq_len(n_spp)) {
    m_s <- species_prop[s] * total_M
    counts <- initial_mosquito_counts(p, s, foim0[s], m_s) / hp
    ME0[s] <- counts[["E"]]; ML0[s] <- counts[["L"]]; MP0[s] <- counts[["P"]]
    Sm0[s] <- counts[["Sm"]]                       # = m_s*mum/(foim+mum)/hp
    Em0[s] <- counts[["Pm"]]                       # incubating stock (ms E state)
    Im0[s] <- counts[["Im"]]                       # = m_s*foim/(foim+mum)*exp(-mum*dem)/hp
    K0[s] <- calculate_carrying_capacity(p, m_s, s) / hp
  }
  # The delay lines start full of the seed's own values, as the IBM's LaggedValue
  # defaults and its incubation queue do: the EIR a * Im, the incubation inflow
  # S * foim, and the human infectivity. When the two clocks agree, the mosquito
  # model's seed is its fixed point at constant foim, day one included.
  eir0 <- a_spp * Im0
  inc0 <- Sm0 * foim0

  # A species with proportion 0 (common in site files: one dominant vector, others
  # at 0) gives K0 = 0, so its larval carrying capacity Kcap = 0. Its larval count
  # is also 0, so the larval update evaluates 0 / Kcap = 0/0 = NaN and the run is
  # lost. Floor K0 to a negligible positive value: the species stays inert (0
  # mosquitoes, 0 biting) but every rate is finite.
  if (max(K0) > 0) K0 <- pmax(K0, max(K0) * 1e-9)

  ## Intervention time series. Each starts at its baseline at t = 0, so the
  ## equilibrium seed above is preserved; interventions act from their scheduled
  ## timesteps. (The vector-control t = 0 row equals the baseline a/mu used above.)
  ##
  ## The daily update that runs from time t to t + 1 is the IBM's timestep t + 1,
  ## and reads every series at time t. Two kinds of schedule meet that clock
  ## differently:
  ##   * a DEPLOYMENT (nets, spraying, vaccination, a carrying-capacity change) is
  ##     an event in the IBM, which fires after the day's processes and takes
  ##     effect the next day. Those series are built with the deployment's own
  ##     timestep as the knot that includes it, so time t (IBM day t + 1) sees a
  ##     deployment made on day t. They are used as built.
  ##   * a SETTING read by a process on the day itself (treatment coverage, the
  ##     drug mix, resistance, death rates: match_timestep(ts, timestep)) takes
  ##     effect on its own timestep. Those knots move back one day, so that IBM
  ##     day t + 1 reads the value scheduled for timestep t + 1.
  trt <- treatment_series(p)
  vc <- vector_control_series(p, timesteps)
  pevs <- pev_series(p, age_mid, timesteps)
  tbvs <- tbv_series(p, age_mid, timesteps)
  ccs <- carrying_capacity_series(p, K0, timesteps)

  pars <- list(
    n_age = n_age, n_het = n_het, n_spp = n_spp,
    n_ph = n_ph, n_phc = n_phc, n_phct = n_phct, n_sub = n_sub,
    de_floor = de_floor, de_frac = p$de - de_floor, n_eirh = de_floor + 1L,
    fl_floor = fl_floor, fl_frac = p$delay_gam - fl_floor, n_infh = fl_floor + 1L,
    n_tau = n_tau, n_inch = max(1L, n_tau),
    r_age = r_age, psi = psi,
    n_mut = ncol(mu_age_z), mu_age_t = mu_age_t - 1, mu_age_z = mu_age_z,
    age_mid = age_mid, mask20 = mask20,
    icm_factor = icm_factor, ivm_factor = ivm_factor,
    zeta = zeta, het_wt = het_wt,
    rA = eqp[["rA"]], rD = eqp[["rD"]], rU = eqp[["rU"]],
    rT = eqp[["rT"]], rP_c = rP_c, rT_slow = res$rT_slow,
    rP_ct = chp$rate_t, rT_slow_c = .chemo_rT_slow(p, eqp[["rT"]]),
    # slow-clearance fraction at t = 0, used to split the Tr seed between the fast
    # and slow treated compartments (the IBM assigns each treated individual to one
    # or the other by a Bernoulli draw, so Tr is a two-component mixture).
    spc0 = if (length(res$spc)) res$spc[1] else 0,
    d_ib = eqp[["db"]], d_ica = eqp[["dc"]], d_id = eqp[["dd"]], d_iva = eqp[["dv"]],
    # integer refractory windows (u_eff, above), and the same as array indices for
    # the refractory closure, which under vivax runs over the empty falciparum
    # block and is kept to a single day
    ub_eff = u_eff[["ub"]], uc_eff = u_eff[["uc"]], ud_eff = u_eff[["ud"]], uv_eff = u_eff[["uv"]],
    n_w = if (is_vivax) 1L else as.integer(max(1, u_eff[c("uc", "ud", "uv")])),
    wc = if (is_vivax) 1L else as.integer(max(1, u_eff[["uc"]])),
    wd = if (is_vivax) 1L else as.integer(max(1, u_eff[["ud"]])),
    wv = if (is_vivax) 1L else as.integer(max(1, u_eff[["uv"]])),
    # Acquired-immunity offset in the Hill calls (acq_offset, above): 0 for
    # falciparum's cell means, 0.5 for vivax's per-person quadrature nodes.
    acq_offset = acq_offset,
    # Reproduce the IBM's per-timestep bite deduplication by default;
    # parameters$bite_dedup = 0 lets every bite infect independently instead.
    bite_dedup = bite_dedup,
    b0 = eqp[["b0"]], b1 = eqp[["b1"]], ib0 = eqp[["IB0"]], kb = eqp[["kb"]],
    phi0 = eqp[["phi0"]], phi1 = eqp[["phi1"]], ic0 = eqp[["IC0"]], kc = eqp[["kc"]],
    d1 = eqp[["d1"]], id0 = eqp[["ID0"]], kd = eqp[["kd"]],
    fd0 = eqp[["fd0"]], ad0 = eqp[["ad0"]], gd = eqp[["gd"]],
    theta0 = eqp[["theta0"]], theta1 = eqp[["theta1"]], iv0 = eqp[["IV0"]], kv = eqp[["kv"]],
    fv0 = eqp[["fv0"]], av = eqp[["av"]], gammav = eqp[["gammav"]],
    cD = eqp[["cD"]], cU = eqp[["cU"]], g_inf = eqp[["g_inf"]],
    PM = eqp[["PM"]], PVM = eqp[["PVM"]],
    n_ftt = length(trt$times), ft_times = trt$times - 1, ft_vals = trt$vals,
    n_dmix = length(dser$times), dmix_times = dser$times - 1,
    cT_vals = dser$cT, drug_eff_vals = dser$drug_eff, rP_vals = dser$rP,
    hyp_vals = dser$hyp, eff_hyp_vals = dser$eff_hyp, mean_ls_vals = dser$mean_ls,
    n_rest = length(res$times), res_times = res$times - 1,
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
    dem = p$dem,
    S0 = S0, D0 = D0, A0 = A0, U0 = U0, Tr0 = Tr0, Ph0 = Ph0, Phc0 = Phc0,
    IB_init = IB_init, ICA_init = ICA_init, ID_init = ID_init, IVA_init = IVA_init,
    ME0 = ME0, ML0 = ML0, MP0 = MP0, Sm0 = Sm0, Em0 = Em0, Im0 = Im0,
    eir0 = eir0, inc0 = inc0, inf0 = Xf0
  )
  pars <- c(pars, if (is_vivax) {
    stopifnot(dim(hs$Sv0)[3] == p$kmax + 1)
    # rc_on sizes the radical-cure transfers: full size only if some drug on
    # the schedule gives radical cure, one inert cell otherwise
    n_bat <- as.integer(p$kmax) + 1L
    pl <- function(x) pad_levels(x, n_ls)
    ua_eff <- max(0, ceiling(p$ua) - 1); uc_eff <- u_eff[["uc"]]
    c(list(pf_on = 0, pv_on = 1, n_age_v = n_age, n_het_v = n_het,
           n_bat = n_bat, n_hyp = n_bat + n_ls, n_phv = n_phv,
           rc_on = as.numeric(any(dser$hyp > 0)),
           n_ra = as.integer(max(1, min(ua_eff, VIVAX_WINDOW_STAGES))),
           n_rc = as.integer(max(1, min(uc_eff, VIVAX_WINDOW_STAGES))),
           Sv0 = pl(hs$Sv0), Dv0 = pl(hs$Dv0), Av0 = pl(hs$Av0), Uv0 = pl(hs$Uv0),
           Trv0 = pl(hs$Trv0), Phv0 = array(0, c(n_age, n_het, n_bat + n_ls, n_phv)),
           JAv0 = pl(hs$JAv0), JCv0 = pl(hs$JCv0), KAv0 = pl(hs$KAv0), KCv0 = pl(hs$KCv0)),
      immunity_spread_rule(p),
      vivax_scalars(p))
  } else inert_vivax_block())

  list(
    pars = pars,
    # `meta` carries only what render_output() and apply_chemoprevention_pulse()
    # read, plus the seeding scalars the tests assert on. Nothing here is kept
    # "for reference": a field nobody reads is a solve nobody needed.
    meta = list(
      n_age = n_age, n_het = n_het,
      age_mid = age_mid, age_lo = age_days, age_hi = age_days + width,
      prop = prop, init_EIR = init_EIR, eir_seed = eir_seed,
      total_M = total_M, total_M_ibm = total_M_ibm, ft = ft,
      human_population = hp, parameters = p, eqp = eqp, parasite = if (is_vivax) "vivax" else "falciparum",
      n_ph = n_ph, n_phc = n_phc, n_phct = n_phct, n_phv = n_phv, n_ls = n_ls,
      # what a chemoprevention pulse needs to find the LM-detectable asymptomatic
      # and their infectivity (the model's q and cA)
      detect = list(d1 = eqp[["d1"]], id0 = eqp[["ID0"]], kd = eqp[["kd"]],
                    fd = 1 - (1 - eqp[["fd0"]]) / (1 + (age_mid / eqp[["ad0"]])^eqp[["gd"]]),
                    cU = eqp[["cU"]], cD = eqp[["cD"]], g_inf = eqp[["g_inf"]])
    )
  )
}
