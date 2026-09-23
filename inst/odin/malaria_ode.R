# =============================================================================
# fleet: mean-field (ODE) twin of malariasimulation (P. falciparum)
# Core transmission model (Phase 1): human S/D/A/U/Tr/Ph over [age, het] with
# the six Griffin immunity functions, coupled to the compartmental mosquito
# model (E/L/P/Sm/[EIP chain]/Im per species).
#
# Lags are implemented as Erlang (linear-chain) filters rather than delay():
# odin2's delay() warm-up seeding is unreliable with multiple delays, whereas
# chains are pure ODE, seed exactly at equilibrium, and are robust. The two
# human transmission lags are loss-free (exact at equilibrium); the mosquito
# EIP chain is a loss-free delay with the incubation survival applied at exit.
# The two prophylaxis compartments (Ph post-treatment, Ph_c chemoprevention) are
# Erlang chains too, so their exit-time distributions approximate the IBM's
# Weibull protection rather than an exponential at its mean.
# =============================================================================

## ---- dimensions -----------------------------------------------------------
n_age <- parameter(constant = TRUE)
n_het <- parameter(constant = TRUE)
n_spp <- parameter(constant = TRUE)
n_eir <- parameter(constant = TRUE)   # EIR-lag chain stages
n_foim <- parameter(constant = TRUE)  # FOIM-lag chain stages
n_eip <- parameter(constant = TRUE)   # EIP chain stages
n_ph <- parameter(constant = TRUE)    # post-treatment prophylaxis chain stages
n_phc <- parameter(constant = TRUE)   # chemoprevention prophylaxis chain stages

## ---- parasite species -------------------------------------------------------
# One parasite per run, as in malariasimulation, where `parasite` is a scalar.
# Both blocks are always compiled; the one not in use is held empty and inert.
#
# The vivax block has its own dimensions so that under falciparum it collapses to
# a single cell -- n_age_v = n_het_v = n_hyp = 1 -- and costs a handful of states.
# The falciparum block is NOT collapsed under vivax. It is dimensioned by n_age
# and n_het, which it shares with the demography, biting and vaccine arrays the
# vivax block also needs at full size; collapsing it would mean giving every
# falciparum equation dimensions of its own. So under vivax it sits at full
# size, empty -- about 3,000 inert states, a small fraction of a vivax run.
#
# An empty block divides by its own zero population, and a NaN derivative
# anywhere fails the whole run even if nothing reaches a shared quantity -- so
# multiplying its output by zero is not enough (0 x NaN is NaN). Divisions by a
# cell population therefore carry pop_floor in the denominator: an empty cell
# gives 0/pop_floor = 0, and in an occupied one it is an exact no-op, since the
# smallest real cell (~1e-5) has an ulp near 1e-21.
#
# NOT `if (pf_on > 0.5) ... else 0`. That was tried first, and although it is
# correct it made every falciparum run 10x slower at an unchanged step count:
# the branch sat in the loop that also computes 1 - exp(-EPS) and 1 - exp(-FOI),
# and stopped the compiler optimising it. Keep species switches out of loops
# that carry transcendentals.
pop_floor <- 1e-300
pf_on <- parameter(constant = TRUE)   # 1 under falciparum, 0 under vivax
pv_on <- parameter(constant = TRUE)   # 1 under vivax, 0 under falciparum
n_age_v <- parameter(constant = TRUE) # n_age under vivax, 1 otherwise
n_het_v <- parameter(constant = TRUE) # n_het under vivax, 1 otherwise
n_bat <- parameter(constant = TRUE)   # batch levels: kmax + 1 under vivax, 1 otherwise
n_hyp <- parameter(constant = TRUE)   # all levels: n_bat plus liver-stage protection

## ---- grid / demography (data) ---------------------------------------------
r_age <- parameter(); psi <- parameter()
age_mid <- parameter(); mask20 <- parameter()
icm_factor <- parameter(); ivm_factor <- parameter()
dim(r_age, psi, age_mid, mask20, icm_factor, ivm_factor) <- n_age
# age-specific mortality, time-varying for custom demography (constant otherwise):
# supplied as [n_age, n_mut] with time as the LAST dimension, like pev_vals.
n_mut <- parameter(constant = TRUE)
mu_age_t <- parameter(); dim(mu_age_t) <- n_mut
mu_age_z <- parameter(); dim(mu_age_z) <- c(n_age, n_mut)
mu_age <- interpolate(mu_age_t, mu_age_z, "constant"); dim(mu_age) <- n_age

## ---- heterogeneity (data) -------------------------------------------------
zeta <- parameter(); het_wt <- parameter()
dim(zeta, het_wt) <- n_het

## ---- human rate constants (from eq_params) --------------------------------
rA <- parameter(); rD <- parameter(); rU <- parameter()
rT <- parameter()
d_ib <- parameter(); d_ica <- parameter(); d_id <- parameter(); d_iva <- parameter()
# effective refractory in whole days: ceil(u) - 1 (see the boosting block below)
ub_eff <- parameter(); uc_eff <- parameter(); ud_eff <- parameter(); uv_eff <- parameter()
# malariasimulation adds +0.5 to POSITIVE acquired immunity before the b/phi/theta
# Hill calls (human_infection.R blood_immunity/clinical_immunity/severe_immunity),
# but NOT before q and NOT on the maternal term. That is a PER-INDIVIDUAL detail: in
# the mean field it does NOT translate to +0.5 on the stratum mean. An A/B against the
# ms 3.0.0 IBM ensemble mean (controlled burn-in over an EIR sweep) shows offset 0 is a
# markedly better match for clinical/severe incidence than 0.5, so the default is 0.
# Set 0.5 to reproduce the IBM's literal per-individual Hill functions instead.
acq_offset <- parameter(0)
# 1 = reproduce the IBM's per-timestep bite deduplication (saturating infection
# hazard, see the FOI block); 0 = the plain linear b*EPS, for which the
# malariaEquilibrium seed drifts least (it is still not an exact fixed point).
bite_dedup <- parameter(1)
b0 <- parameter(); b1 <- parameter(); ib0 <- parameter(); kb <- parameter()
phi0 <- parameter(); phi1 <- parameter(); ic0 <- parameter(); kc <- parameter()
d1 <- parameter(); id0 <- parameter(); kd <- parameter()
fd0 <- parameter(); ad0 <- parameter(); gd <- parameter()
theta0 <- parameter(); theta1 <- parameter(); iv0 <- parameter(); kv <- parameter()
fv0 <- parameter(); av <- parameter(); gammav <- parameter()
cD <- parameter(); cU <- parameter(); g_inf <- parameter()
PM <- parameter(); PVM <- parameter()

## ---- clinical treatment coverage ft(t) (time-varying, step) ---------------
n_ftt <- parameter(constant = TRUE)
ft_times <- parameter(); ft_vals <- parameter()
dim(ft_times, ft_vals) <- n_ftt
ft <- interpolate(ft_times, ft_vals, "constant")
# drug-linked efficacy / treated infectivity (cT) / prophylaxis rate (rP), all
# time-varying so a first-line drug switch changes the mix, not just total ft.
n_dmix <- parameter(constant = TRUE)
dmix_times <- parameter(); dim(dmix_times) <- n_dmix
cT_vals <- parameter(); drug_eff_vals <- parameter(); rP_vals <- parameter()
dim(cT_vals, drug_eff_vals, rP_vals) <- n_dmix
cT <- interpolate(dmix_times, cT_vals, "constant")
drug_eff <- interpolate(dmix_times, drug_eff_vals, "constant")
rP <- interpolate(dmix_times, rP_vals, "constant")
# P. vivax radical cure on the same schedule: the radically cured share of
# treatment, the share also cleared in the blood, and the mean liver-stage
# protection. All zero or unused under falciparum.
hyp_vals <- parameter(); eff_hyp_vals <- parameter(); mean_ls_vals <- parameter()
dim(hyp_vals, eff_hyp_vals, mean_ls_vals) <- n_dmix
hyp_mix <- interpolate(dmix_times, hyp_vals, "constant")
eff_hyp <- interpolate(dmix_times, eff_hyp_vals, "constant")
mean_ls <- interpolate(dmix_times, mean_ls_vals, "constant")
# Prophylaxis is an Erlang chain of n_ph stages, each left at rate n_ph*rP, so the
# chain's mean is 1/rP and its exit-time distribution approximates the IBM's
# Weibull (build_inputs picks n_ph from the Weibull shape). The IBM applies the
# Weibull survival W(t - t_drug) to each treated person's infection probability;
# the cohort's mean protection at lag t is therefore exactly W(t), which a chain
# whose survival matches W reproduces in the mean field. A single exponential
# stage at the Weibull mean leaks protection early (SP-AQ: 42% still protected
# at day 30 against the Weibull's 70%) and under-estimates SMC impact.
rPk <- rP * n_ph

## ---- antimalarial resistance (time-varying, step) -------------------------
# etf_t = fraction of would-be-treated that fail early (-> stay clinical D);
# spc_t = fraction of treated with slow parasite clearance (longer Tr).
res_times <- parameter(); etf_vals <- parameter(); spc_vals <- parameter()
n_rest <- parameter(constant = TRUE)
dim(res_times, etf_vals, spc_vals) <- n_rest
etf <- interpolate(res_times, etf_vals, "constant")
spc <- interpolate(res_times, spc_vals, "constant")
rT_slow <- parameter()          # 1/dt_slow_parasite_clearance
spc0 <- parameter(0)            # slow-clearance fraction at t = 0 (splits the Tr seed)
ft_eff <- ft * drug_eff * (1 - etf)   # coverage * efficacy * (1 - early-treatment-failure)

## ---- mosquito parameters (data + time-varying vector control) -------------
del <- parameter(); dl <- parameter(); dpl <- parameter()
me <- parameter(); ml <- parameter(); mup <- parameter(); mosq_gamma <- parameter()
# larval carrying capacity, time-varying under seasonality / set_carrying_capacity
n_cct <- parameter(constant = TRUE)
cc_times <- parameter(); dim(cc_times) <- n_cct
K_vals <- parameter(); dim(K_vals) <- c(n_spp, n_cct)
Kcap <- interpolate(cc_times, K_vals, "linear"); dim(Kcap) <- n_spp
# a (human blood-meal rate) and mu (adult death) vary in time under nets/IRS;
# supplied per-species. odin2 array interpolate() requires time as the LAST
# dimension, so a_vals/mum_vals are [n_spp, n_time].
n_vct <- parameter(constant = TRUE)
vc_times <- parameter(); dim(vc_times) <- n_vct
a_vals <- parameter(); mum_vals <- parameter()
dim(a_vals) <- c(n_spp, n_vct)
dim(mum_vals) <- c(n_spp, n_vct)
# linear (not constant): the IBM recomputes a/mu every timestep, so IRS logistic
# phase decay within a spray round is a ramp, not a step. Flat at baseline (equal
# knots), so the equilibrium is preserved.
a_spp <- interpolate(vc_times, a_vals, "linear")
mum <- interpolate(vc_times, mum_vals, "linear")
dim(a_spp) <- n_spp
dim(mum) <- n_spp
# oviposition rate: eggs_laid(beta, mu, f) is algebraically identical to beta for
# all mu, f, so beta_eff is constant (vector control acts via mu, not fecundity).
beta_eff <- parameter(); dim(beta_eff) <- n_spp

## ---- vaccines: per-age FOI (PEV) and infectivity (TBV) multipliers ---------
n_pevt <- parameter(constant = TRUE)
pev_times <- parameter(); dim(pev_times) <- n_pevt
pev_vals <- parameter(); dim(pev_vals) <- c(n_age, n_pevt)
pev_factor <- interpolate(pev_times, pev_vals, "linear"); dim(pev_factor) <- n_age
n_tbvt <- parameter(constant = TRUE)
tbv_times <- parameter(); dim(tbv_times) <- n_tbvt
tbv_fU_vals <- parameter(); tbv_fA_vals <- parameter()
tbv_fD_vals <- parameter(); tbv_fT_vals <- parameter()
dim(tbv_fU_vals) <- c(n_age, n_tbvt); dim(tbv_fA_vals) <- c(n_age, n_tbvt)
dim(tbv_fD_vals) <- c(n_age, n_tbvt); dim(tbv_fT_vals) <- c(n_age, n_tbvt)
tbv_fU <- interpolate(tbv_times, tbv_fU_vals, "linear"); dim(tbv_fU) <- n_age
tbv_fA <- interpolate(tbv_times, tbv_fA_vals, "linear"); dim(tbv_fA) <- n_age
tbv_fD <- interpolate(tbv_times, tbv_fD_vals, "linear"); dim(tbv_fD) <- n_age
tbv_fT <- interpolate(tbv_times, tbv_fT_vals, "linear"); dim(tbv_fT) <- n_age

## ---- lag time constants (data) --------------------------------------------
de <- parameter(); tl <- parameter(); dem <- parameter()
reir <- n_eir / de     # per-stage rate of the EIR-lag chain
# per-stage rate of the FOIM-lag chain. P. vivax has no gametocyte lag
# (delay_gam = 0), which malariasimulation reads as the current infectivity;
# the chain is then held static and bypassed (inf_lag below).
rfoim <- if (tl > 0) n_foim / tl else 0
reip <- n_eip / dem    # per-stage rate of the EIP chain

## ---- derived per-age quantities -------------------------------------------
re[] <- r_age[i] + mu_age[i]
fd[] <- 1 - (1 - fd0) / (1 + (age_mid[i] / ad0)^gd)
fv[] <- 1 - (1 - fv0) / (1 + (age_mid[i] / av)^gammav)
dim(re, fd, fv) <- n_age

## ---- maternal immunity (algebraic; matches malariaEquilibrium) ------------
ICAmask[, ] <- ICA[i, j] * mask20[i]
IVAmask[, ] <- IVA[i, j] * mask20[i]
dim(ICAmask, IVAmask) <- c(n_age, n_het)
ICA20[] <- sum(ICAmask[, i])
IVA20[] <- sum(IVAmask[, i])
dim(ICA20, IVA20) <- n_het
ICM[, ] <- PM * ICA20[j] * icm_factor[i]
IVM[, ] <- PVM * IVA20[j] * ivm_factor[i]
dim(ICM, IVM) <- c(n_age, n_het)

## ---- immunity -> probability (Hill) functions -----------------------------
b[, ] <- b0 * (b1 + (1 - b1) / (1 + ((IB[i, j] + acq_offset) / ib0)^kb))
phi[, ] <- phi0 * (phi1 + (1 - phi1) /
  (1 + (((ICA[i, j] + acq_offset) + ICM[i, j]) / ic0)^kc))
q[, ] <- d1 + (1 - d1) / (1 + (ID[i, j] / id0)^kd * fd[i])       # ms adds NO +0.5 to q
cA[, ] <- cU + (cD - cU) * q[i, j]^g_inf
theta[, ] <- theta0 * (theta1 + (1 - theta1) /
  (1 + fv[i] * (((IVA[i, j] + acq_offset) + IVM[i, j]) / iv0)^kv))
dim(b, phi, q, cA, theta) <- c(n_age, n_het)

## ---- mosquito -> human EIR (via EIR-lag chain) ----------------------------
aIm[] <- a_spp[i] * Im[i]
dim(aIm) <- n_spp
eir_now <- sum(aIm)
# Erlang chain tracking eir_now; eir_lag = last stage
deriv(Xe[1]) <- reir * (eir_now - Xe[1])
deriv(Xe[2:n_eir]) <- reir * (Xe[i - 1] - Xe[i])
dim(Xe) <- n_eir
eir_lag <- Xe[n_eir]
EPS[, ] <- eir_lag * zeta[j] * psi[i]              # expected infectious bites/person/day
# malariasimulation draws the day's bites from a Poisson, scatters them over individuals
# and collects the bitten in a BITSET (biting_process.R): a person bitten several times
# in one timestep counts ONCE, so can be infected at most once per day. Its daily
# infection probability is therefore p = (1 - exp(-EPS)) * b * pev, which SATURATES,
# converted to a hazard by -log(1 - p) (ms prob_to_rate / rate_to_prob are exact
# inverses). The plain product b*EPS instead grows without bound, over-predicting
# infection where EPS approaches 1 -- at seasonal peaks and in high-zeta strata (the
# clinical incidence "peaks too high" bias; severe is spared because young children
# have small psi). bite_dedup = 0 restores the linear form, which is what
# malariaEquilibrium assumes, which halves the drift off the seed without
# q_b, defined with the immunity boosting below, is this same 1 - exp(-EPS):
# the probability of at least one bite in the day. It is reused here rather than
# written out again because odin2 does no common-subexpression elimination, so
# the two spellings cost two exp() per cell per evaluation. The RHS is bound by
# its ~2,600 transcendentals, and this was one in ten of them.
p_inf[, ] <- q_b[i, j] * b[i, j] * pev_factor[i]
FOI[, ] <- bite_dedup * (-log(1 - p_inf[i, j])) +
  (1 - bite_dedup) * (b[i, j] * EPS[i, j] * pev_factor[i])
# Split FOI into the CLINICAL and NON-CLINICAL hazards separately. In the IBM each
# person resolves at most ONE infection outcome per day (competing_hazards.R), so
# the daily clinical count is exactly phi*p*N. Clinical infections leave the
# at-risk pool (to D/Tr) while non-clinical ones do not (S/U -> A, A -> A), so the
# pool depletes only through the clinical route: taking h_c = -log(1 - phi*p)
# integrates to exactly phi*p*N over a day. Using phi*FOI instead over-counts by
# FOI/(1-exp(-FOI)) because the ODE re-exposes within the same day.
# h_c + h_a = FOI exactly, so occupancies and the immunity boosting are unchanged.
h_c[, ] <- bite_dedup * (-log(1 - phi[i, j] * p_inf[i, j])) +
  (1 - bite_dedup) * (phi[i, j] * FOI[i, j])
h_a[, ] <- FOI[i, j] - h_c[i, j]
dim(EPS, FOI, p_inf, h_c, h_a) <- c(n_age, n_het)

## ---- human -> mosquito infectivity (via FOIM-lag chain) -------------------
# TBV reduces onward infectivity, per infection state (state-specific TBA)
inf[, ] <- cD * D[i, j] * tbv_fD[i] + cA[i, j] * A[i, j] * tbv_fA[i] +
  cU * U[i, j] * tbv_fU[i] + cT * (Tr[i, j] + Tr_slow[i, j]) * tbv_fT[i]
infw[, ] <- zeta[j] * psi[i] * inf[i, j]
dim(inf, infw) <- c(n_age, n_het)
# malariasimulation normalises by sum(zeta*psi) over the LIVE population
# (biting_process.R, human_pi), which is the population mean of psi: zeta has
# mean 1 in every age band, because births are spread by het_wt and mortality
# is age-only. That mean moves whenever the age structure does -- under
# set_demography() with time-varying death rates it drifts for decades -- so it
# is taken from the current Npop rather than frozen at the seed's stationary
# structure, as it was. At the seed the two agree to the last bit. Weighted by
# psi alone, not zeta*psi: the Gauss-Hermite nodes carry a mean zeta of 0.9997,
# and folding that in would move the seed by 0.03% for no reason the IBM has.
psi_n[] <- psi[i] * n_g[i]
dim(psi_n) <- n_age
mpsi <- sum(psi_n) / sum(n_g)
inf_sum <- (sum(infw) + sum(infw_v)) / mpsi   # both blocks; one is empty
deriv(Xf[1]) <- rfoim * (inf_sum - Xf[1])
deriv(Xf[2:n_foim]) <- rfoim * (Xf[i - 1] - Xf[i])
dim(Xf) <- n_foim
inf_lag <- if (tl > 0) Xf[n_foim] else inf_sum
foim[] <- a_spp[i] * inf_lag
dim(foim) <- n_spp

## ---- demography bookkeeping (constant population) -------------------------
Ph_tot[, ] <- sum(Ph[i, j, ])
Phc_tot[, ] <- sum(Ph_c[i, j, ])
dim(Ph_tot, Phc_tot) <- c(n_age, n_het)
Npop[, ] <- S[i, j] + D[i, j] + A[i, j] + U[i, j] + Tr[i, j] + Tr_slow[i, j] + Ph_tot[i, j] + Phc_tot[i, j]
deaths[, ] <- mu_age[i] * Npop[i, j]
dim(Npop, deaths) <- c(n_age, n_het)
# Births replace deaths from whichever block holds the population. Under
# falciparum the vivax block is empty, so its deaths are exactly 0 and this is
# sum(deaths) + 0. Each block then receives births only while it is the one in
# use; births_f keeps the falciparum equations below in their original shape.
births <- sum(deaths) + sum(deaths_v)
births_f <- births * pf_on
births_v <- births * pv_on

## ---- human ODEs -----------------------------------------------------------
deriv(S[1, ]) <- births_f * het_wt[j] + rU * U[1, j] + rPk * Ph[1, j, n_ph] +
  rPck * Ph_c[1, j, n_phc] - FOI[1, j] * S[1, j] - re[1] * S[1, j]
deriv(S[2:n_age, ]) <- r_age[i - 1] * S[i - 1, j] + rU * U[i, j] + rPk * Ph[i, j, n_ph] +
  rPck * Ph_c[i, j, n_phc] - FOI[i, j] * S[i, j] - re[i] * S[i, j]

# Treated. malariasimulation assigns each successfully-treated individual to slow
# parasite clearance by a Bernoulli draw with probability `spc` (= artemisinin
# resistance proportion x slow_parasite_clearance_probability) and gives them
# dt_slow, everyone else dt (human_infection.R calculate_successful_treatments ->
# dt_spc_combined). That is a MIXTURE OF TWO EXPONENTIALS, not one exponential at
# the blended mean, so Tr is split into two parallel compartments exactly as the
# IBM splits the treated.
trt_in[, ] <- ft_eff * h_c[i, j] * (S[i, j] + A[i, j] + U[i, j])
dim(trt_in) <- c(n_age, n_het)
deriv(Tr[1, ]) <- (1 - spc) * trt_in[1, j] - rT * Tr[1, j] - re[1] * Tr[1, j]
deriv(Tr[2:n_age, ]) <- r_age[i - 1] * Tr[i - 1, j] +
  (1 - spc) * trt_in[i, j] - rT * Tr[i, j] - re[i] * Tr[i, j]
deriv(Tr_slow[1, ]) <- spc * trt_in[1, j] - rT_slow * Tr_slow[1, j] - re[1] * Tr_slow[1, j]
deriv(Tr_slow[2:n_age, ]) <- r_age[i - 1] * Tr_slow[i - 1, j] +
  spc * trt_in[i, j] - rT_slow * Tr_slow[i, j] - re[i] * Tr_slow[i, j]

# Untreated clinical (D): (1 - ft_eff) includes early-treatment-failure diversions
deriv(D[1, ]) <- (1 - ft_eff) * h_c[1, j] * (S[1, j] + A[1, j] + U[1, j]) -
  rD * D[1, j] - re[1] * D[1, j]
deriv(D[2:n_age, ]) <- r_age[i - 1] * D[i - 1, j] +
  (1 - ft_eff) * h_c[i, j] * (S[i, j] + A[i, j] + U[i, j]) -
  rD * D[i, j] - re[i] * D[i, j]

deriv(A[1, ]) <- h_a[1, j] * (S[1, j] + U[1, j]) -
  h_c[1, j] * A[1, j] + rD * D[1, j] - rA * A[1, j] - re[1] * A[1, j]
deriv(A[2:n_age, ]) <- r_age[i - 1] * A[i - 1, j] +
  h_a[i, j] * (S[i, j] + U[i, j]) -
  h_c[i, j] * A[i, j] + rD * D[i, j] - rA * A[i, j] - re[i] * A[i, j]

deriv(U[1, ]) <- rA * A[1, j] - FOI[1, j] * U[1, j] - rU * U[1, j] - re[1] * U[1, j]
deriv(U[2:n_age, ]) <- r_age[i - 1] * U[i - 1, j] + rA * A[i, j] -
  FOI[i, j] * U[i, j] - rU * U[i, j] - re[i] * U[i, j]

# Post-treatment prophylaxis chain: Tr/Tr_slow feed stage 1, each stage is left at
# rPk = n_ph*rP, the last stage returns to S. Aging runs along every stage.
deriv(Ph[1, , 1]) <- rT * Tr[1, j] + rT_slow * Tr_slow[1, j] - rPk * Ph[1, j, 1] -
  re[1] * Ph[1, j, 1]
deriv(Ph[2:n_age, , 1]) <- r_age[i - 1] * Ph[i - 1, j, 1] + rT * Tr[i, j] +
  rT_slow * Tr_slow[i, j] - rPk * Ph[i, j, 1] - re[i] * Ph[i, j, 1]
deriv(Ph[1, , 2:n_ph]) <- rPk * (Ph[1, j, k - 1] - Ph[1, j, k]) - re[1] * Ph[1, j, k]
deriv(Ph[2:n_age, , 2:n_ph]) <- r_age[i - 1] * Ph[i - 1, j, k] +
  rPk * (Ph[i, j, k - 1] - Ph[i, j, k]) - re[i] * Ph[i, j, k]

# Chemoprevention prophylaxis chain (filled at stage 1 by the MDA/SMC/PMC pulses;
# each stage left at rPck = n_phc*rP_c, the chemoprevention drug's chain rate; the
# last stage returns to S). No FOI, no infectivity.
rP_c <- parameter()
rPck <- rP_c * n_phc
deriv(Ph_c[1, , 1]) <- -rPck * Ph_c[1, j, 1] - re[1] * Ph_c[1, j, 1]
deriv(Ph_c[2:n_age, , 1]) <- r_age[i - 1] * Ph_c[i - 1, j, 1] - rPck * Ph_c[i, j, 1] -
  re[i] * Ph_c[i, j, 1]
deriv(Ph_c[1, , 2:n_phc]) <- rPck * (Ph_c[1, j, k - 1] - Ph_c[1, j, k]) - re[1] * Ph_c[1, j, k]
deriv(Ph_c[2:n_age, , 2:n_phc]) <- r_age[i - 1] * Ph_c[i - 1, j, k] +
  rPck * (Ph_c[i, j, k - 1] - Ph_c[i, j, k]) - re[i] * Ph_c[i, j, k]

## ---- P. vivax human block -------------------------------------------------
# Indexed [age, heterogeneity, hypnozoite level]. malariasimulation carries an
# integer batch count per person, k in 0..kmax, and relapse is a hazard of k*f.
# A mean-field model cannot hold only the mean count: the relapse hazard is
# linear in k, but immunity, detectability and the U->S rate are not, and
# malariaEquilibriumVivax returns every state as [age, het, batch] for the same
# reason. kk[k] is the batch count at level k, so kk = 0 at the first level.
#
# Levels 1..n_bat are the batch counts 0..kmax. Levels after n_bat hold people
# with no batches whose liver stage is still drug-protected after radical cure:
# an Erlang chain from the dose, left for level 1. Liver-stage protection is a
# property of the person, not of their disease state, so it has to sit on the
# dimension every compartment and immunity stock already carries. On those
# levels kk = 0 (no relapse, nothing to clear) and a bite infects but forms no
# batch, which is what ls_prophylaxis does in the IBM. Without a radical-cure
# drug there are no such levels and n_hyp = n_bat.
#
# Boundaries are handled with `if` guards on the index, not with index ranges.
# odin2 ranges are C loops, so an empty one such as 2:1 is harmless; the trouble
# with ranges here is that at n_hyp = 1 -- the falciparum collapse -- the batch-1
# and batch-n_hyp boundary equations would both target the same element, and the
# second reads k - 1 = 0, out of bounds. Guards need one equation per
# compartment rather than six. odin2 warns that it cannot validate accesses such
# as X[i, j, k + 1]; they are safe, because the generated code is a ternary and
# the out-of-range branch is never evaluated. Keep guards out of the arrays that
# compute exp/log -- see pop_floor for what a branch there costs.
kk[] <- if (i <= n_bat) as.numeric(i) - 1 else 0
shmask[] <- if (i < n_bat) 1 else 0       # a bite here moves the person up a batch
lsmask[] <- if (i > n_bat) 1 else 0       # liver-stage-protected levels
dim(kk, shmask, lsmask) <- n_hyp
k_rc <- if (n_hyp > n_bat) n_bat + 1 else 1   # where radical cure puts people
rls_k <- (n_hyp - n_bat) / mean_ls            # per-stage rate of the liver-stage chain

gammal <- parameter()     # per-batch hypnozoite clearance (k -> k-1 at k*gammal)
ff <- parameter()         # per-batch relapse rate (hazard k*ff)
bv <- parameter()         # infection probability per infectious bite; no IB for vivax
philm_min <- parameter(); philm_max <- parameter(); alm50 <- parameter(); klm <- parameter()
phi0_v <- parameter(); phi1_v <- parameter(); ic0_v <- parameter(); kc_v <- parameter()
dpcr_min <- parameter(); dpcr_max <- parameter(); apcr50 <- parameter(); kpcr <- parameter()
# progression D -> A and A -> U at malariasimulation's own rates 1/dd and 1/da;
# they enter the daily competing-hazard draw below. Tr, which no infection
# competes with, leaves at the whole-day rate 1 - exp(-1/dt).
rD_ms <- parameter(); rA_ms <- parameter(); rT_v <- parameter()
cD_v <- parameter(); cA_v <- parameter(); cU_v <- parameter()
d_iaa <- parameter(); d_ica_v <- parameter()
ua_eff <- parameter(); uc_v_eff <- parameter()
PM_v <- parameter()
mat_factor_v <- parameter(); dim(mat_factor_v) <- n_age
n_phv <- parameter(constant = TRUE)   # vivax post-treatment prophylaxis stages
rPk_v <- rP * n_phv                   # rP is the vivax drug's under vivax
# within-cell immunity spread: on/off, and the Gauss-Hermite rule it is read at
spread_on <- parameter()
n_q <- parameter(constant = TRUE)
qz <- parameter(); qw <- parameter()
dim(qz, qw) <- n_q

dim(Sv, Dv, Av, Uv, Trv, Trv_slow, JAv, JCv, KAv, KCv) <- c(n_age_v, n_het_v, n_hyp)
dim(Phv) <- c(n_age_v, n_het_v, n_hyp, n_phv)

## population, and who is exposed to infection. Treated and drug-protected
## people are not: at treatment the IBM's prophylaxis is 1, which fleet
## represents by holding them in Tr and Phv. Everyone else is, including D --
## an infection does not change a D person's disease state, but it does boost
## them and, if it came from a bite, give them a new batch.
Phv_tot[, , ] <- sum(Phv[i, j, k, ])
dim(Phv_tot) <- c(n_age_v, n_het_v, n_hyp)
expo_v[, , ] <- Sv[i, j, k] + Dv[i, j, k] + Av[i, j, k] + Uv[i, j, k]
Npop_v[, , ] <- expo_v[i, j, k] + Trv[i, j, k] + Trv_slow[i, j, k] + Phv_tot[i, j, k]
deaths_v[, , ] <- mu_age[i] * Npop_v[i, j, k]
dim(expo_v, Npop_v, deaths_v) <- c(n_age_v, n_het_v, n_hyp)

# Vivax population by age, mapped onto the shared age index so it can join n_g
# (and through it mpsi, the population-weighted mean biting rate). Under
# falciparum n_age_v = 1 and this is exactly 0 everywhere.
Nv_age[] <- sum(Npop_v[i, , ])
dim(Nv_age) <- n_age_v
nv_g[] <- if (i <= n_age_v) Nv_age[i] else 0
dim(nv_g) <- n_age

## Immunity is held as a STOCK, J = N * I -- the total in a cell -- not the mean
## that the falciparum block holds. The two are the same mathematics. The
## difference is numerical: in mean form every inflow is divided by the
## receiving cell's population, and the vivax ladder has cells that are exactly
## empty (a newborn cannot carry ten batches), so an empty cell with a populated
## neighbour gets a rate near 1e291 and the solver fails on its first step. In
## stock form nothing divides by a population except turning a stock back into
## the mean the Hill functions need, and an empty cell there simply reads 0.
## The falciparum block keeps its mean form: no age group there is ever empty.
IAAv[, , ] <- JAv[i, j, k] / (Npop_v[i, j, k] + pop_floor)
ICAv[, , ] <- JCv[i, j, k] / (Npop_v[i, j, k] + pop_floor)
dim(IAAv, ICAv) <- c(n_age_v, n_het_v, n_hyp)

## maternal immunity (algebraic). Inherited at birth -- at batch 0, since
## newborns carry no hypnozoites -- and decaying with age regardless of the
## batches acquired later, so it is [age, het] and applies at every batch. The
## mother's value is the mean over batches at age 20, which the stocks give
## directly. malariasimulation uses the same pcm and waning rate rm for both.
Nv_ij[, ] <- sum(Npop_v[i, j, ])
ICAv_ij[, ] <- sum(JCv[i, j, ]) / (Nv_ij[i, j] + pop_floor)
IAAv_ij[, ] <- sum(JAv[i, j, ]) / (Nv_ij[i, j] + pop_floor)
ICAv_m20[, ] <- ICAv_ij[i, j] * mask20[i]
IAAv_m20[, ] <- IAAv_ij[i, j] * mask20[i]
dim(Nv_ij, ICAv_ij, IAAv_ij, ICAv_m20, IAAv_m20) <- c(n_age_v, n_het_v)
ICAv20[] <- sum(ICAv_m20[, i])
IAAv20[] <- sum(IAAv_m20[, i])
dim(ICAv20, IAAv20) <- n_het_v
ICMv[, ] <- PM_v * ICAv20[j] * mat_factor_v[i]
IAMv[, ] <- PM_v * IAAv20[j] * mat_factor_v[i]
dim(ICMv, IAMv) <- c(n_age_v, n_het_v)

## immunity -> probability. anti_parasite_immunity() in malariasimulation adds
## NO +0.5 to acquired immunity; clinical_immunity() does, which fleet carries as
## acq_offset exactly as for falciparum.
##
## The IBM evaluates these per person, and people in one cell do not share one
## immunity: someone who carried few batches for years was infected far less
## often than a neighbour with many, and still differs from them after both
## land at the same age, heterogeneity and batch. Measured in the IBM at EIR 20,
## children's immunity has a CV of 0.5-0.6 within a cell, and because the vivax
## curves are steep (kc = 5.4) the mean probability is up to 5x the probability
## at the mean immunity. So each cell carries the second moment too (K = sum of
## I^2, beside J = sum of I) and the curves are averaged over a gamma with that
## mean and variance -- the shape the IBM's within-cell distribution takes. The
## gamma is sampled at Gauss-Hermite nodes through the Wilson-Hilferty map,
## X = m (1 - c/9 + z sqrt(c)/3)^3 with c = CV^2, capped at 9 where the map
## stops being monotone. With spread_on = 0 the variance is ignored and every
## node sits at the mean.
##
## The two immunities share their nodes: a person's IAA and ICA are boosted by
## the same infections, and within a cell they correlate at 0.92-0.97 below age
## 40 in the IBM. So node l is one kind of person -- low in both, or high in
## both -- and anything that needs both at once (LM-detectable AND clinical) is
## averaged over nodes jointly, not multiplied out of two separate averages.
##
## Means are clamped at zero. The vivax ladder holds cells as small as 1e-36
## people, and to choose its first step the solver evaluates the derivatives at a
## trial state one tiny Euler step away -- which can leave such a cell
## fractionally negative, so that J / N is negative and a fractional power of it
## is NaN. max() compiles to one branchless instruction, so this does not
## reintroduce the cost of a branch in a loop that computes powers.
mA[, , ] <- max(IAAv[i, j, k], 0)
mC[, , ] <- max(ICAv[i, j, k], 0)
cv2A[, , ] <- spread_on * min(max(KAv[i, j, k] / (Npop_v[i, j, k] + pop_floor) -
  mA[i, j, k] * mA[i, j, k], 0) / (mA[i, j, k] * mA[i, j, k] + 1e-12), 9)
cv2C[, , ] <- spread_on * min(max(KCv[i, j, k] / (Npop_v[i, j, k] + pop_floor) -
  mC[i, j, k] * mC[i, j, k], 0) / (mC[i, j, k] * mC[i, j, k] + 1e-12), 9)
whaA[, , ] <- 1 - cv2A[i, j, k] / 9
whaC[, , ] <- 1 - cv2C[i, j, k] / 9
whbA[, , ] <- sqrt(cv2A[i, j, k]) / 3
whbC[, , ] <- sqrt(cv2C[i, j, k]) / 3
dim(mA, mC, cv2A, cv2C, whaA, whaC, whbA, whbC) <- c(n_age_v, n_het_v, n_hyp)
uA[, , , ] <- max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0)
uC[, , , ] <- max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0)
xA[, , , ] <- mA[i, j, k] * uA[i, j, k, l] * uA[i, j, k, l] * uA[i, j, k, l] + IAMv[i, j]
xC[, , , ] <- mC[i, j, k] * uC[i, j, k, l] * uC[i, j, k, l] * uC[i, j, k, l] +
  acq_offset + ICMv[i, j]
dim(uA, uC, xA, xC) <- c(n_age_v, n_het_v, n_hyp, n_q)

## infection. Vivax counts every bite, 1 - (1 - b)^n, rather than deduplicating
## them; for Poisson bites at rate EPS that averages to 1 - exp(-b*EPS), whose
## rate is exactly b*EPS. Relapse adds k*ff. Both are summed and THEN reduced by
## PEV, as in calculate_vivax_infections(), so the vaccine blocks relapses too;
## hv is the infection rate that function hands to the competing hazards. Only
## an infection that came from a bite forms a new batch, which is what the
## relative_rates draw in relapse_bite_infection_hazard_resolution() decides;
## sh_v is that bite share wherever a batch can form -- below the top batch,
## which stays at kmax (pmin(k + 1, kmax)), and outside liver-stage protection.
lam_bv[, ] <- bv * EPS[i, j]
dim(lam_bv) <- c(n_age_v, n_het_v)
r_tot_v[, , ] <- lam_bv[i, j] + kk[k] * ff
hv[, , ] <- -log(1 - (1 - exp(-r_tot_v[i, j, k])) * pev_factor[i])
sh_v[, , ] <- shmask[k] * lam_bv[i, j] / (r_tot_v[i, j, k] + pop_floor)
dim(r_tot_v, hv, sh_v) <- c(n_age_v, n_het_v, n_hyp)

## The day, as malariasimulation resolves it. Each person faces infection (hv)
## and their own state's progression (U -> S, A -> U, D -> A) in ONE competing
## hazard draw a day (CompetingHazard$resolve): some event happens with
## probability 1 - exp(-(hv + r)), and it is infection with share hv / (hv + r).
## So infection and recovery each make the other less likely that day, and a
## person is infected at most once. fleet uses those daily probabilities as its
## rates. Without infection the recovery probability is 1 - exp(-r), the
## whole-day rate the falciparum model uses; S has no progression, so its
## infection probability is 1 - exp(-hv). U's progression rate depends on IAA,
## so U's two probabilities are taken per node.
pLMn[, , , ] <- philm_min + (philm_max - philm_min) / (1 + (xA[i, j, k, l] / alm50)^klm)
pDn[, , , ] <- phi0_v * (phi1_v + (1 - phi1_v) / (1 + (xC[i, j, k, l] / ic0_v)^kc_v))
rUn[, , , ] <- 1 / (dpcr_min + (dpcr_max - dpcr_min) / (1 + (xA[i, j, k, l] / apcr50)^kpcr))
eUn[, , , ] <- 1 - exp(-(hv[i, j, k] + rUn[i, j, k, l]))
iUn[, , , ] <- eUn[i, j, k, l] * hv[i, j, k] / (hv[i, j, k] + rUn[i, j, k, l])
## weighted node terms, summed below: LM-detectable, LM-detectable AND
## clinical, clinical; and for U, infected, infected and LM-detectable,
## infected and LM-detectable and clinical, and any event at all
w_lm[, , , ] <- qw[l] * pLMn[i, j, k, l]
w_lmc[, , , ] <- w_lm[i, j, k, l] * pDn[i, j, k, l]
w_c[, , , ] <- qw[l] * pDn[i, j, k, l]
w_iu[, , , ] <- qw[l] * iUn[i, j, k, l]
w_iulm[, , , ] <- w_iu[i, j, k, l] * pLMn[i, j, k, l]
w_iulmc[, , , ] <- w_iulm[i, j, k, l] * pDn[i, j, k, l]
w_eu[, , , ] <- qw[l] * eUn[i, j, k, l]
dim(pLMn, pDn, rUn, eUn, iUn, w_lm, w_lmc, w_c, w_iu, w_iulm, w_iulmc, w_eu) <-
  c(n_age_v, n_het_v, n_hyp, n_q)
s_lm[, , ] <- sum(w_lm[i, j, k, ])
s_lmc[, , ] <- sum(w_lmc[i, j, k, ])
s_c[, , ] <- sum(w_c[i, j, k, ])
iU_v[, , ] <- sum(w_iu[i, j, k, ])
s_iulm[, , ] <- sum(w_iulm[i, j, k, ])
s_iulmc[, , ] <- sum(w_iulmc[i, j, k, ])
recU_v[, , ] <- sum(w_eu[i, j, k, ]) - iU_v[i, j, k]
iS_v[, , ] <- 1 - exp(-hv[i, j, k])
eA_v[, , ] <- 1 - exp(-(hv[i, j, k] + rA_ms))
iA_v[, , ] <- eA_v[i, j, k] * hv[i, j, k] / (hv[i, j, k] + rA_ms)
recA_v[, , ] <- eA_v[i, j, k] - iA_v[i, j, k]
eD_v[, , ] <- 1 - exp(-(hv[i, j, k] + rD_ms))
iD_v[, , ] <- eD_v[i, j, k] * hv[i, j, k] / (hv[i, j, k] + rD_ms)
recD_v[, , ] <- eD_v[i, j, k] - iD_v[i, j, k]
dim(s_lm, s_lmc, s_c, iU_v, s_iulm, s_iulmc, recU_v, iS_v, eA_v, iA_v, recA_v,
    eD_v, iD_v, recD_v) <- c(n_age_v, n_het_v, n_hyp)

## where infections go. vivax_infection_outcome_process(): S and U pass through
## LM-detectability first (the rest of S go to U, the rest of U stay U), then
## LM-detectable S and U, and ALL infected A, go through the clinical draw; an
## infected D stays D. to*_n are those fluxes, each routed at its SOURCE -- the
## probabilities belong to the person's own cell -- and then landing in place if
## it was a relapse, or one batch up if it was a bite: a bite changes disease
## state and batch at once.
toU_n[, , ] <- iS_v[i, j, k] * (1 - s_lm[i, j, k]) * Sv[i, j, k] +
  (iU_v[i, j, k] - s_iulm[i, j, k]) * Uv[i, j, k]
toA_n[, , ] <- iS_v[i, j, k] * (s_lm[i, j, k] - s_lmc[i, j, k]) * Sv[i, j, k] +
  (s_iulm[i, j, k] - s_iulmc[i, j, k]) * Uv[i, j, k] +
  iA_v[i, j, k] * (1 - s_c[i, j, k]) * Av[i, j, k]
toC_n[, , ] <- iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] + s_iulmc[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k]
toD_n[, , ] <- iD_v[i, j, k] * Dv[i, j, k]
inf_n[, , ] <- iS_v[i, j, k] * Sv[i, j, k] + iU_v[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * Av[i, j, k] + iD_v[i, j, k] * Dv[i, j, k]
dim(toU_n, toA_n, toC_n, toD_n, inf_n) <- c(n_age_v, n_het_v, n_hyp)

arr_U[, , ] <- (1 - sh_v[i, j, k]) * toU_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * toU_n[i, j, k - 1] else 0)
arr_A[, , ] <- (1 - sh_v[i, j, k]) * toA_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * toA_n[i, j, k - 1] else 0)
arr_C[, , ] <- (1 - sh_v[i, j, k]) * toC_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * toC_n[i, j, k - 1] else 0)
arr_D[, , ] <- (1 - sh_v[i, j, k]) * toD_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * toD_n[i, j, k - 1] else 0)
dim(arr_U, arr_A, arr_C, arr_D) <- c(n_age_v, n_het_v, n_hyp)

## Treatment splits the clinical flux exactly as for falciparum: ft_eff to
## treatment (split fast/slow clearance by spc), the rest -- untreated and early
## treatment failures -- to D.
trt_v[, , ] <- ft_eff * arr_C[i, j, k]
## Radical cure. Of everyone who is treated, drug_hypnozoite_efficacy lose every
## batch -- whether or not the blood stage cleared, since calculate_treated()
## draws the two independently -- and start liver-stage protection, so they move
## to level k_rc. rcT_v are cleared in the blood too (to Tr), rcD_v are not (to
## D). The treated flux is summed over levels and delivered to k_rc. rcT_v
## carries ft_eff's (1 - etf): an early treatment failure leaves the blood stage
## uncleared, but its liver-stage draw still stands, so it falls to rcD_v.
rcT_v[, , ] <- ft * eff_hyp * (1 - etf) * arr_C[i, j, k]
rcD_v[, , ] <- ft * hyp_mix * arr_C[i, j, k] - rcT_v[i, j, k]
dim(trt_v, rcT_v, rcD_v) <- c(n_age_v, n_het_v, n_hyp)
rcT_tot[, ] <- sum(rcT_v[i, j, ])
rcD_tot[, ] <- sum(rcD_v[i, j, ])
dim(rcT_tot, rcD_tot) <- c(n_age_v, n_het_v)

## Ageing, death and batch clearance are common to every compartment; newborns
## enter S at batch 0. Batch clearance moves everyone down one level whatever
## their disease state -- hypnozoites die in the liver regardless of the blood.
## The liver-stage chain is common to every compartment too, since protection
## runs from the dose whatever happens to the blood stage: each protected level
## passes on at rls_k and the last one returns to level 1.
deriv(Sv[, , ]) <-
  (if (i > 1) r_age[i - 1] * Sv[i - 1, j, k] else (if (k == 1) births_v * het_wt[j] else 0)) -
  re[i] * Sv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Sv[i, j, k + 1] else 0) - gammal * kk[k] * Sv[i, j, k] +
  (if (k > n_bat + 1) rls_k * Sv[i, j, k - 1] else 0) - rls_k * lsmask[k] * Sv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Sv[i, j, n_hyp] else 0) +
  recU_v[i, j, k] * Uv[i, j, k] + rPk_v * Phv[i, j, k, n_phv] -
  iS_v[i, j, k] * Sv[i, j, k]
deriv(Uv[, , ]) <-
  (if (i > 1) r_age[i - 1] * Uv[i - 1, j, k] else 0) - re[i] * Uv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Uv[i, j, k + 1] else 0) - gammal * kk[k] * Uv[i, j, k] +
  (if (k > n_bat + 1) rls_k * Uv[i, j, k - 1] else 0) - rls_k * lsmask[k] * Uv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Uv[i, j, n_hyp] else 0) +
  recA_v[i, j, k] * Av[i, j, k] - recU_v[i, j, k] * Uv[i, j, k] +
  arr_U[i, j, k] - iU_v[i, j, k] * Uv[i, j, k]
deriv(Av[, , ]) <-
  (if (i > 1) r_age[i - 1] * Av[i - 1, j, k] else 0) - re[i] * Av[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Av[i, j, k + 1] else 0) - gammal * kk[k] * Av[i, j, k] +
  (if (k > n_bat + 1) rls_k * Av[i, j, k - 1] else 0) - rls_k * lsmask[k] * Av[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Av[i, j, n_hyp] else 0) +
  recD_v[i, j, k] * Dv[i, j, k] - recA_v[i, j, k] * Av[i, j, k] +
  arr_A[i, j, k] - iA_v[i, j, k] * Av[i, j, k]
deriv(Dv[, , ]) <-
  (if (i > 1) r_age[i - 1] * Dv[i - 1, j, k] else 0) - re[i] * Dv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Dv[i, j, k + 1] else 0) - gammal * kk[k] * Dv[i, j, k] +
  (if (k > n_bat + 1) rls_k * Dv[i, j, k - 1] else 0) - rls_k * lsmask[k] * Dv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Dv[i, j, n_hyp] else 0) -
  recD_v[i, j, k] * Dv[i, j, k] +
  arr_D[i, j, k] + arr_C[i, j, k] - trt_v[i, j, k] - iD_v[i, j, k] * Dv[i, j, k] -
  rcD_v[i, j, k] + (if (k == k_rc) rcD_tot[i, j] else 0)
deriv(Trv[, , ]) <-
  (if (i > 1) r_age[i - 1] * Trv[i - 1, j, k] else 0) - re[i] * Trv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Trv[i, j, k + 1] else 0) - gammal * kk[k] * Trv[i, j, k] +
  (if (k > n_bat + 1) rls_k * Trv[i, j, k - 1] else 0) - rls_k * lsmask[k] * Trv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Trv[i, j, n_hyp] else 0) -
  rT_v * Trv[i, j, k] +
  (1 - spc) * (trt_v[i, j, k] - rcT_v[i, j, k] + (if (k == k_rc) rcT_tot[i, j] else 0))
deriv(Trv_slow[, , ]) <-
  (if (i > 1) r_age[i - 1] * Trv_slow[i - 1, j, k] else 0) - re[i] * Trv_slow[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * Trv_slow[i, j, k + 1] else 0) -
  gammal * kk[k] * Trv_slow[i, j, k] +
  (if (k > n_bat + 1) rls_k * Trv_slow[i, j, k - 1] else 0) - rls_k * lsmask[k] * Trv_slow[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Trv_slow[i, j, n_hyp] else 0) -
  rT_slow * Trv_slow[i, j, k] +
  spc * (trt_v[i, j, k] - rcT_v[i, j, k] + (if (k == k_rc) rcT_tot[i, j] else 0))

## Post-treatment prophylaxis: an Erlang chain like falciparum's Ph, but carrying
## the batch dimension, because people keep losing batches while drug-protected
## and must return to the right one. Capped at a few stages for vivax (see
## build_inputs), which keeps the mean and widens the spread.
deriv(Phv[, , , ]) <-
  (if (l == 1) rT_v * Trv[i, j, k] + rT_slow * Trv_slow[i, j, k] else rPk_v * Phv[i, j, k, l - 1]) -
  rPk_v * Phv[i, j, k, l] +
  (if (i > 1) r_age[i - 1] * Phv[i - 1, j, k, l] else 0) - re[i] * Phv[i, j, k, l] +
  (if (k < n_hyp) gammal * kk[k + 1] * Phv[i, j, k + 1, l] else 0) -
  gammal * kk[k] * Phv[i, j, k, l] +
  (if (k > n_bat + 1) rls_k * Phv[i, j, k - 1, l] else 0) - rls_k * lsmask[k] * Phv[i, j, k, l] +
  (if (k == 1) rls_k * lsmask[n_hyp] * Phv[i, j, n_hyp, l] else 0)

## Immunity dynamics, in stock form (see the note at IAAv). Each flow of people
## carries its cell's mean with it, which in stock terms is just the stock times
## the per-capita rate: ageing in r_age * J, ageing out and death re * J, batch
## clearance kk * gammal * J, the liver-stage chain rls_k * J, newborns nothing
## acquired. A bite moves the infected up a batch, so the stock carried up is
## sh_v * (inf_n / N) * J, where inf_n / N is a daily probability and cannot
## blow up.
##
## Boosting is by infection, with the same renewal approximation to the integer
## refractory window as falciparum, p / (p * u + 1) for a person infected with
## daily probability p -- which now differs by state, since infection competes
## with each state's own progression. A boost belongs to the cell the person
## LANDS in: a relapse leaves them where they are, a bite moves them up a batch,
## so boosts are split by the shift share and shifted exactly as the people are.
## t*_ = 1 / (p u + 1): boosts per infection in each state, for each window
tSA[, , ] <- 1 / (iS_v[i, j, k] * ua_eff + 1)
tUA[, , ] <- 1 / (iU_v[i, j, k] * ua_eff + 1)
tAA[, , ] <- 1 / (iA_v[i, j, k] * ua_eff + 1)
tDA[, , ] <- 1 / (iD_v[i, j, k] * ua_eff + 1)
tSC[, , ] <- 1 / (iS_v[i, j, k] * uc_v_eff + 1)
tUC[, , ] <- 1 / (iU_v[i, j, k] * uc_v_eff + 1)
tAC[, , ] <- 1 / (iA_v[i, j, k] * uc_v_eff + 1)
tDC[, , ] <- 1 / (iD_v[i, j, k] * uc_v_eff + 1)
bstA_n[, , ] <- iS_v[i, j, k] * tSA[i, j, k] * Sv[i, j, k] + iU_v[i, j, k] * tUA[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * tAA[i, j, k] * Av[i, j, k] + iD_v[i, j, k] * tDA[i, j, k] * Dv[i, j, k]
bstC_n[, , ] <- iS_v[i, j, k] * tSC[i, j, k] * Sv[i, j, k] + iU_v[i, j, k] * tUC[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * tAC[i, j, k] * Av[i, j, k] + iD_v[i, j, k] * tDC[i, j, k] * Dv[i, j, k]
finf_v[, , ] <- inf_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)
dim(tSA, tUA, tAA, tDA, tSC, tUC, tAC, tDC, bstA_n, bstC_n, finf_v) <- c(n_age_v, n_het_v, n_hyp)

bstA_in[, , ] <- bstA_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bstA_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
bstC_in[, , ] <- bstC_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bstC_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
upA_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * JAv[i, j, k]
upC_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * JCv[i, j, k]
dim(bstA_in, bstC_in, upA_v, upC_v) <- c(n_age_v, n_het_v, n_hyp)

## Radical cure takes its people's immunity to k_rc with them: what they carried
## from the cell they were infected in -- in place, or a batch down for a bite,
## exactly as upA_v carries it -- plus the boosts that infection gave them, each
## clinical source counted at its own 1 / (p * u + 1) boosts per infection. fC_v
## is the clinical share of a cell, so as with finf_v nothing divides a flow by a
## population that could be empty.
fC_v[, , ] <- toC_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)
cBA_v[, , ] <- iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * tSA[i, j, k] +
  s_iulmc[i, j, k] * Uv[i, j, k] * tUA[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * tAA[i, j, k]
cBC_v[, , ] <- iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * tSC[i, j, k] +
  s_iulmc[i, j, k] * Uv[i, j, k] * tUC[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * tAC[i, j, k]
cJA_n[, , ] <- fC_v[i, j, k] * JAv[i, j, k] + cBA_v[i, j, k]
cJC_n[, , ] <- fC_v[i, j, k] * JCv[i, j, k] + cBC_v[i, j, k]
dim(fC_v, cBA_v, cBC_v, cJA_n, cJC_n) <- c(n_age_v, n_het_v, n_hyp)
rcJA_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cJA_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cJA_n[i, j, k - 1] else 0))
rcJC_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cJC_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cJC_n[i, j, k - 1] else 0))
dim(rcJA_v, rcJC_v) <- c(n_age_v, n_het_v, n_hyp)
rcJA_tot[, ] <- sum(rcJA_v[i, j, ])
rcJC_tot[, ] <- sum(rcJC_v[i, j, ])
dim(rcJA_tot, rcJC_tot) <- c(n_age_v, n_het_v)

deriv(JAv[, , ]) <- bstA_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * JAv[i - 1, j, k] else 0) - re[i] * JAv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * JAv[i, j, k + 1] else 0) - gammal * kk[k] * JAv[i, j, k] +
  (if (k > n_bat + 1) rls_k * JAv[i, j, k - 1] else 0) - rls_k * lsmask[k] * JAv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * JAv[i, j, n_hyp] else 0) +
  (if (k > 1) upA_v[i, j, k - 1] else 0) - upA_v[i, j, k] -
  rcJA_v[i, j, k] + (if (k == k_rc) rcJA_tot[i, j] else 0) -
  JAv[i, j, k] / d_iaa
deriv(JCv[, , ]) <- bstC_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * JCv[i - 1, j, k] else 0) - re[i] * JCv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * JCv[i, j, k + 1] else 0) - gammal * kk[k] * JCv[i, j, k] +
  (if (k > n_bat + 1) rls_k * JCv[i, j, k - 1] else 0) - rls_k * lsmask[k] * JCv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * JCv[i, j, n_hyp] else 0) +
  (if (k > 1) upC_v[i, j, k - 1] else 0) - upC_v[i, j, k] -
  rcJC_v[i, j, k] + (if (k == k_rc) rcJC_tot[i, j] else 0) -
  JCv[i, j, k] / d_ica_v

## Second moments, K = sum of I^2 over a cell. Every flow of people carries its
## share of K exactly as it carries its share of J -- which is exact, not an
## approximation, because nothing that moves a person between cells or infects
## them depends on their own immunity once their cell is known. What differs is
## what a boost and a day of decay do to I^2: a boost takes I to I + 1, adding
## 2 I + 1, and decay by exp(-1/d) a day shrinks I^2 at rate 2 / d. Boosted people
## are a random draw of their cell, so b boosts carry b * J / N of immunity and
## add 2 b J / N to K; b / N is a daily probability, so nothing blows up.
##
## The rest of a boost's contribution is the spread it injects, and that is NOT
## one unit per boost. The refractory window makes a person's boosts a renewal
## process -- ceil(u) - 1 blocked days, then a geometric wait at daily
## probability p -- far more regular than independent events: over a long
## stretch its count varies by the gap's squared CV times its mean, the Fano
## factor (1 - p) / (p u + 1)^2. For IAA's 44-day window that is ~0.06, so
## counting each boost as a unit of variance, as a Markov jump would, overstates
## the IBM's within-cell spread of IAA by 20-45%. The injection is therefore
## boosts x Fano, i.e. p (1 - p) t^3 per person with t = 1 / (p u + 1).
bvA_n[, , ] <-
  iS_v[i, j, k] * (1 - iS_v[i, j, k]) * tSA[i, j, k] * tSA[i, j, k] * tSA[i, j, k] * Sv[i, j, k] +
  iU_v[i, j, k] * (1 - iU_v[i, j, k]) * tUA[i, j, k] * tUA[i, j, k] * tUA[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * (1 - iA_v[i, j, k]) * tAA[i, j, k] * tAA[i, j, k] * tAA[i, j, k] * Av[i, j, k] +
  iD_v[i, j, k] * (1 - iD_v[i, j, k]) * tDA[i, j, k] * tDA[i, j, k] * tDA[i, j, k] * Dv[i, j, k]
bvC_n[, , ] <-
  iS_v[i, j, k] * (1 - iS_v[i, j, k]) * tSC[i, j, k] * tSC[i, j, k] * tSC[i, j, k] * Sv[i, j, k] +
  iU_v[i, j, k] * (1 - iU_v[i, j, k]) * tUC[i, j, k] * tUC[i, j, k] * tUC[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * (1 - iA_v[i, j, k]) * tAC[i, j, k] * tAC[i, j, k] * tAC[i, j, k] * Av[i, j, k] +
  iD_v[i, j, k] * (1 - iD_v[i, j, k]) * tDC[i, j, k] * tDC[i, j, k] * tDC[i, j, k] * Dv[i, j, k]
dim(bvA_n, bvC_n) <- c(n_age_v, n_het_v, n_hyp)
bKA_n[, , ] <- 2 * bstA_n[i, j, k] / (Npop_v[i, j, k] + pop_floor) * JAv[i, j, k] +
  bvA_n[i, j, k]
bKC_n[, , ] <- 2 * bstC_n[i, j, k] / (Npop_v[i, j, k] + pop_floor) * JCv[i, j, k] +
  bvC_n[i, j, k]
bKA_in[, , ] <- bKA_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bKA_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
bKC_in[, , ] <- bKC_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bKC_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
upKA_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * KAv[i, j, k]
upKC_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * KCv[i, j, k]
## radical cure moves the spread its people's boosts injected with them, by
## the same Fano rule, from each clinical source
cvA_v[, , ] <-
  iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * (1 - iS_v[i, j, k]) * tSA[i, j, k] * tSA[i, j, k] * tSA[i, j, k] +
  s_iulmc[i, j, k] * Uv[i, j, k] * (1 - iU_v[i, j, k]) * tUA[i, j, k] * tUA[i, j, k] * tUA[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * (1 - iA_v[i, j, k]) * tAA[i, j, k] * tAA[i, j, k] * tAA[i, j, k]
cvC_v[, , ] <-
  iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * (1 - iS_v[i, j, k]) * tSC[i, j, k] * tSC[i, j, k] * tSC[i, j, k] +
  s_iulmc[i, j, k] * Uv[i, j, k] * (1 - iU_v[i, j, k]) * tUC[i, j, k] * tUC[i, j, k] * tUC[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * (1 - iA_v[i, j, k]) * tAC[i, j, k] * tAC[i, j, k] * tAC[i, j, k]
dim(cvA_v, cvC_v) <- c(n_age_v, n_het_v, n_hyp)
cKA_n[, , ] <- fC_v[i, j, k] * KAv[i, j, k] +
  2 * cBA_v[i, j, k] / (Npop_v[i, j, k] + pop_floor) * JAv[i, j, k] + cvA_v[i, j, k]
cKC_n[, , ] <- fC_v[i, j, k] * KCv[i, j, k] +
  2 * cBC_v[i, j, k] / (Npop_v[i, j, k] + pop_floor) * JCv[i, j, k] + cvC_v[i, j, k]
rcKA_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cKA_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cKA_n[i, j, k - 1] else 0))
rcKC_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cKC_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cKC_n[i, j, k - 1] else 0))
dim(bKA_n, bKC_n, bKA_in, bKC_in, upKA_v, upKC_v, cKA_n, cKC_n, rcKA_v, rcKC_v) <-
  c(n_age_v, n_het_v, n_hyp)
rcKA_tot[, ] <- sum(rcKA_v[i, j, ])
rcKC_tot[, ] <- sum(rcKC_v[i, j, ])
dim(rcKA_tot, rcKC_tot) <- c(n_age_v, n_het_v)

deriv(KAv[, , ]) <- bKA_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * KAv[i - 1, j, k] else 0) - re[i] * KAv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * KAv[i, j, k + 1] else 0) - gammal * kk[k] * KAv[i, j, k] +
  (if (k > n_bat + 1) rls_k * KAv[i, j, k - 1] else 0) - rls_k * lsmask[k] * KAv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * KAv[i, j, n_hyp] else 0) +
  (if (k > 1) upKA_v[i, j, k - 1] else 0) - upKA_v[i, j, k] -
  rcKA_v[i, j, k] + (if (k == k_rc) rcKA_tot[i, j] else 0) -
  2 * KAv[i, j, k] / d_iaa
deriv(KCv[, , ]) <- bKC_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * KCv[i - 1, j, k] else 0) - re[i] * KCv[i, j, k] +
  (if (k < n_hyp) gammal * kk[k + 1] * KCv[i, j, k + 1] else 0) - gammal * kk[k] * KCv[i, j, k] +
  (if (k > n_bat + 1) rls_k * KCv[i, j, k - 1] else 0) - rls_k * lsmask[k] * KCv[i, j, k] +
  (if (k == 1) rls_k * lsmask[n_hyp] * KCv[i, j, n_hyp] else 0) +
  (if (k > 1) upKC_v[i, j, k - 1] else 0) - upKC_v[i, j, k] -
  rcKC_v[i, j, k] + (if (k == k_rc) rcKC_tot[i, j] else 0) -
  2 * KCv[i, j, k] / d_ica_v

## onward infectivity. A is LM-detectable by definition for vivax and carries
## the constant ca, where falciparum's asymptomatic infectivity depends on ID.
inf_v[, , ] <- cD_v * Dv[i, j, k] * tbv_fD[i] + cA_v * Av[i, j, k] * tbv_fA[i] +
  cU_v * Uv[i, j, k] * tbv_fU[i] + cT * (Trv[i, j, k] + Trv_slow[i, j, k]) * tbv_fT[i]
infw_v[, , ] <- zeta[j] * psi[i] * inf_v[i, j, k]
dim(inf_v, infw_v) <- c(n_age_v, n_het_v, n_hyp)

## ---- immunity ODEs (aging uses re[i]; boundary I0 = 0) --------------------
# Immunity boosting, matching the IBM's refractory renewal process exactly.
# malariasimulation boosts on a per-DAY EVENT with an integer refractory window:
# boost_immunity() fires only if (timestep - last_boosted) >= u, and timestep is an
# integer, so the effective wait is ceil(u) days and the person then boosts on the
# first subsequent day carrying an event -- geometric with mean 1/q. Mean inter-boost
# gap = (ceil(u) - 1) + 1/q, i.e. boost rate = q / (q * u_eff + 1) with
# u_eff = ceil(u) - 1. The event probability q is the DEDUPLICATED per-day one:
#   IB          -> bitten at least once, q = 1 - exp(-EPS)  (bitten_humans is a Bitset)
#   ICA/ID/IVA  -> infected that day,    q = 1 - exp(-FOI)  (one outcome per person/day)
# fleet previously used the raw rates EPS/FOI, over-boosting where exposure is high.
# ICA/ID/IVA are boosted only for individuals eligible to be INFECTED (ms restricts
# source_humans to S/A/U), so their boost is scaled by the at-risk fraction; IB is
# boosted for everyone bitten, whatever their state, so it is not scaled.
q_b[, ] <- 1 - exp(-EPS[i, j])
q_f[, ] <- 1 - exp(-FOI[i, j])
at_risk[, ] <- (S[i, j] + A[i, j] + U[i, j]) / (Npop[i, j] + pop_floor)
bst_b[, ] <- q_b[i, j] / (q_b[i, j] * ub_eff + 1)
bst_c[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * uc_eff + 1)
bst_d[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * ud_eff + 1)
bst_v[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * uv_eff + 1)
dim(q_b, q_f, at_risk, bst_b, bst_c, bst_d, bst_v) <- c(n_age, n_het)
# IB/ICA/ID/IVA hold the MEAN immunity of the people in a cell, not a stock, so
# ageing is a difference ain*(I[i-1] - I[i]) rather than a flux -- do not
# "fix" it to r_age[i-1]*I[i-1] by analogy with the disease compartments.
# Writing J = N*I for the extensive total and differentiating I = J/N,
#     dI_i/dt = (r_age[i-1] N_{i-1,j} / N_{i,j}) (I_{i-1} - I_i) + boost - decay,
# and mu_age cancels out of it entirely: the dead carry the cell mean, so death
# leaves the mean alone. At i = 1 the inflow is the newborns, births*het_wt[j],
# arriving with I = 0. That coefficient, ain, is the per-capita inflow into the
# cell, and it is taken from the live Npop.
# An earlier version used re[i] = r_age[i] + mu_age[i] instead. The two are
# equal at a stationary age structure, where r_age[i-1] N_{i-1} = re[i] N_i,
# and agree to 2e-17 at the seed -- but under set_demography() with
# time-varying death rates the structure drifts for decades, and re[i], built
# from the NEW death rate while N still had the OLD shape, was the wrong
# coefficient for the whole transient. The site files carry exactly such
# series (yearly rates, 2000 onward).
# births_f rather than births -- the same number under falciparum, and the
# correct one by construction, since it is what actually enters S. The
# pop_floor keeps an empty block finite; see its definition.
ain[1, ] <- births_f * het_wt[j] / (Npop[1, j] + pop_floor)
ain[2:n_age, ] <- r_age[i - 1] * Npop[i - 1, j] / (Npop[i, j] + pop_floor)
dim(ain) <- c(n_age, n_het)
deriv(IB[1, ]) <- bst_b[1, j] - IB[1, j] / d_ib - ain[1, j] * IB[1, j]
deriv(IB[2:n_age, ]) <- bst_b[i, j] - IB[i, j] / d_ib +
  ain[i, j] * (IB[i - 1, j] - IB[i, j])
deriv(ICA[1, ]) <- bst_c[1, j] - ICA[1, j] / d_ica - ain[1, j] * ICA[1, j]
deriv(ICA[2:n_age, ]) <- bst_c[i, j] - ICA[i, j] / d_ica +
  ain[i, j] * (ICA[i - 1, j] - ICA[i, j])
deriv(ID[1, ]) <- bst_d[1, j] - ID[1, j] / d_id - ain[1, j] * ID[1, j]
deriv(ID[2:n_age, ]) <- bst_d[i, j] - ID[i, j] / d_id +
  ain[i, j] * (ID[i - 1, j] - ID[i, j])
deriv(IVA[1, ]) <- bst_v[1, j] - IVA[1, j] / d_iva - ain[1, j] * IVA[1, j]
deriv(IVA[2:n_age, ]) <- bst_v[i, j] - IVA[i, j] / d_iva +
  ain[i, j] * (IVA[i - 1, j] - IVA[i, j])

## ---- mosquito ODEs (per species) ------------------------------------------
Mtot[] <- Sm[i] + Em_tot[i] + Im[i]
nL[] <- ME[i] + ML[i]
Em_tot[] <- Em_inc[i]                # total incubating per species (ms E state)
dim(Mtot, nL, Em_tot) <- n_spp

deriv(ME[]) <- beta_eff[i] * Mtot[i] - ME[i] / del - ME[i] * me * (1 + nL[i] / Kcap[i])
deriv(ML[]) <- ME[i] / del - ML[i] / dl - ML[i] * ml * (1 + mosq_gamma * nL[i] / Kcap[i])
deriv(MP[]) <- ML[i] / dl - MP[i] / dpl - MP[i] * mup
dim(ME, ML, MP) <- n_spp

# The 0.5 is the SEX RATIO, not a loss: emerging pupae are half female, and only
# females are modelled from here on (as in malariasimulation). MP therefore
# empties at MP/dpl while Sm gains half of it; the males are intentionally
# untracked, not leaked.
deriv(Sm[]) <- 0.5 * MP[i] / dpl - Sm[i] * foim[i] - Sm[i] * mum[i]
# EIP. malariasimulation uses a FIXED dem-day delay and applies the incubation
# survival exp(-mu*dem) at the exit, with the CURRENT mu
# (src/adult_mosquito_eqs.cpp: incubation_survival = exp(-model.mu * model.tau);
#  dE/dt = S*foim - lagged_incubating*incubation_survival - E*mu;
#  dI/dt = lagged_incubating*incubation_survival - I*mu).
# fleet substitutes an Erlang chain for the delay (odin2's delay() is unreliable
# here), but the chain must be LOSS-FREE: putting death in every stage gives
# through-survival (reip/(reip+mum))^n_eip instead of exp(-mum*dem), which is
# ~4% too high at baseline mum and worsens as vector control raises mum (the
# static total_M compensation is exact only at the seed). So Xi is a pure delay
# of Sm*foim, survival is applied once at the exit, and the incubating stock
# Em_tot follows malariasimulation's own E equation.
deriv(Xi[, 1]) <- reip * (Sm[i] * foim[i] - Xi[i, 1])
deriv(Xi[, 2:n_eip]) <- reip * (Xi[i, j - 1] - Xi[i, j])
eip_surv[] <- exp(-mum[i] * dem)
dim(eip_surv) <- n_spp
matured[] <- Xi[i, n_eip] * eip_surv[i]        # survived the incubation period
dim(matured) <- n_spp
# Xi and Em_inc are not two tallies of the same thing: Xi carries the infection
# RATE through the delay, Em_inc is the STOCK of mosquitoes currently incubating.
# Losing `matured` (proportional to Xi, not to Em_inc) is correct, because for a
# fixed delay with continuous mortality the exact stock is
#     E(t) = int_0^dem F(t-a) e^{-mum a} da,
# and differentiating that gives exactly
#     dE/dt = F(t) - F(t-dem) e^{-mum dem} - mum E(t),
# which is this line, with Xi[n_eip] standing in for F(t-dem). At steady state it
# yields E = (F/mum)(1 - e^{-mum dem}); the seeded value matches that to 2e-16.
# The only approximation is the n_eip-stage Erlang standing in for a fixed delay.
deriv(Em_inc[]) <- Sm[i] * foim[i] - matured[i] - mum[i] * Em_inc[i]
deriv(Im[]) <- matured[i] - mum[i] * Im[i]
dim(Sm, Im) <- n_spp
dim(Xi) <- c(n_spp, n_eip)
dim(Em_inc) <- n_spp

## ---- initial conditions ---------------------------------------------------
S0 <- parameter(); D0 <- parameter(); A0 <- parameter()
U0 <- parameter(); Tr0 <- parameter(); Ph0 <- parameter(); Phc0 <- parameter()
IB_init <- parameter(); ICA_init <- parameter(); ID_init <- parameter(); IVA_init <- parameter()
dim(S0, D0, A0, U0, Tr0, IB_init, ICA_init, ID_init, IVA_init) <- c(n_age, n_het)
dim(Ph0) <- c(n_age, n_het, n_ph)
dim(Phc0) <- c(n_age, n_het, n_phc)
ME0 <- parameter(); ML0 <- parameter(); MP0 <- parameter(); Sm0 <- parameter(); Im0 <- parameter()
dim(ME0, ML0, MP0, Sm0, Im0) <- n_spp
Xi0 <- parameter(); dim(Xi0) <- c(n_spp, n_eip)
Em_inc0 <- parameter(); dim(Em_inc0) <- n_spp
Xe0 <- parameter(); Xf0 <- parameter()   # scalar equilibrium lag values

initial(S[, ]) <- S0[i, j]
initial(D[, ]) <- D0[i, j]
initial(A[, ]) <- A0[i, j]
initial(U[, ]) <- U0[i, j]
initial(Tr[, ]) <- (1 - spc0) * Tr0[i, j]
initial(Tr_slow[, ]) <- spc0 * Tr0[i, j]
initial(Ph[, , ]) <- Ph0[i, j, k]
initial(Ph_c[, , ]) <- Phc0[i, j, k]
initial(IB[, ]) <- IB_init[i, j]
initial(ICA[, ]) <- ICA_init[i, j]
initial(ID[, ]) <- ID_init[i, j]
initial(IVA[, ]) <- IVA_init[i, j]
initial(ME[]) <- ME0[i]
initial(ML[]) <- ML0[i]
initial(MP[]) <- MP0[i]
initial(Sm[]) <- Sm0[i]
initial(Xi[, ]) <- Xi0[i, j]
initial(Em_inc[]) <- Em_inc0[i]
initial(Im[]) <- Im0[i]
initial(Xe[]) <- Xe0
initial(Xf[]) <- Xf0

# vivax block: seeded from malariaEquilibriumVivax under vivax, all zero (in a
# single cell) under falciparum
Sv0 <- parameter(); Dv0 <- parameter(); Av0 <- parameter()
Uv0 <- parameter(); Trv0 <- parameter()
dim(Sv0, Dv0, Av0, Uv0, Trv0) <- c(n_age_v, n_het_v, n_hyp)
initial(Sv[, , ]) <- Sv0[i, j, k]
initial(Dv[, , ]) <- Dv0[i, j, k]
initial(Av[, , ]) <- Av0[i, j, k]
initial(Uv[, , ]) <- Uv0[i, j, k]
initial(Trv[, , ]) <- (1 - spc0) * Trv0[i, j, k]
initial(Trv_slow[, , ]) <- spc0 * Trv0[i, j, k]
Phv0 <- parameter(); dim(Phv0) <- c(n_age_v, n_het_v, n_hyp, n_phv)
JAv0 <- parameter(); JCv0 <- parameter(); KAv0 <- parameter(); KCv0 <- parameter()
dim(JAv0, JCv0, KAv0, KCv0) <- c(n_age_v, n_het_v, n_hyp)
initial(Phv[, , , ]) <- Phv0[i, j, k, l]
initial(JAv[, , ]) <- JAv0[i, j, k]
initial(JCv[, , ]) <- JCv0[i, j, k]
initial(KAv[, , ]) <- KAv0[i, j, k]
initial(KCv[, , ]) <- KCv0[i, j, k]

dim(S, D, A, U, Tr, Tr_slow) <- c(n_age, n_het)
dim(Ph) <- c(n_age, n_het, n_ph)
dim(Ph_c) <- c(n_age, n_het, n_phc)
dim(IB, ICA, ID, IVA) <- c(n_age, n_het)

## ---- outputs (per age group; aggregated to bands in R) --------------------
# All as fractions of the total human population (x human_population -> counts).
#
# These are RATES that R integrates over a day, so the rate to use for a
# compartment depends on how fast that compartment drains. The IBM counts events
# per person per day, i.e. the DEDUPLICATED probability p_inf (or phi*p_inf), not
# the hazard.
#   S and U drain at the full FOI, so int_0^1 FOI*X exp(-FOI t) dt = X*(1-exp(-FOI))
#     = X*p_inf exactly. The hazard form is already the IBM's count for them.
#   A does NOT: a sub-clinical re-infection of an A leaves them in A, so A drains
#     only at h_c and stays ~flat over the day. int FOI*A dt is then ~ FOI*A, which
#     over-counts by FOI/p_inf = -log(1-p_inf)/p_inf. So A takes p_inf directly.
# Clinical is unaffected: it is counted with h_c = -log(1 - phi*p_inf), and
# int h_c*A dt = A*(1-exp(-h_c)) = A*phi*p_inf, the IBM's clinical count.
# Severe is drawn from the SAME infected set as clinical in the IBM
# (update_severe_disease takes infected_humans), so it is theta x the all-infection
# count, not theta x the clinical one.
# Measured against 3 IBM replicates of 10,000 people: before this correction
# all-age n_inc_* ran +2.6% (EIR 20) and +4.2% (EIR 50) above the IBM median, both
# outside its replicate spread; after, +0.6% and +0.4%. The bias was concentrated
# in adults, where A is a large share of the at-risk pool (under-5 n_inc_* was only
# +0.3% / +1.7%), which is the signature of an A-only defect.
clin_inc_a[, ] <- h_c[i, j] * (S[i, j] + A[i, j] + U[i, j])
sev_inc_a[, ] <- theta[i, j] * (FOI[i, j] * (S[i, j] + U[i, j]) + p_inf[i, j] * A[i, j])
inc_a[, ] <- FOI[i, j] * (S[i, j] + U[i, j]) + p_inf[i, j] * A[i, j]   # all new infections
detlm[, ] <- D[i, j] + Tr[i, j] + Tr_slow[i, j] + q[i, j] * A[i, j]
# PCR follows the malariasimulation IBM convention (all D/Tr/A/U count as
# PCR-positive), NOT malariaEquilibrium's sub-patent-weighted pos_PCR
# (D+Tr+A*q^aA+U*q^aU). We target the IBM, and this matches it to <0.3%.
detpcr[, ] <- D[i, j] + Tr[i, j] + Tr_slow[i, j] + A[i, j] + U[i, j]
dim(clin_inc_a, sev_inc_a, inc_a, detlm, detpcr) <- c(n_age, n_het)

## The vivax block's share of the same outputs (plus its own two), by age.
## vivax A is LM-detectable by definition, so LM counts D, Tr and all of A, as
## create_prevalence_renderer() does for vivax. The infection fluxes are the
## daily probabilities of the competing-hazard draw, so they are the IBM's daily
## counts as they stand. Relapses are the relapse share of that count,
## relapse_rates / infection_rates as in calculate_vivax_infections(). Severe
## disease does not exist for vivax.
hypmask[] <- if (i > 1 && i <= n_bat) 1 else 0     # levels holding at least one batch
dim(hypmask) <- n_hyp
detlm_v[, , ] <- Dv[i, j, k] + Trv[i, j, k] + Trv_slow[i, j, k] + Av[i, j, k]
detpcr_v[, , ] <- detlm_v[i, j, k] + Uv[i, j, k]
rel_v[, , ] <- inf_n[i, j, k] * kk[k] * ff / (r_tot_v[i, j, k] + pop_floor)
hyp_v[, , ] <- hypmask[k] * Npop_v[i, j, k]
dim(detlm_v, detpcr_v, rel_v, hyp_v) <- c(n_age_v, n_het_v, n_hyp)
dlm_va[] <- sum(detlm_v[i, , ]); dpcr_va[] <- sum(detpcr_v[i, , ])
clin_va[] <- sum(arr_C[i, , ]); inc_va[] <- sum(inf_n[i, , ])
rel_va[] <- sum(rel_v[i, , ]); hyp_va[] <- sum(hyp_v[i, , ])
S_va[] <- sum(Sv[i, , ]); D_va[] <- sum(Dv[i, , ]); A_va[] <- sum(Av[i, , ])
U_va[] <- sum(Uv[i, , ]); Tr_va[] <- sum(Trv[i, , ]) + sum(Trv_slow[i, , ])
Ph_va[] <- sum(Phv_tot[i, , ])
dim(dlm_va, dpcr_va, clin_va, inc_va, rel_va, hyp_va) <- n_age_v
dim(S_va, D_va, A_va, U_va, Tr_va, Ph_va) <- n_age_v

n_g[] <- sum(Npop[i, ]) + nv_g[i]    # both blocks; the one not in use is empty
det_lm_g[] <- sum(detlm[i, ]) + (if (i <= n_age_v) dlm_va[i] else 0)
det_pcr_g[] <- sum(detpcr[i, ]) + (if (i <= n_age_v) dpcr_va[i] else 0)
clin_g[] <- sum(clin_inc_a[i, ]) + (if (i <= n_age_v) clin_va[i] else 0)
sev_g[] <- sum(sev_inc_a[i, ])
inc_g[] <- sum(inc_a[i, ]) + (if (i <= n_age_v) inc_va[i] else 0)
relapse_g[] <- if (i <= n_age_v) rel_va[i] else 0
hyp_g[] <- if (i <= n_age_v) hyp_va[i] else 0
S_g[] <- sum(S[i, ]) + (if (i <= n_age_v) S_va[i] else 0)
D_g[] <- sum(D[i, ]) + (if (i <= n_age_v) D_va[i] else 0)
A_g[] <- sum(A[i, ]) + (if (i <= n_age_v) A_va[i] else 0)
U_g[] <- sum(U[i, ]) + (if (i <= n_age_v) U_va[i] else 0)
Tr_g[] <- sum(Tr[i, ]) + sum(Tr_slow[i, ]) + (if (i <= n_age_v) Tr_va[i] else 0)
Ph_g[] <- sum(Ph_tot[i, ]) + sum(Phc_tot[i, ]) + (if (i <= n_age_v) Ph_va[i] else 0)
dim(n_g, det_lm_g, det_pcr_g, clin_g, sev_g, inc_g, relapse_g, hyp_g) <- n_age
dim(S_g, D_g, A_g, U_g, Tr_g, Ph_g) <- n_age

output(n_g) <- TRUE
output(det_lm_g) <- TRUE
output(det_pcr_g) <- TRUE
output(clin_g) <- TRUE
output(sev_g) <- TRUE
output(inc_g) <- TRUE
output(relapse_g) <- TRUE
output(hyp_g) <- TRUE
output(S_g) <- TRUE; output(D_g) <- TRUE; output(A_g) <- TRUE
output(U_g) <- TRUE; output(Tr_g) <- TRUE; output(Ph_g) <- TRUE
output(EIR_yr) <- eir_now * 365
output(FOIM) <- foim[1]
output(ft_out) <- ft
