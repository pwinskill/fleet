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

## ---- grid / demography (data) ---------------------------------------------
r_age <- parameter(); psi <- parameter()
age_mid <- parameter(); mask20 <- parameter()
icm_factor <- parameter(); ivm_factor <- parameter(); mean_psi <- parameter()
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
rfoim <- n_foim / tl   # per-stage rate of the FOIM-lag chain
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
p_inf[, ] <- (1 - exp(-EPS[i, j])) * b[i, j] * pev_factor[i]
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
inf_sum <- sum(infw) / mean_psi
deriv(Xf[1]) <- rfoim * (inf_sum - Xf[1])
deriv(Xf[2:n_foim]) <- rfoim * (Xf[i - 1] - Xf[i])
dim(Xf) <- n_foim
inf_lag <- Xf[n_foim]
foim[] <- a_spp[i] * inf_lag
dim(foim) <- n_spp

## ---- demography bookkeeping (constant population) -------------------------
Ph_tot[, ] <- sum(Ph[i, j, ])
Phc_tot[, ] <- sum(Ph_c[i, j, ])
dim(Ph_tot, Phc_tot) <- c(n_age, n_het)
Npop[, ] <- S[i, j] + D[i, j] + A[i, j] + U[i, j] + Tr[i, j] + Tr_slow[i, j] + Ph_tot[i, j] + Phc_tot[i, j]
deaths[, ] <- mu_age[i] * Npop[i, j]
dim(Npop, deaths) <- c(n_age, n_het)
births <- sum(deaths)

## ---- human ODEs -----------------------------------------------------------
deriv(S[1, ]) <- births * het_wt[j] + rU * U[1, j] + rPk * Ph[1, j, n_ph] +
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
at_risk[, ] <- (S[i, j] + A[i, j] + U[i, j]) / Npop[i, j]
bst_b[, ] <- q_b[i, j] / (q_b[i, j] * ub_eff + 1)
bst_c[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * uc_eff + 1)
bst_d[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * ud_eff + 1)
bst_v[, ] <- at_risk[i, j] * q_f[i, j] / (q_f[i, j] * uv_eff + 1)
dim(q_b, q_f, at_risk, bst_b, bst_c, bst_d, bst_v) <- c(n_age, n_het)
deriv(IB[1, ]) <- bst_b[1, j] - IB[1, j] / d_ib - re[1] * IB[1, j]
deriv(IB[2:n_age, ]) <- bst_b[i, j] - IB[i, j] / d_ib +
  re[i] * (IB[i - 1, j] - IB[i, j])
deriv(ICA[1, ]) <- bst_c[1, j] - ICA[1, j] / d_ica - re[1] * ICA[1, j]
deriv(ICA[2:n_age, ]) <- bst_c[i, j] - ICA[i, j] / d_ica +
  re[i] * (ICA[i - 1, j] - ICA[i, j])
deriv(ID[1, ]) <- bst_d[1, j] - ID[1, j] / d_id - re[1] * ID[1, j]
deriv(ID[2:n_age, ]) <- bst_d[i, j] - ID[i, j] / d_id +
  re[i] * (ID[i - 1, j] - ID[i, j])
deriv(IVA[1, ]) <- bst_v[1, j] - IVA[1, j] / d_iva - re[1] * IVA[1, j]
deriv(IVA[2:n_age, ]) <- bst_v[i, j] - IVA[i, j] / d_iva +
  re[i] * (IVA[i - 1, j] - IVA[i, j])

## ---- mosquito ODEs (per species) ------------------------------------------
Mtot[] <- Sm[i] + Em_tot[i] + Im[i]
nL[] <- ME[i] + ML[i]
Em_tot[] <- Em_inc[i]                # total incubating per species (ms E state)
dim(Mtot, nL, Em_tot) <- n_spp

deriv(ME[]) <- beta_eff[i] * Mtot[i] - ME[i] / del - ME[i] * me * (1 + nL[i] / Kcap[i])
deriv(ML[]) <- ME[i] / del - ML[i] / dl - ML[i] * ml * (1 + mosq_gamma * nL[i] / Kcap[i])
deriv(MP[]) <- ML[i] / dl - MP[i] / dpl - MP[i] * mup
dim(ME, ML, MP) <- n_spp

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

n_g[] <- sum(Npop[i, ])
det_lm_g[] <- sum(detlm[i, ])
det_pcr_g[] <- sum(detpcr[i, ])
clin_g[] <- sum(clin_inc_a[i, ])
sev_g[] <- sum(sev_inc_a[i, ])
inc_g[] <- sum(inc_a[i, ])
S_g[] <- sum(S[i, ]); D_g[] <- sum(D[i, ]); A_g[] <- sum(A[i, ])
U_g[] <- sum(U[i, ]); Tr_g[] <- sum(Tr[i, ]) + sum(Tr_slow[i, ]); Ph_g[] <- sum(Ph_tot[i, ]) + sum(Phc_tot[i, ])
dim(n_g, det_lm_g, det_pcr_g, clin_g, sev_g, inc_g) <- n_age
dim(S_g, D_g, A_g, U_g, Tr_g, Ph_g) <- n_age

output(n_g) <- TRUE
output(det_lm_g) <- TRUE
output(det_pcr_g) <- TRUE
output(clin_g) <- TRUE
output(sev_g) <- TRUE
output(inc_g) <- TRUE
output(S_g) <- TRUE; output(D_g) <- TRUE; output(A_g) <- TRUE
output(U_g) <- TRUE; output(Tr_g) <- TRUE; output(Ph_g) <- TRUE
output(EIR_yr) <- eir_now * 365
output(FOIM) <- foim[1]
output(ft_out) <- ft
