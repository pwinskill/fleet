# =============================================================================
# fleet: mean-field twin of malariasimulation (P. falciparum), on its daily clock
# =============================================================================
# The population is advanced one day at a time, in the order malariasimulation
# resolves a day. Every process reads the state at the start of the day and
# every change lands at the end of it, as the IBM queues its updates:
#
#   1. Immunity decays by exp(-1/d). A person boosted today skips the decay:
#      the IBM queues the decay first and the boost (start-of-day value + 1)
#      after it, and the later update overwrites the earlier.
#   2. Biting. The EIR saved `de` days ago reaches each age x heterogeneity
#      stratum; the day's bites are deduplicated (a person bitten several times
#      counts once), IB is boosted for everyone bitten, and a bitten person in
#      S, A or U is infected with probability b x PEV.
#   3. Infection and each state's progression are ONE competing draw a day
#      (CompetingHazard$resolve): an infection pre-empts that day's A -> U or
#      U -> S. An infected person is clinical with probability phi, and a
#      clinical case is treated with probability ft; every infection counts
#      towards severe incidence with probability theta.
#   4. The mosquito model is stepped across the day with the day's biting rate,
#      mortality, force of infection and larval carrying capacity held fixed.
#      The EIR, human-infectivity and incubation lags are the IBM's own delay
#      lines, not smoothed stand-ins.
#   5. Deaths are replaced by births into the youngest age group.
#
# Chemoprevention rounds (MDA, SMC, PMC) are applied between days, in R.
#
# All state is expressed as fractions of the human population.
# =============================================================================

pop_floor <- 1e-300   # keeps a division by an empty stratum finite

## ---- dimensions -----------------------------------------------------------
n_age <- parameter(constant = TRUE)
n_het <- parameter(constant = TRUE)
n_spp <- parameter(constant = TRUE)
n_ph <- parameter(constant = TRUE)    # post-treatment prophylaxis chain stages
n_phc <- parameter(constant = TRUE)   # chemoprevention prophylaxis chain stages
n_phct <- parameter(constant = TRUE)  # its counterpart after a chemoprevention Tr stay
n_sub <- parameter(constant = TRUE)   # mosquito sub-steps per day
# delay lines, in days: the EIR lag de, the infectivity lag delay_gam and the
# incubation lag ceiling(dem) - 1. A fractional lag is read as the IBM's
# LaggedValue reads it, linearly between the two saved days either side.
de_floor <- parameter(type = "integer", constant = TRUE)   # floor(de)
de_frac <- parameter()                                     # de - floor(de)
n_eirh <- parameter(constant = TRUE)                       # de_floor + 1 saved days
fl_floor <- parameter(type = "integer", constant = TRUE)   # floor(delay_gam)
fl_frac <- parameter()
n_infh <- parameter(constant = TRUE)                       # fl_floor + 1 saved days
n_tau <- parameter(type = "integer", constant = TRUE)      # ceiling(dem) - 1
n_inch <- parameter(constant = TRUE)                       # max(1, n_tau) saved days

## ---- grid / demography (data) ---------------------------------------------
r_age <- parameter(); psi <- parameter()
age_mid <- parameter(); mask20 <- parameter()
icm_factor <- parameter(); ivm_factor <- parameter()
dim(r_age, psi, age_mid, mask20, icm_factor, ivm_factor) <- n_age
# age-specific mortality, time-varying for custom demography (constant otherwise):
# supplied as [n_age, n_mut] with time as the LAST dimension
n_mut <- parameter(constant = TRUE)
mu_age_t <- parameter(); dim(mu_age_t) <- n_mut
mu_age_z <- parameter(); dim(mu_age_z) <- c(n_age, n_mut)
mu_age <- interpolate(mu_age_t, mu_age_z, "constant"); dim(mu_age) <- n_age

## ---- heterogeneity (data) -------------------------------------------------
zeta <- parameter(); het_wt <- parameter()
dim(zeta, het_wt) <- n_het

## ---- human constants --------------------------------------------------------
# rA, rD, rU, rT are the IBM's whole-day exit probabilities, 1 - exp(-1/d)
rA <- parameter(); rD <- parameter(); rU <- parameter(); rT <- parameter()
d_ib <- parameter(); d_ica <- parameter(); d_id <- parameter(); d_iva <- parameter()
# refractory windows in whole days, ceiling(u) - 1: boost_immunity() fires when
# (timestep - last_boosted) >= u, and timestep is an integer
ub_eff <- parameter(); uc_eff <- parameter(); ud_eff <- parameter(); uv_eff <- parameter()
# the same windows as array indices, at least 1, for the refractory closure below
n_w <- parameter(type = "integer", constant = TRUE)   # max(1, uc_eff, ud_eff, uv_eff)
wc <- parameter(type = "integer", constant = TRUE)
wd <- parameter(type = "integer", constant = TRUE)
wv <- parameter(type = "integer", constant = TRUE)
# malariasimulation adds +0.5 to POSITIVE acquired immunity before the b/phi/theta
# Hill calls, but NOT before q and NOT on the maternal term. That is a
# PER-INDIVIDUAL detail: in the mean field it does NOT translate to +0.5 on the
# stratum mean, and an A/B against the IBM ensemble mean favours 0, the default.
# Set 0.5 to reproduce the IBM's literal per-individual Hill functions instead.
acq_offset <- parameter(0)
# 1 = the IBM's deduplicated bites (a person bitten several times in a day is
# infected at most once); 0 = independent bites, each infecting with probability b
bite_dedup <- parameter(1)
b0 <- parameter(); b1 <- parameter(); ib0 <- parameter(); kb <- parameter()
phi0 <- parameter(); phi1 <- parameter(); ic0 <- parameter(); kc <- parameter()
d1 <- parameter(); id0 <- parameter(); kd <- parameter()
fd0 <- parameter(); ad0 <- parameter(); gd <- parameter()
theta0 <- parameter(); theta1 <- parameter(); iv0 <- parameter(); kv <- parameter()
fv0 <- parameter(); av <- parameter(); gammav <- parameter()
cD <- parameter(); cU <- parameter(); g_inf <- parameter()
PM <- parameter(); PVM <- parameter()

## ---- clinical treatment, drugs and resistance (time-varying steps) -----------
n_ftt <- parameter(constant = TRUE)
ft_times <- parameter(); ft_vals <- parameter()
dim(ft_times, ft_vals) <- n_ftt
ft <- interpolate(ft_times, ft_vals, "constant")
# drug-linked efficacy, treated infectivity (cT) and prophylaxis (rP), so a
# first-line drug switch changes the mix, not just total coverage
n_dmix <- parameter(constant = TRUE)
dmix_times <- parameter(); dim(dmix_times) <- n_dmix
cT_vals <- parameter(); drug_eff_vals <- parameter(); rP_vals <- parameter()
dim(cT_vals, drug_eff_vals, rP_vals) <- n_dmix
cT <- interpolate(dmix_times, cT_vals, "constant")
drug_eff <- interpolate(dmix_times, drug_eff_vals, "constant")
rP <- interpolate(dmix_times, rP_vals, "constant")
# etf = fraction of would-be-treated that fail early (-> stay clinical, D);
# spc = fraction of the treated with slow parasite clearance (longer Tr)
res_times <- parameter(); etf_vals <- parameter(); spc_vals <- parameter()
n_rest <- parameter(constant = TRUE)
dim(res_times, etf_vals, spc_vals) <- n_rest
etf <- interpolate(res_times, etf_vals, "constant")
spc <- interpolate(res_times, spc_vals, "constant")
rT_slow <- parameter()          # whole-day exit probability for slow clearance
spc0 <- parameter(0)            # slow-clearance fraction at the seed (splits Tr)
ft_eff <- ft * drug_eff * (1 - etf)   # coverage x efficacy x (1 - early failure)
# Prophylaxis. The IBM applies each treated person's Weibull protection
# W(t - t_drug) to their infection probability, so a treated cohort's mean
# protection at lag t is W(t). fleet carries it as a chain of fully protected
# stages whose exit-time distribution matches W's mean and variance: a stage is
# left with probability p a day, so it lasts a geometric time, and build_inputs
# picks the stage count for that (see erlang_stages()). rP and rP_c are
# 1 / the chain means, so the per-stage probabilities are n x rP.
rPk <- rP * n_ph
rP_c <- parameter()
rPck <- rP_c * n_phc
# A chemoprevention round sends the clinical and the LM-detectable asymptomatic
# it clears through Tr first, as the IBM does (update_mass_drug_admin()): they
# stay detectable and infectious at their infectivity x drug_rel_c for the Tr
# stay, then join a chain carrying the protection left after it.
rP_ct <- parameter()
rPctk <- rP_ct * n_phct
rT_slow_c <- parameter()        # slow clearance for the chemoprevention drug(s)

## ---- vaccines: per-age infection (PEV) and infectivity (TBV) multipliers ----
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

## ---- derived per-age quantities -------------------------------------------
re[] <- r_age[i] + mu_age[i]          # fraction leaving an age group each day
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

## ---- the day's bites ----------------------------------------------------------
# lagged_eir$get(t - de), per species: a * Im saved de days ago
eir_used[] <- (1 - de_frac) * (if (de_floor == 0) aIm[i] else eir_hist[i, de_floor]) +
  de_frac * eir_hist[i, de_floor + 1]
dim(eir_used) <- n_spp
eir_lag <- sum(eir_used)
# Expected infectious bites per person today. The IBM draws the day's bites as a
# Poisson total of EIR x mean(psi) and shares them in proportion to zeta x psi
# over the live population (biting_process.R), so a person's share is
# zeta psi mean(psi) / mean(zeta psi); mean(zeta) is the quadrature's
# sum(w zeta), 0.99972 at 5 nodes, not 1.
EPS[, ] <- eir_lag * zeta[j] * psi[i] * mpsi / mzp
q_b[, ] <- 1 - exp(-EPS[i, j])               # bitten at least once today
# Daily infection probability for a person in S, A or U. The IBM collects the
# day's bitten in a bitset, so repeated bites count once: a person is bitten at
# least once with probability q_b, and a bitten person is infected with
# probability b x PEV. bite_dedup = 0 lets every bite infect independently
# instead, a hazard b x EPS x PEV.
pbit[, ] <- b[i, j] * pev_factor[i]
hnd[, ] <- b[i, j] * EPS[i, j] * pev_factor[i]
p_inf[, ] <- bite_dedup * q_b[i, j] * pbit[i, j] + (1 - bite_dedup) * (1 - exp(-hnd[i, j]))
dim(EPS, q_b, pbit, hnd, p_inf) <- c(n_age, n_het)

## ---- one competing draw per person ---------------------------------------------
# CompetingHazard$resolve: a person facing infection at hazard h and progression
# at hazard hA has an event with probability 1 - exp(-(h + hA)), and it is the
# infection with probability h / (h + hA). In the IBM only the day's bitten
# carry an infection hazard, -log(1 - b x PEV); the unbitten face progression
# alone. S has no progression, D and Tr no infection.
hA <- -log(1 - rA)
hU <- -log(1 - rU)
hb[, ] <- -log(1 - pbit[i, j])
eAb[, ] <- 1 - (1 - pbit[i, j]) * (1 - rA)      # bitten A: infected, or A -> U
eUb[, ] <- 1 - (1 - pbit[i, j]) * (1 - rU)
eAn[, ] <- 1 - exp(-(hnd[i, j] + hA))            # every A, bites independent
eUn[, ] <- 1 - exp(-(hnd[i, j] + hU))
iA[, ] <- bite_dedup * q_b[i, j] * eAb[i, j] * hb[i, j] / (hb[i, j] + hA) +
  (1 - bite_dedup) * eAn[i, j] * hnd[i, j] / (hnd[i, j] + hA)
pAU[, ] <- bite_dedup * (q_b[i, j] * eAb[i, j] * hA / (hb[i, j] + hA) + (1 - q_b[i, j]) * rA) +
  (1 - bite_dedup) * eAn[i, j] * hA / (hnd[i, j] + hA)
iU[, ] <- bite_dedup * q_b[i, j] * eUb[i, j] * hb[i, j] / (hb[i, j] + hU) +
  (1 - bite_dedup) * eUn[i, j] * hnd[i, j] / (hnd[i, j] + hU)
pUS[, ] <- bite_dedup * (q_b[i, j] * eUb[i, j] * hU / (hb[i, j] + hU) + (1 - q_b[i, j]) * rU) +
  (1 - bite_dedup) * eUn[i, j] * hU / (hnd[i, j] + hU)
infS[, ] <- p_inf[i, j] * S[i, j]
infA[, ] <- iA[i, j] * A[i, j]
infU[, ] <- iU[i, j] * U[i, j]
inf_tot[, ] <- infS[i, j] + infA[i, j] + infU[i, j]
clin_n[, ] <- phi[i, j] * inf_tot[i, j]
trt_in[, ] <- ft_eff * clin_n[i, j]
dim(hb, eAb, eUb, eAn, eUn, iA, pAU, iU, pUS, infS, infA, infU, inf_tot, clin_n, trt_in) <-
  c(n_age, n_het)

## ---- population ------------------------------------------------------------------
Ph_tot[, ] <- sum(Ph[i, j, ])
Phc_tot[, ] <- sum(Ph_c[i, j, ]) + sum(Ph_ct[i, j, ])
Trc_tot[, ] <- Tr_c[i, j] + Tr_cs[i, j]
Npop[, ] <- S[i, j] + D[i, j] + A[i, j] + U[i, j] + Tr[i, j] + Tr_slow[i, j] +
  Trc_tot[i, j] + Ph_tot[i, j] + Phc_tot[i, j]
deaths[, ] <- mu_age[i] * Npop[i, j]
dim(Ph_tot, Phc_tot, Trc_tot, Npop, deaths) <- c(n_age, n_het)
births <- sum(deaths)
n_g_now[] <- sum(Npop[i, ])
dim(n_g_now) <- n_age

## ---- the day's transitions -----------------------------------------------------------
# Each compartment keeps what does not leave it, gains what ages in from the
# group below (births at the bottom) and loses what ages out or dies. A clinical
# case goes to Tr if treated successfully and to D otherwise; a non-clinical
# infection goes to A (an A stays A).
update(S[, ]) <- S[i, j] + (if (i > 1) r_age[i - 1] * S[i - 1, j] else births * het_wt[j]) -
  re[i] * S[i, j] + pUS[i, j] * U[i, j] + rPk * Ph[i, j, n_ph] + rPck * Ph_c[i, j, n_phc] +
  rPctk * Ph_ct[i, j, n_phct] - infS[i, j]
update(D[, ]) <- D[i, j] + (if (i > 1) r_age[i - 1] * D[i - 1, j] else 0) - re[i] * D[i, j] +
  (1 - ft_eff) * clin_n[i, j] - rD * D[i, j]
update(A[, ]) <- A[i, j] + (if (i > 1) r_age[i - 1] * A[i - 1, j] else 0) - re[i] * A[i, j] +
  (1 - phi[i, j]) * (infS[i, j] + infU[i, j]) - phi[i, j] * infA[i, j] + rD * D[i, j] -
  pAU[i, j] * A[i, j]
update(U[, ]) <- U[i, j] + (if (i > 1) r_age[i - 1] * U[i - 1, j] else 0) - re[i] * U[i, j] +
  pAU[i, j] * A[i, j] - infU[i, j] - pUS[i, j] * U[i, j]
# The IBM gives each successfully treated person slow parasite clearance with
# probability spc and dt_slow, everyone else dt: a mixture of two sojourns, so
# two parallel compartments.
update(Tr[, ]) <- Tr[i, j] + (if (i > 1) r_age[i - 1] * Tr[i - 1, j] else 0) - re[i] * Tr[i, j] +
  (1 - spc) * trt_in[i, j] - rT * Tr[i, j]
update(Tr_slow[, ]) <- Tr_slow[i, j] + (if (i > 1) r_age[i - 1] * Tr_slow[i - 1, j] else 0) -
  re[i] * Tr_slow[i, j] + spc * trt_in[i, j] - rT_slow * Tr_slow[i, j]
# post-treatment prophylaxis: fed by the treated, the last stage returns to S
update(Ph[, , ]) <- Ph[i, j, k] + (if (i > 1) r_age[i - 1] * Ph[i - 1, j, k] else 0) -
  re[i] * Ph[i, j, k] +
  (if (k == 1) rT * Tr[i, j] + rT_slow * Tr_slow[i, j] else rPk * Ph[i, j, k - 1]) -
  rPk * Ph[i, j, k]
# chemoprevention prophylaxis: filled at stage 1 by the MDA/SMC/PMC pulses
update(Ph_c[, , ]) <- Ph_c[i, j, k] + (if (i > 1) r_age[i - 1] * Ph_c[i - 1, j, k] else 0) -
  re[i] * Ph_c[i, j, k] + (if (k > 1) rPck * Ph_c[i, j, k - 1] else 0) - rPck * Ph_c[i, j, k]
# The chemoprevention treated phase, filled by the pulses: fast and slow
# clearance, and beside each the infectivity its occupants carry (J = the sum of
# their infectivities x drug_rel_c), which leaves with them.
update(Tr_c[, ]) <- Tr_c[i, j] + (if (i > 1) r_age[i - 1] * Tr_c[i - 1, j] else 0) -
  re[i] * Tr_c[i, j] - rT * Tr_c[i, j]
update(Tr_cs[, ]) <- Tr_cs[i, j] + (if (i > 1) r_age[i - 1] * Tr_cs[i - 1, j] else 0) -
  re[i] * Tr_cs[i, j] - rT_slow_c * Tr_cs[i, j]
update(J_c[, ]) <- J_c[i, j] + (if (i > 1) r_age[i - 1] * J_c[i - 1, j] else 0) -
  re[i] * J_c[i, j] - rT * J_c[i, j]
update(J_cs[, ]) <- J_cs[i, j] + (if (i > 1) r_age[i - 1] * J_cs[i - 1, j] else 0) -
  re[i] * J_cs[i, j] - rT_slow_c * J_cs[i, j]
update(Ph_ct[, , ]) <- Ph_ct[i, j, k] + (if (i > 1) r_age[i - 1] * Ph_ct[i - 1, j, k] else 0) -
  re[i] * Ph_ct[i, j, k] +
  (if (k == 1) rT * Tr_c[i, j] + rT_slow_c * Tr_cs[i, j] else rPctk * Ph_ct[i, j, k - 1]) -
  rPctk * Ph_ct[i, j, k]
dim(S, D, A, U, Tr, Tr_slow, Tr_c, Tr_cs, J_c, J_cs) <- c(n_age, n_het)
dim(Ph) <- c(n_age, n_het, n_ph)
dim(Ph_c) <- c(n_age, n_het, n_phc)
dim(Ph_ct) <- c(n_age, n_het, n_phct)

## ---- immunity: boosts, decay, and transport with ageing ------------------------------
# The chance a person is boosted today. The IBM boosts on a day carrying an
# event once the refractory window since the last boost has passed. IB's event
# is being bitten, whatever the person's state, so its boosts are a renewal
# process with mean gap u_eff + 1/p for a daily bite probability p: a boost
# probability p / (p u_eff + 1).
pbB[, ] <- q_b[i, j] / (q_b[i, j] * ub_eff + 1)
# ICA, ID and IVA's event is being infected, and an infection also moves the
# person: someone boosted is in D, A or Tr the next day, and whether they are
# still refractory when next infected depends on where they have got to. A
# refractory window is days and the sojourns in A and U months, so this is a
# quasi-steady closure: Qr[d, c, .] is the chance that someone boosted d days ago
# is in A, U or D (c = 1..3) today, from today's rates, and K_c, its sum over
# the window, is the refractory days spent in c per boost. With B boosts a day
# the refractory stock in c is K_c B; those people are infected at c's own rate
# without boosting, so
#   B = infections / (1 + K_A iA + K_U iU).
# The renewal form, the same p / (p u_eff + 1) in every state, treats someone in
# S or U as though reinfected as often as someone in A, and under-boosts. The
# closure leaves out the refractory who are back in S within the window: from U
# they are of order 1e-3 of the refractory days in A, and a treated case returns
# to S only through the prophylaxis chain, fully protected until then.
qAA[, ] <- 1 - pAU[i, j] - phi[i, j] * iA[i, j] - mu_age[i]   # refractory A stays A
qUU[, ] <- 1 - iU[i, j] - pUS[i, j] - mu_age[i]
qDD[] <- 1 - rD - mu_age[i]
qcl[, ] <- phi[i, j] * (1 - ft_eff)                            # clinical, untreated: D
dim(qAA, qUU, qcl) <- c(n_age, n_het)
dim(qDD) <- n_age
# The window index comes first so each day's step runs over contiguous memory,
# the component branch outside the loops over the cells.
Qr[1, , , ] <- if (j == 1) 1 - phi[k, l] else if (j == 3) qcl[k, l] else 0
Qr[2:n_w, , , ] <- (if (j == 1)
  Qr[i - 1, 1, k, l] * qAA[k, l] + Qr[i - 1, 3, k, l] * rD +
    (1 - phi[k, l]) * Qr[i - 1, 2, k, l] * iU[k, l]
  else if (j == 2)
  Qr[i - 1, 2, k, l] * qUU[k, l] + Qr[i - 1, 1, k, l] * pAU[k, l]
  else
  Qr[i - 1, 3, k, l] * qDD[k] +
    qcl[k, l] * (Qr[i - 1, 1, k, l] * iA[k, l] + Qr[i - 1, 2, k, l] * iU[k, l]))
Kr[1, , , ] <- Qr[1, j, k, l]
Kr[2:n_w, , , ] <- Kr[i - 1, j, k, l] + Qr[i, j, k, l]
dim(Qr) <- c(n_w, 3, n_age, n_het)
dim(Kr) <- c(n_w, 2, n_age, n_het)
# min(u_eff, 1) switches the correction off for a zero window
pbC[, ] <- inf_tot[i, j] / (1 + min(uc_eff, 1) * (Kr[wc, 1, i, j] * iA[i, j] +
  Kr[wc, 2, i, j] * iU[i, j])) / (Npop[i, j] + pop_floor)
pbD[, ] <- inf_tot[i, j] / (1 + min(ud_eff, 1) * (Kr[wd, 1, i, j] * iA[i, j] +
  Kr[wd, 2, i, j] * iU[i, j])) / (Npop[i, j] + pop_floor)
pbV[, ] <- inf_tot[i, j] / (1 + min(uv_eff, 1) * (Kr[wv, 1, i, j] * iA[i, j] +
  Kr[wv, 2, i, j] * iU[i, j])) / (Npop[i, j] + pop_floor)
# IB/ICA/ID/IVA hold the MEAN immunity of the people in a cell, not a stock, so
# ageing moves it as a difference ain * (I[i-1] - I[i]), ain being the day's
# inflow per head (newborns arrive with none); the dead carry the cell mean, so
# death leaves it alone. Do not "fix" this to r_age[i-1] * I[i-1] by analogy
# with the disease compartments.
ain[, ] <- if (i == 1) births * het_wt[j] / (Npop[1, j] + pop_floor) else
  r_age[i - 1] * Npop[i - 1, j] / (Npop[i, j] + pop_floor)
dim(pbB, pbC, pbD, pbV, ain) <- c(n_age, n_het)
# A boosted person gains 1 and skips the day's decay (their update replaces it);
# everyone else decays by exp(-1/d). Averaged over the cell:
#   I' = exp(-1/d) I + p_boost (1 + I (1 - exp(-1/d))).
dec_ib <- 1 - exp(-1 / d_ib); dec_ica <- 1 - exp(-1 / d_ica)
dec_id <- 1 - exp(-1 / d_id); dec_iva <- 1 - exp(-1 / d_iva)
update(IB[, ]) <- IB[i, j] + pbB[i, j] - (1 - pbB[i, j]) * dec_ib * IB[i, j] +
  ain[i, j] * ((if (i > 1) IB[i - 1, j] else 0) - IB[i, j])
update(ICA[, ]) <- ICA[i, j] + pbC[i, j] - (1 - pbC[i, j]) * dec_ica * ICA[i, j] +
  ain[i, j] * ((if (i > 1) ICA[i - 1, j] else 0) - ICA[i, j])
update(ID[, ]) <- ID[i, j] + pbD[i, j] - (1 - pbD[i, j]) * dec_id * ID[i, j] +
  ain[i, j] * ((if (i > 1) ID[i - 1, j] else 0) - ID[i, j])
update(IVA[, ]) <- IVA[i, j] + pbV[i, j] - (1 - pbV[i, j]) * dec_iva * IVA[i, j] +
  ain[i, j] * ((if (i > 1) IVA[i - 1, j] else 0) - IVA[i, j])
dim(IB, ICA, ID, IVA) <- c(n_age, n_het)

## ---- onward infectivity and the force of infection on mosquitoes -----------------------
# TBV reduces onward infectivity, per infection state
inf[, ] <- cD * D[i, j] * tbv_fD[i] + cA[i, j] * A[i, j] * tbv_fA[i] +
  cU * U[i, j] * tbv_fU[i] + (cT * (Tr[i, j] + Tr_slow[i, j]) + J_c[i, j] + J_cs[i, j]) * tbv_fT[i]
infw[, ] <- zeta[j] * psi[i] * inf[i, j]
dim(inf, infw) <- c(n_age, n_het)
# malariasimulation weights each person by zeta x psi, normalised over the live
# population (biting_process.R, human_pi), so the infectivity reaching mosquitoes
# is sum(zeta psi inf) / sum(zeta psi). Taken from the current population, so it
# follows a shifting age structure.
psi_n[] <- psi[i] * n_g_now[i]
zp_n[, ] <- zeta[j] * psi[i] * Npop[i, j]
dim(psi_n) <- n_age
dim(zp_n) <- c(n_age, n_het)
mpsi <- sum(psi_n) / sum(n_g_now)
mzp <- sum(zp_n) / sum(n_g_now)
inf_now <- sum(infw) / mzp
# lagged_infectivity$get(t - delay_gam)
inf_used <- (1 - fl_frac) * (if (fl_floor == 0) inf_now else inf_hist[fl_floor]) +
  fl_frac * inf_hist[fl_floor + 1]
aIm[] <- a_spp[i] * Im[i]
dim(aIm) <- n_spp
foim[] <- a_spp[i] * inf_used
dim(foim) <- n_spp

## ---- mosquitoes: one day of malariasimulation's ODE, per species -------------------------
del <- parameter(); dl <- parameter(); dpl <- parameter()
me <- parameter(); ml <- parameter(); mup <- parameter(); mosq_gamma <- parameter()
# oviposition: eggs_laid(beta, mu, f) is algebraically beta for every mu and f,
# so it is constant (vector control acts through mu, not fecundity)
beta_eff <- parameter(); dim(beta_eff) <- n_spp
dem <- parameter()
# Larval carrying capacity. The IBM's aquatic model evaluates it at the solver's
# time truncated to a whole day (carrying_capacity() takes a size_t timestep), so
# it is held at the start-of-day value across the day.
n_cct <- parameter(constant = TRUE)
cc_times <- parameter(); dim(cc_times) <- n_cct
K_vals <- parameter(); dim(K_vals) <- c(n_spp, n_cct)
Kcap <- interpolate(cc_times, K_vals, "linear"); dim(Kcap) <- n_spp
# a (human blood-meal rate) and mu (adult death) under nets and IRS, per species
n_vct <- parameter(constant = TRUE)
vc_times <- parameter(); dim(vc_times) <- n_vct
a_vals <- parameter(); mum_vals <- parameter()
dim(a_vals) <- c(n_spp, n_vct)
dim(mum_vals) <- c(n_spp, n_vct)
a_spp <- interpolate(vc_times, a_vals, "linear"); dim(a_spp) <- n_spp
mum <- interpolate(vc_times, mum_vals, "linear"); dim(mum) <- n_spp
# Incubation. adult_mosquito_model_update pushes today's S * foim onto a queue of
# ceiling(dem) days and pops the oldest, so the flux completing incubation today
# is the one pushed ceiling(dem) - 1 days ago, surviving at exp(-mu * dem) with
# today's mu (adult_mosquito_eqs.cpp).
lag_in[] <- if (n_tau == 0) Sm[i] * foim[i] else inc_hist[i, n_tau]
lagsurv[] <- lag_in[i] * exp(-mum[i] * dem)
dim(lag_in, lagsurv) <- n_spp

# The day is integrated in n_sub sub-steps of an exponential midpoint rule. Each
# compartment obeys dX/dt = I(state) - R(state) X, and a step of length h takes
# it to X e^{-R h} + (I / R)(1 - e^{-R h}): exact for frozen I and R, stable and
# positive however large R gets. That matters because the late-larval
# density-dependent mortality is stiff (about 5 a day at equilibrium), and the
# larvae can start a run thousands of times over their carrying capacity, the
# seed being the annual mean and K seasonal. Each sub-step predicts the
# half-step state and then takes the full step with I and R taken from it:
# second order.
#
# Every stage is a slot of one array, so a stage can read the whole state
# before it: slot 1 is the day's start, then per sub-step s the half-step state
# (slot 2s) and the next state (slot 2s + 1). An even slot steps half a
# sub-step from the slot before it, with rates from that slot; an odd slot takes
# a full sub-step from two slots back, with rates from the half-step slot.
# Components: 1 early larvae, 2 late larvae, 3 pupae, 4 S, 5 E, 6 I.
nq <- 2 * n_sub + 1
hs <- 1 / n_sub
sbase[] <- if (i == 1) 1 else if (i %% 2 == 0) i - 1 else i - 2
shh[] <- if (i %% 2 == 0) hs / 2 else hs
dim(sbase, shh) <- nq
RP <- 1 / dpl + mup
RS[] <- foim[i] + mum[i]
dim(RS) <- n_spp
# The larval rates depend on the state they are taken from, so they are written
# out in place: the array may refer only to itself.
Z[1, , ] <- if (j == 1) ME[k] else if (j == 2) ML[k] else if (j == 3) MP[k] else
  if (j == 4) Sm[k] else if (j == 5) Em[k] else Im[k]
Z[2:nq, , ] <- (if (j == 1)
  Z[sbase[i], 1, k] * exp(-(1 / del + me * (1 + (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) * shh[i]) +
  beta_eff[k] * (Z[i - 1, 4, k] + Z[i - 1, 5, k] + Z[i - 1, 6, k]) /
    (1 / del + me * (1 + (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) *
    (1 - exp(-(1 / del + me * (1 + (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) * shh[i]))
  else if (j == 2)
  Z[sbase[i], 2, k] * exp(-(1 / dl + ml * (1 + mosq_gamma * (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) * shh[i]) +
  Z[i - 1, 1, k] / del /
    (1 / dl + ml * (1 + mosq_gamma * (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) *
    (1 - exp(-(1 / dl + ml * (1 + mosq_gamma * (Z[i - 1, 1, k] + Z[i - 1, 2, k]) / Kcap[k])) * shh[i]))
  else if (j == 3)
  Z[sbase[i], 3, k] * exp(-RP * shh[i]) +
  Z[i - 1, 2, k] / dl / RP * (1 - exp(-RP * shh[i]))
  else if (j == 4)
  Z[sbase[i], 4, k] * exp(-RS[k] * shh[i]) +
  0.5 * Z[i - 1, 3, k] / dpl / RS[k] * (1 - exp(-RS[k] * shh[i]))
  else if (j == 5)
  Z[sbase[i], 5, k] * exp(-mum[k] * shh[i]) +
  (Z[i - 1, 4, k] * foim[k] - lagsurv[k]) / mum[k] * (1 - exp(-mum[k] * shh[i]))
  else
  Z[sbase[i], 6, k] * exp(-mum[k] * shh[i]) + lagsurv[k] / mum[k] * (1 - exp(-mum[k] * shh[i])))
dim(Z) <- c(nq, 6, n_spp)
# The 0.5 above is the SEX RATIO, not a loss: emerging pupae are half female,
# and only females are modelled from there on (as in malariasimulation).

update(ME[]) <- Z[nq, 1, i]
update(ML[]) <- Z[nq, 2, i]
update(MP[]) <- Z[nq, 3, i]
update(Sm[]) <- Z[nq, 4, i]
update(Em[]) <- Z[nq, 5, i]
update(Im[]) <- Z[nq, 6, i]
dim(ME, ML, MP, Sm, Em, Im) <- n_spp

## ---- the three delay lines -----------------------------------------------------------------
# column m holds the value saved m days ago
update(eir_hist[, 1]) <- aIm[i]
update(eir_hist[, 2:n_eirh]) <- eir_hist[i, j - 1]
update(inc_hist[, 1]) <- Sm[i] * foim[i]
update(inc_hist[, 2:n_inch]) <- inc_hist[i, j - 1]
update(inf_hist[1]) <- inf_now
update(inf_hist[2:n_infh]) <- inf_hist[i - 1]
dim(eir_hist) <- c(n_spp, n_eirh)
dim(inc_hist) <- c(n_spp, n_inch)
dim(inf_hist) <- n_infh

## ---- outputs: the day's counts, recorded at the end of the day ---------------------------
# State counts and prevalence are of the day's start, as the IBM renders them;
# incidence is the day's. Severe disease is drawn from the same infected set as
# clinical in the IBM (update_severe_disease takes infected_humans), so it is
# theta x all infections. PCR follows the IBM convention: every D/Tr/A/U counts.
detlm[, ] <- D[i, j] + Tr[i, j] + Tr_slow[i, j] + Trc_tot[i, j] + q[i, j] * A[i, j]
detpcr[, ] <- D[i, j] + Tr[i, j] + Tr_slow[i, j] + Trc_tot[i, j] + A[i, j] + U[i, j]
sev_n[, ] <- theta[i, j] * inf_tot[i, j]
dim(detlm, detpcr, sev_n) <- c(n_age, n_het)
update(n_g[]) <- n_g_now[i]
update(det_lm_g[]) <- sum(detlm[i, ])
update(det_pcr_g[]) <- sum(detpcr[i, ])
update(clin_g[]) <- sum(clin_n[i, ])
update(sev_g[]) <- sum(sev_n[i, ])
update(inc_g[]) <- sum(inf_tot[i, ])
update(S_g[]) <- sum(S[i, ])
update(D_g[]) <- sum(D[i, ])
update(A_g[]) <- sum(A[i, ])
update(U_g[]) <- sum(U[i, ])
update(Tr_g[]) <- sum(Tr[i, ]) + sum(Tr_slow[i, ]) + sum(Trc_tot[i, ])
update(Ph_g[]) <- sum(Ph_tot[i, ]) + sum(Phc_tot[i, ])
update(EIR_yr) <- eir_lag * 365            # the EIR biting humans today, as EIR_<species>
update(FOIM) <- foim[1]
update(ft_out) <- ft
dim(n_g, det_lm_g, det_pcr_g, clin_g, sev_g, inc_g) <- n_age
dim(S_g, D_g, A_g, U_g, Tr_g, Ph_g) <- n_age

## ---- initial conditions ---------------------------------------------------
S0 <- parameter(); D0 <- parameter(); A0 <- parameter()
U0 <- parameter(); Tr0 <- parameter(); Ph0 <- parameter(); Phc0 <- parameter()
IB_init <- parameter(); ICA_init <- parameter(); ID_init <- parameter(); IVA_init <- parameter()
dim(S0, D0, A0, U0, Tr0, IB_init, ICA_init, ID_init, IVA_init) <- c(n_age, n_het)
dim(Ph0) <- c(n_age, n_het, n_ph)
dim(Phc0) <- c(n_age, n_het, n_phc)
ME0 <- parameter(); ML0 <- parameter(); MP0 <- parameter(); Sm0 <- parameter()
Em0 <- parameter(); Im0 <- parameter()
dim(ME0, ML0, MP0, Sm0, Em0, Im0) <- n_spp
# the delay lines start full of the seed's values, as the IBM's LaggedValue
# defaults and its incubation queue do
eir0 <- parameter(); inc0 <- parameter(); inf0 <- parameter()
dim(eir0, inc0) <- n_spp

initial(S[, ]) <- S0[i, j]
initial(D[, ]) <- D0[i, j]
initial(A[, ]) <- A0[i, j]
initial(U[, ]) <- U0[i, j]
initial(Tr[, ]) <- (1 - spc0) * Tr0[i, j]
initial(Tr_slow[, ]) <- spc0 * Tr0[i, j]
initial(Ph[, , ]) <- Ph0[i, j, k]
initial(Ph_c[, , ]) <- Phc0[i, j, k]
initial(Tr_c[, ]) <- 0
initial(Tr_cs[, ]) <- 0
initial(J_c[, ]) <- 0
initial(J_cs[, ]) <- 0
initial(Ph_ct[, , ]) <- 0
initial(IB[, ]) <- IB_init[i, j]
initial(ICA[, ]) <- ICA_init[i, j]
initial(ID[, ]) <- ID_init[i, j]
initial(IVA[, ]) <- IVA_init[i, j]
initial(ME[]) <- ME0[i]
initial(ML[]) <- ML0[i]
initial(MP[]) <- MP0[i]
initial(Sm[]) <- Sm0[i]
initial(Em[]) <- Em0[i]
initial(Im[]) <- Im0[i]
initial(eir_hist[, ]) <- eir0[i]
initial(inc_hist[, ]) <- inc0[i]
initial(inf_hist[]) <- inf0
initial(n_g[]) <- 0
initial(det_lm_g[]) <- 0
initial(det_pcr_g[]) <- 0
initial(clin_g[]) <- 0
initial(sev_g[]) <- 0
initial(inc_g[]) <- 0
initial(S_g[]) <- 0
initial(D_g[]) <- 0
initial(A_g[]) <- 0
initial(U_g[]) <- 0
initial(Tr_g[]) <- 0
initial(Ph_g[]) <- 0
initial(EIR_yr) <- 0
initial(FOIM) <- 0
initial(ft_out) <- 0
