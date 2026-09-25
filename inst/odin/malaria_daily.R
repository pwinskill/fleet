# =============================================================================
# fleet: mean-field twin of malariasimulation, on its daily clock
# =============================================================================
# P. falciparum: human S/D/A/U/Tr over [age, het] with the six Griffin
# immunity functions. P. vivax: the same disease states over
# [age, het, hypnozoite level], with IAA/ICA immunity and its within-cell spread
# (the "P. vivax human block" below). One parasite per run, as in
# malariasimulation; the other block is held empty. Both are coupled to one
# mosquito model.
#
# The population is advanced one day at a time, in the order malariasimulation
# resolves a day. Every process reads the state at the start of the day and
# every change lands at the end of it, as the IBM queues its updates:
#
#   1. Immunity decays by exp(-1/d). A person boosted today skips the decay:
#      the IBM queues the decay first and the boost (start-of-day value + 1)
#      after it, and the later update overwrites the earlier.
#   2. Biting. The EIR saved `de` days ago reaches each age x heterogeneity
#      stratum. Falciparum deduplicates the day's bites (a person bitten several
#      times counts once), boosts IB for everyone bitten, and infects a bitten
#      person in S, A or U with probability b x PEV. Vivax counts every bite,
#      adds the relapse hazard of each person's hypnozoite batches, and reduces
#      the total by PEV.
#   3. Infection and each state's progression are ONE competing draw a day
#      (CompetingHazard$resolve): an infection pre-empts that day's progression.
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

# Keeps a division by an empty stratum finite. An empty block divides by its own
# zero population, and 0/0 anywhere would reach a shared quantity as NaN (0 x NaN
# is NaN), so divisions by a cell population carry pop_floor in the denominator:
# an empty cell gives 0/pop_floor = 0, and in an occupied one it is an exact
# no-op. Not `if (N > 0) ... else 0`: a branch in a loop that also computes exp or
# log stops the compiler optimising it, and made a run 10x slower at the same
# number of steps. Written, like every constant here that is not a short binary
# fraction, as a quotient of whole numbers: odin2 prints a fractional literal at
# 17 digits, which differs between platforms, and a whole number the same on all.
pop_floor <- 1 / 1e300

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

## ---- parasite ---------------------------------------------------------------
# Both blocks are always compiled; the one not in use is held empty and inert.
# The vivax block has dimensions of its own, so under falciparum it collapses to
# one cell -- n_age_v = n_het_v = n_hyp = 1 -- and costs a handful of states. The
# falciparum block is not collapsed under vivax: it shares n_age and n_het with
# the demography, biting and vaccine arrays the vivax block needs at full size,
# and sits there empty, a small fraction of a vivax day.
pf_on <- parameter(constant = TRUE)    # 1 under falciparum, 0 under vivax
pv_on <- parameter(constant = TRUE)    # 1 under vivax, 0 under falciparum
n_age_v <- parameter(constant = TRUE)  # n_age under vivax, 1 otherwise
n_het_v <- parameter(constant = TRUE)  # n_het under vivax, 1 otherwise
n_bat <- parameter(constant = TRUE)    # batch levels: kmax + 1 under vivax, 1 otherwise
n_hyp <- parameter(constant = TRUE)    # all levels: n_bat plus liver-stage protection
n_phv <- parameter(constant = TRUE)    # vivax post-treatment prophylaxis stages
n_q <- parameter(constant = TRUE)      # quadrature nodes for vivax immunity spread

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
# rA, rD, rU, rT are the IBM's whole-day exit probabilities, 1 - exp(-1/d). A
# vivax list carries da, dd and dt under the falciparum names, so under vivax
# these are the vivax ones.
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
# PER-INDIVIDUAL detail. The falciparum block reads its curves at the stratum
# mean, where it does NOT translate to +0.5, and an A/B against the IBM ensemble
# mean favours 0, its default. The vivax block reads them at quadrature nodes that
# stand for people, where it does, and takes 0.5 (build_inputs()).
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
# first-line drug switch changes the mix, not just total coverage; and under
# vivax the radically cured share of treatment (hyp_mix), the share of it whose
# blood stage clears too (eff_hyp) and the mean liver-stage protection (mean_ls)
n_dmix <- parameter(constant = TRUE)
dmix_times <- parameter(); dim(dmix_times) <- n_dmix
cT_vals <- parameter(); drug_eff_vals <- parameter(); rP_vals <- parameter()
hyp_vals <- parameter(); eff_hyp_vals <- parameter(); mean_ls_vals <- parameter()
dim(cT_vals, drug_eff_vals, rP_vals, hyp_vals, eff_hyp_vals, mean_ls_vals) <- n_dmix
cT <- interpolate(dmix_times, cT_vals, "constant")
drug_eff <- interpolate(dmix_times, drug_eff_vals, "constant")
rP <- interpolate(dmix_times, rP_vals, "constant")
hyp_mix <- interpolate(dmix_times, hyp_vals, "constant")
eff_hyp <- interpolate(dmix_times, eff_hyp_vals, "constant")
mean_ls <- interpolate(dmix_times, mean_ls_vals, "constant")
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
## A newborn's mother is drawn from those aged [20, 21) in its heterogeneity
## group, so it inherits their population-weighted mean; mask20 is each age
## group's share of that year.
ICAmask[, ] <- ICA[i, j] * Npop[i, j] * mask20[i]
IVAmask[, ] <- IVA[i, j] * Npop[i, j] * mask20[i]
Nmask20[, ] <- Npop[i, j] * mask20[i]
dim(ICAmask, IVAmask, Nmask20) <- c(n_age, n_het)
ICA20[] <- sum(ICAmask[, i]) / (sum(Nmask20[, i]) + pop_floor)
IVA20[] <- sum(IVAmask[, i]) / (sum(Nmask20[, i]) + pop_floor)
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
# sum(w zeta), 0.99972 at 5 nodes, not 1. The same draw serves both parasites.
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
# Both blocks' deaths are replaced, and the births go to whichever block is in
# use: under falciparum the vivax block is empty and adds exactly 0.
births <- sum(deaths) + sum(deaths_v)
births_f <- births * pf_on
births_v <- births * pv_on
# the living population by age and stratum, both blocks (one of them is empty)
n_g_now[] <- sum(Npop[i, ]) + nv_g[i]
dim(n_g_now) <- n_age

## ---- the day's transitions -----------------------------------------------------------
# Each compartment keeps what does not leave it, gains what ages in from the
# group below (births at the bottom) and loses what ages out or dies. A clinical
# case goes to Tr if treated successfully and to D otherwise; a non-clinical
# infection goes to A (an A stays A).
update(S[, ]) <- S[i, j] + (if (i > 1) r_age[i - 1] * S[i - 1, j] else births_f * het_wt[j]) -
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
ain[, ] <- if (i == 1) births_f * het_wt[j] / (Npop[1, j] + pop_floor) else
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

## =====================================================================================
## ---- P. vivax human block --------------------------------------------------------------
## =====================================================================================
# Indexed [age, heterogeneity, hypnozoite level]. malariasimulation carries an
# integer batch count per person, k in 0..kmax, and relapse is a hazard of k*f.
# A mean-field model cannot hold only the mean count: the relapse hazard is
# linear in k, but immunity, detectability and the U->S rate are not, and
# malariaEquilibriumVivax returns every state as [age, het, batch] for the same
# reason. kk[k] is the batch count at level k, so kk = 0 at the first level.
#
# Levels 1..n_bat are the batch counts 0..kmax. Levels after n_bat hold people
# with no batches whose liver stage is still drug-protected after radical cure:
# a chain of stages from the dose, left for level 1. Liver-stage protection is a
# property of the person, not of their disease state, so it has to sit on the
# dimension every compartment and immunity stock already carries. On those
# levels kk = 0 (no relapse, nothing to clear) and a bite infects but forms no
# batch, which is what ls_prophylaxis does in the IBM. Without a radical-cure
# drug there are no such levels and n_hyp = n_bat.
#
# Boundaries are handled with `if` guards on the index, not with index ranges:
# at n_hyp = 1 -- the falciparum collapse -- the batch-1 and batch-n_hyp
# boundary equations of a range would both target the same element, and the
# second would read k - 1 = 0. odin2 warns that it cannot validate accesses such
# as X[i, j, k + 1]; they are safe, because the generated code is a ternary and
# the out-of-range branch is never evaluated. Keep guards out of the arrays that
# compute exp/log -- see pop_floor for what a branch there costs.
gammal <- parameter()     # per-batch hypnozoite clearance (k -> k-1 at k*gammal)
ff <- parameter()         # per-batch relapse rate (hazard k*ff)
bv <- parameter()         # infection probability per infectious bite; no IB for vivax
philm_min <- parameter(); philm_max <- parameter(); alm50 <- parameter(); klm <- parameter()
dpcr_min <- parameter(); dpcr_max <- parameter(); apcr50 <- parameter(); kpcr <- parameter()
cA_v <- parameter()       # infectivity of A, constant for vivax (ca)
d_iaa <- parameter()      # anti-parasite immunity decay (ra)
ua_eff <- parameter()     # its refractory window, ceiling(ua) - 1
# The refractory windows as chains of stages (see the refractory stocks below):
# n stages over a window of u days, each passed on with probability n / u a day.
n_ra <- parameter(constant = TRUE)    # IAA's window (ua)
n_rc <- parameter(constant = TRUE)    # ICA's window (uc)
# within-cell immunity spread: on/off, and the Gauss-Hermite rule it is read at
spread_on <- parameter()
qz <- parameter(); qw <- parameter()
dim(qz, qw) <- n_q

kk[] <- if (i <= n_bat) as.numeric(i) - 1 else 0
shmask[] <- if (i < n_bat) 1 else 0       # a bite here moves the person up a batch
lsmask[] <- if (i > n_bat) 1 else 0       # liver-stage-protected levels
dim(kk, shmask, lsmask) <- n_hyp
k_rc <- if (n_hyp > n_bat) n_bat + 1 else 1   # where radical cure puts people
rls_k <- (n_hyp - n_bat) / mean_ls            # per-stage exit of the liver-stage chain
rPk_v <- rP * n_phv                           # rP is the vivax drug's under vivax
# hypnozoite_batch_decay_process: a carrier of k batches loses one today with
# probability rate_to_prob(k * gammal); create_exponential_decay_process
# multiplies immunity by exp(-1/d) a day, so its square by exp(-2/d)
pbl[] <- 1 - exp(-gammal * kk[i])
dim(pbl) <- n_hyp
e1_iaa <- exp(-1 / d_iaa); e2_iaa <- exp(-2 / d_iaa)
e1_ica <- exp(-1 / d_ica); e2_ica <- exp(-2 / d_ica)

dim(Sv, Dv, Av, Uv, Trv, Trv_slow, JAv, JCv, KAv, KCv) <- c(n_age_v, n_het_v, n_hyp)
dim(Phv) <- c(n_age_v, n_het_v, n_hyp, n_phv)

## Population. Treated and drug-protected people are not exposed to infection:
## at treatment the IBM's prophylaxis is 1, which fleet represents by holding them
## in Tr and Phv. Everyone else is, including D -- an infection does not change a
## D person's disease state, but it does boost them and, if it came from a bite,
## give them a new batch.
Phv_tot[, , ] <- sum(Phv[i, j, k, ])
Npop_v[, , ] <- Sv[i, j, k] + Dv[i, j, k] + Av[i, j, k] + Uv[i, j, k] + Trv[i, j, k] +
  Trv_slow[i, j, k] + Phv_tot[i, j, k]
deaths_v[, , ] <- mu_age[i] * Npop_v[i, j, k]
dim(Phv_tot, Npop_v, deaths_v) <- c(n_age_v, n_het_v, n_hyp)
Nv_ij[, ] <- sum(Npop_v[i, j, ])
dim(Nv_ij) <- c(n_age_v, n_het_v)
Nv_age[] <- sum(Npop_v[i, , ])
dim(Nv_age) <- n_age_v
# The vivax population on the shared age and heterogeneity indices, so it joins
# the population the day's bites are shared over. Under falciparum n_age_v =
# n_het_v = 1, the one vivax cell is empty, and these are exactly 0.
nv_g[] <- if (i <= n_age_v) Nv_age[i] else 0
dim(nv_g) <- n_age
nv_ij[, ] <- if (i <= n_age_v && j <= n_het_v) Nv_ij[i, j] else 0
dim(nv_ij) <- c(n_age, n_het)

## Immunity is held as a STOCK, J = N * I -- the total in a cell -- not the mean
## that the falciparum block holds. The two are the same mathematics. The
## difference is numerical: in mean form every inflow is divided by the
## receiving cell's population, and the vivax ladder has cells that are exactly
## empty (a newborn cannot carry ten batches), so an empty cell with a populated
## neighbour would take an inflow of order 1e291. In stock form nothing divides
## by a population except turning a stock back into the mean the Hill functions
## need, and an empty cell there simply reads 0.

## Maternal immunity (algebraic). Inherited at birth -- at batch 0, since
## newborns carry no hypnozoites -- from a mother of the same heterogeneity
## group aged 20 to 21, and decaying with age regardless of the batches acquired
## later, so it is [age, het] and applies at every batch. The mother's value is
## the population-weighted mean over the year and over batches, which the stocks
## give directly. malariasimulation uses the same pcm and waning rate rm for
## both, which a vivax list carries as PM and the decay behind icm_factor.
ICAv_m20[, ] <- sum(JCv[i, j, ]) * mask20[i]
IAAv_m20[, ] <- sum(JAv[i, j, ]) * mask20[i]
Nv_m20[, ] <- Nv_ij[i, j] * mask20[i]
dim(ICAv_m20, IAAv_m20, Nv_m20) <- c(n_age_v, n_het_v)
ICAv20[] <- sum(ICAv_m20[, i]) / (sum(Nv_m20[, i]) + pop_floor)
IAAv20[] <- sum(IAAv_m20[, i]) / (sum(Nv_m20[, i]) + pop_floor)
dim(ICAv20, IAAv20) <- n_het_v
ICMv[, ] <- PM * ICAv20[j] * icm_factor[i]
IAMv[, ] <- PM * IAAv20[j] * icm_factor[i]
dim(ICMv, IAMv) <- c(n_age_v, n_het_v)

## Immunity -> probability. anti_parasite_immunity() in malariasimulation adds
## NO +0.5 to acquired immunity; clinical_immunity() does, to each person's, and
## the nodes below are people, so acq_offset is 0.5 here by default (without it
## school-age clinical incidence ran 2 points high at EIR 10). A vivax list
## carries its clinical curve (phi0, phi1, ic0, kc) under the falciparum names.
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
## Means are clamped at zero: the ladder holds cells as small as 1e-36 people,
## rounding can leave such a cell fractionally negative, and a fractional power
## of a negative mean is NaN. max() compiles to one branchless instruction.
mA[, , ] <- max(JAv[i, j, k] / (Npop_v[i, j, k] + pop_floor), 0)
mC[, , ] <- max(JCv[i, j, k] / (Npop_v[i, j, k] + pop_floor), 0)
cv2A[, , ] <- spread_on * min(max(KAv[i, j, k] / (Npop_v[i, j, k] + pop_floor) -
  mA[i, j, k] * mA[i, j, k], 0) / (mA[i, j, k] * mA[i, j, k] + 1 / 1e12), 9)
cv2C[, , ] <- spread_on * min(max(KCv[i, j, k] / (Npop_v[i, j, k] + pop_floor) -
  mC[i, j, k] * mC[i, j, k], 0) / (mC[i, j, k] * mC[i, j, k] + 1 / 1e12), 9)
whaA[, , ] <- 1 - cv2A[i, j, k] / 9
whaC[, , ] <- 1 - cv2C[i, j, k] / 9
whbA[, , ] <- sqrt(cv2A[i, j, k]) / 3
whbC[, , ] <- sqrt(cv2C[i, j, k]) / 3
dim(mA, mC, cv2A, cv2C, whaA, whaC, whbA, whbC) <- c(n_age_v, n_het_v, n_hyp)

## Infection. Vivax counts every bite, 1 - (1 - b)^n, rather than deduplicating
## them; for Poisson bites at rate EPS that averages to 1 - exp(-b EPS). Relapse
## adds k ff. Both are summed and THEN reduced by PEV, as in
## calculate_vivax_infections(), so the vaccine blocks relapses too; iS_v is the
## day's infection probability, and hv = -log(1 - iS_v) the infection rate that
## function hands to the competing hazards. Only an infection that came from a
## bite forms a new batch, which is what the relative_rates draw in
## relapse_bite_infection_hazard_resolution() decides; sh_v is that bite share
## wherever a batch can form -- below the top batch, which stays at kmax
## (pmin(k + 1, kmax)), and outside liver-stage protection.
lam_bv[, ] <- bv * EPS[i, j]
dim(lam_bv) <- c(n_age_v, n_het_v)
r_tot_v[, , ] <- lam_bv[i, j] + kk[k] * ff
iS_v[, , ] <- (1 - exp(-r_tot_v[i, j, k])) * pev_factor[i]
hv[, , ] <- -log(1 - iS_v[i, j, k])
sh_v[, , ] <- shmask[k] * lam_bv[i, j] / (r_tot_v[i, j, k] + pop_floor)
dim(r_tot_v, iS_v, hv, sh_v) <- c(n_age_v, n_het_v, n_hyp)

## The day, as malariasimulation resolves it. Each person faces infection (hv)
## and their own state's progression (U -> S, A -> U, D -> A) in ONE competing
## hazard draw a day (CompetingHazard$resolve): some event happens with
## probability 1 - exp(-(hv + h)), and it is infection with share hv / (hv + h).
## So infection and recovery each make the other less likely that day, and a
## person is infected at most once. S has no progression, so its infection
## probability is iS_v; A and D progress at the whole-day probabilities rA and
## rD, whose hazards are hA and hD, so exp(-(hv + hA)) is (1 - iS_v) (1 - rA)
## and needs no exp of its own. U's progression rate depends on IAA, so U's two
## probabilities are taken per node.
# The three Hill curves as exp(k (log x - log x50)) rather than (x / x50)^k:
# the same function, but the two curves of IAA share one log, and a pow() is a
# log and an exp, so this saves a log per node. x is never negative (the means
# are clamped), and log(0) = -Inf gives exp() = 0, the value pow(0, k) has. x is
# the Wilson-Hilferty node, m u^3 plus maternal immunity with u = max(a + b z,
# 0), written out in place rather than held in an array of its own.
lxA[, , , ] <- log(mA[i, j, k] * max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0) *
  max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0) * max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0) +
  IAMv[i, j])
lxC[, , , ] <- log(mC[i, j, k] * max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0) *
  max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0) * max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0) +
  acq_offset + ICMv[i, j])
pLMn[, , , ] <- philm_min + (philm_max - philm_min) / (1 + exp(klm * (lxA[i, j, k, l] - l_alm50)))
pDn[, , , ] <- phi0 * (phi1 + (1 - phi1) / (1 + exp(kc * (lxC[i, j, k, l] - l_ic0))))
rUn[, , , ] <- 1 / (dpcr_min + (dpcr_max - dpcr_min) / (1 + exp(kpcr * (lxA[i, j, k, l] - l_apcr50))))
eUnv[, , , ] <- 1 - exp(-(hv[i, j, k] + rUn[i, j, k, l]))
## weighted node terms, summed below: LM-detectable, LM-detectable AND
## clinical, clinical; and for U, infected, infected and LM-detectable,
## infected and LM-detectable and clinical, and any event at all
w_lm[, , , ] <- qw[l] * pLMn[i, j, k, l]
w_lmc[, , , ] <- w_lm[i, j, k, l] * pDn[i, j, k, l]
w_c[, , , ] <- qw[l] * pDn[i, j, k, l]
w_iu[, , , ] <- qw[l] * eUnv[i, j, k, l] * hv[i, j, k] / (hv[i, j, k] + rUn[i, j, k, l])
w_iulm[, , , ] <- w_iu[i, j, k, l] * pLMn[i, j, k, l]
w_iulmc[, , , ] <- w_iulm[i, j, k, l] * pDn[i, j, k, l]
w_eu[, , , ] <- qw[l] * eUnv[i, j, k, l]
dim(lxA, lxC, pLMn, pDn, rUn, eUnv, w_lm, w_lmc, w_c, w_iu, w_iulm, w_iulmc, w_eu) <-
  c(n_age_v, n_het_v, n_hyp, n_q)
l_alm50 <- log(alm50); l_ic0 <- log(ic0); l_apcr50 <- log(apcr50)
hD <- -log(1 - rD)
s_lm[, , ] <- sum(w_lm[i, j, k, ])
s_lmc[, , ] <- sum(w_lmc[i, j, k, ])
s_c[, , ] <- sum(w_c[i, j, k, ])
iU_v[, , ] <- sum(w_iu[i, j, k, ])
s_iulm[, , ] <- sum(w_iulm[i, j, k, ])
s_iulmc[, , ] <- sum(w_iulmc[i, j, k, ])
recU_v[, , ] <- sum(w_eu[i, j, k, ]) - iU_v[i, j, k]
eA_v[, , ] <- 1 - (1 - iS_v[i, j, k]) * (1 - rA)
iA_v[, , ] <- eA_v[i, j, k] * hv[i, j, k] / (hv[i, j, k] + hA)
recA_v[, , ] <- eA_v[i, j, k] - iA_v[i, j, k]
eD_v[, , ] <- 1 - (1 - iS_v[i, j, k]) * (1 - rD)
iD_v[, , ] <- eD_v[i, j, k] * hv[i, j, k] / (hv[i, j, k] + hD)
recD_v[, , ] <- eD_v[i, j, k] - iD_v[i, j, k]
dim(s_lm, s_lmc, s_c, iU_v, s_iulm, s_iulmc, recU_v, eA_v, iA_v, recA_v,
    eD_v, iD_v, recD_v) <- c(n_age_v, n_het_v, n_hyp)

## Where infections go. vivax_infection_outcome_process(): S and U pass through
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

## The day's batch-formers keep their batches. A bite that forms a batch queues
## pmin(k + 1, kmax) from the start-of-day k after
## hypnozoite_batch_decay_process() has queued k - 1, and the later update
## overwrites the earlier -- at the top batch too, where the bite writes kmax
## back without moving anyone up. Everyone else's batch decay is a draw of its
## own, independent of the day's other events, so it is taken after them, from
## the day's outcome less its batch-formers (below). ex_* are the batch-formers
## among each flux's arrivals: bitten one level down, and at the top batch
## bitten in place.
shu_v[, , ] <- if (k == n_bat) lam_bv[i, j] / (r_tot_v[i, j, k] + pop_floor) else 0
ex_U[, , ] <- (if (k > 1) sh_v[i, j, k - 1] * toU_n[i, j, k - 1] else 0) +
  shu_v[i, j, k] * toU_n[i, j, k]
ex_A[, , ] <- (if (k > 1) sh_v[i, j, k - 1] * toA_n[i, j, k - 1] else 0) +
  shu_v[i, j, k] * toA_n[i, j, k]
ex_C[, , ] <- (if (k > 1) sh_v[i, j, k - 1] * toC_n[i, j, k - 1] else 0) +
  shu_v[i, j, k] * toC_n[i, j, k]
ex_D[, , ] <- (if (k > 1) sh_v[i, j, k - 1] * toD_n[i, j, k - 1] else 0) +
  shu_v[i, j, k] * toD_n[i, j, k]
fex_C[, , ] <- ex_C[i, j, k] / (arr_C[i, j, k] + pop_floor)
dim(shu_v, ex_U, ex_A, ex_C, ex_D, fex_C) <- c(n_age_v, n_het_v, n_hyp)

## Treatment splits the clinical flux exactly as for falciparum: ft_eff to
## treatment (split fast/slow clearance by spc), the rest -- untreated and early
## treatment failures -- to D.
trt_v[, , ] <- ft_eff * arr_C[i, j, k]
dim(trt_v) <- c(n_age_v, n_het_v, n_hyp)
## Radical cure. Of everyone who is treated, drug_hypnozoite_efficacy lose every
## batch -- whether or not the blood stage cleared, since
## calculate_successful_treatments() draws the two independently -- and start
## liver-stage protection, so they move to level k_rc. rcT_v are cleared in the
## blood too (to Tr), rcD_v are not (to D). The treated flux is summed over
## levels and delivered to k_rc. rcT_v carries ft_eff's (1 - etf): an early
## treatment failure leaves the blood stage uncleared, but its liver-stage draw
## still stands, so it falls to rcD_v.
##
## Without a radical-cure drug these transfers, and the immunity they carry
## below, are all exactly zero, and evaluating them anyway costs 6-10% of a
## vivax day. So they have dimensions of their own: full size when radical cure
## is given at any point in the run (rc_on = 1), a single cell otherwise, and
## their terms in the updates are guarded by rc_on -- safe for the same reason as
## the batch-index guards, since the branch that reads a full-size index is
## never evaluated.
rc_on <- parameter(constant = TRUE)
n_age_rc <- if (rc_on == 1) n_age_v else 1
n_het_rc <- if (rc_on == 1) n_het_v else 1
n_hyp_rc <- if (rc_on == 1) n_hyp else 1
rcT_v[, , ] <- ft * eff_hyp * (1 - etf) * arr_C[i, j, k]
## the treated who lose their batches but keep the blood stage (rcD_v), and the
## treated who clear the blood stage but keep their batches (trt_keep). Each is a
## difference of two shares that can be equal -- full coverage of a fully
## effective drug -- so each is one clamped scalar rather than a difference of
## two fluxes, which fused arithmetic can leave a hair below zero.
rcD_frac <- max(hyp_mix - eff_hyp * (1 - etf), 0)
trt_keep <- max(ft_eff - ft * eff_hyp * (1 - etf), 0)
rcD_v[, , ] <- ft * rcD_frac * arr_C[i, j, k]
dim(rcT_v, rcD_v) <- c(n_age_rc, n_het_rc, n_hyp_rc)
rcT_tot[, ] <- sum(rcT_v[i, j, ])
rcD_tot[, ] <- sum(rcD_v[i, j, ])
dim(rcT_tot, rcD_tot) <- c(n_age_rc, n_het_rc)

## The day for each compartment: what stays, ageing in from the group below
## (newborns enter S at batch 0) and out or dying, and the day's disease
## transitions. Radical cure's departures are here, its arrivals at k_rc below.
qS[, , ] <- Sv[i, j, k] +
  (if (i > 1) r_age[i - 1] * Sv[i - 1, j, k] else (if (k == 1) births_v * het_wt[j] else 0)) -
  re[i] * Sv[i, j, k] +
  recU_v[i, j, k] * Uv[i, j, k] + rPk_v * Phv[i, j, k, n_phv] -
  iS_v[i, j, k] * Sv[i, j, k]
qU[, , ] <- Uv[i, j, k] +
  (if (i > 1) r_age[i - 1] * Uv[i - 1, j, k] else 0) - re[i] * Uv[i, j, k] +
  recA_v[i, j, k] * Av[i, j, k] - recU_v[i, j, k] * Uv[i, j, k] +
  arr_U[i, j, k] - iU_v[i, j, k] * Uv[i, j, k]
qA[, , ] <- Av[i, j, k] +
  (if (i > 1) r_age[i - 1] * Av[i - 1, j, k] else 0) - re[i] * Av[i, j, k] +
  recD_v[i, j, k] * Dv[i, j, k] - recA_v[i, j, k] * Av[i, j, k] +
  arr_A[i, j, k] - iA_v[i, j, k] * Av[i, j, k]
qD[, , ] <- Dv[i, j, k] +
  (if (i > 1) r_age[i - 1] * Dv[i - 1, j, k] else 0) - re[i] * Dv[i, j, k] -
  recD_v[i, j, k] * Dv[i, j, k] +
  arr_D[i, j, k] + arr_C[i, j, k] - trt_v[i, j, k] - iD_v[i, j, k] * Dv[i, j, k] -
  (if (rc_on == 1) rcD_v[i, j, k] else 0)
qT[, , ] <- Trv[i, j, k] +
  (if (i > 1) r_age[i - 1] * Trv[i - 1, j, k] else 0) - re[i] * Trv[i, j, k] -
  rT * Trv[i, j, k] +
  (1 - spc) * (if (rc_on == 1) trt_keep * arr_C[i, j, k] else trt_v[i, j, k])
qTs[, , ] <- Trv_slow[i, j, k] +
  (if (i > 1) r_age[i - 1] * Trv_slow[i - 1, j, k] else 0) - re[i] * Trv_slow[i, j, k] -
  rT_slow * Trv_slow[i, j, k] +
  spc * (if (rc_on == 1) trt_keep * arr_C[i, j, k] else trt_v[i, j, k])
## Post-treatment prophylaxis: a chain like falciparum's Ph, but carrying the
## batch dimension, because people keep losing batches while drug-protected and
## must return to the right one.
qPh[, , , ] <- Phv[i, j, k, l] +
  (if (l == 1) rT * Trv[i, j, k] + rT_slow * Trv_slow[i, j, k] else rPk_v * Phv[i, j, k, l - 1]) -
  rPk_v * Phv[i, j, k, l] +
  (if (i > 1) r_age[i - 1] * Phv[i - 1, j, k, l] else 0) - re[i] * Phv[i, j, k, l]
dim(qS, qU, qA, qD, qT, qTs) <- c(n_age_v, n_het_v, n_hyp)
dim(qPh) <- c(n_age_v, n_het_v, n_hyp, n_phv)

## Then the day's batch decay, which moves everyone down a level whatever their
## disease state -- hypnozoites die in the liver regardless of the blood -- at
## their start-of-day level's probability, but for the day's batch-formers: in
## D and in treatment, the clinical arrivals' batch-forming share (fex_C) of
## what stays. The dead were taken above and the radically cured leave for
## k_rc, which the IBM's later updates also write over the decay.
eD_x[, , ] <- ex_D[i, j, k] + fex_C[i, j, k] * (arr_C[i, j, k] - trt_v[i, j, k] -
  (if (rc_on == 1) rcD_v[i, j, k] else 0))
eT_x[, , ] <- fex_C[i, j, k] * (1 - spc) *
  (if (rc_on == 1) trt_keep * arr_C[i, j, k] else trt_v[i, j, k])
eTs_x[, , ] <- fex_C[i, j, k] * spc *
  (if (rc_on == 1) trt_keep * arr_C[i, j, k] else trt_v[i, j, k])
dim(eD_x, eT_x, eTs_x) <- c(n_age_v, n_het_v, n_hyp)
pS[, , ] <- qS[i, j, k] * (1 - pbl[k]) + (if (k < n_hyp) pbl[k + 1] * qS[i, j, k + 1] else 0)
pU[, , ] <- qU[i, j, k] - pbl[k] * (qU[i, j, k] - ex_U[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * (qU[i, j, k + 1] - ex_U[i, j, k + 1]) else 0)
pA[, , ] <- qA[i, j, k] - pbl[k] * (qA[i, j, k] - ex_A[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * (qA[i, j, k + 1] - ex_A[i, j, k + 1]) else 0)
pD[, , ] <- qD[i, j, k] - pbl[k] * (qD[i, j, k] - eD_x[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * (qD[i, j, k + 1] - eD_x[i, j, k + 1]) else 0)
pT[, , ] <- qT[i, j, k] - pbl[k] * (qT[i, j, k] - eT_x[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * (qT[i, j, k + 1] - eT_x[i, j, k + 1]) else 0)
pTs[, , ] <- qTs[i, j, k] - pbl[k] * (qTs[i, j, k] - eTs_x[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * (qTs[i, j, k + 1] - eTs_x[i, j, k + 1]) else 0)
pPh[, , , ] <- qPh[i, j, k, l] * (1 - pbl[k]) +
  (if (k < n_hyp) pbl[k + 1] * qPh[i, j, k + 1, l] else 0)
dim(pS, pU, pA, pD, pT, pTs) <- c(n_age_v, n_het_v, n_hyp)
dim(pPh) <- c(n_age_v, n_het_v, n_hyp, n_phv)
## The immunity and refractory stocks lose a batch with their people: the
## batch-formers' share of the day's outcome keeps its level, and the rest of
## a cell's stock decays at its cell's per-head rate, as the bite carries it
## up (upA_v).
qN_v[, , ] <- qS[i, j, k] + qU[i, j, k] + qA[i, j, k] + qD[i, j, k] + qT[i, j, k] +
  qTs[i, j, k] + sum(qPh[i, j, k, ])
fdec_v[, , ] <- 1 - (ex_U[i, j, k] + ex_A[i, j, k] + eD_x[i, j, k] + eT_x[i, j, k] +
  eTs_x[i, j, k]) / (qN_v[i, j, k] + pop_floor)
dim(qN_v, fdec_v) <- c(n_age_v, n_het_v, n_hyp)

## The liver-stage clock. Protection runs from the dose whatever happens to the
## blood stage, so it is common to every compartment: each protected level
## passes rls_k of its people on, and the last returns them to level 1, batch 0
## and unprotected. It acts on the day's outcome above rather than beside it,
## because its stages are short (primaquine's four span five days, so each
## passes on 0.85 a day) and a Tr or prophylaxis stage empties fast too: taken
## side by side from the start of the day the two could remove more than a cell
## holds. In this order someone who leaves Tr and moves a stage on the same day
## does both, as in the IBM. Radical cure's arrivals join level k_rc after the
## clock, so a stage is left no earlier than the day after the dose, as W(1) is
## the first day of protection.
update(Sv[, , ]) <- pS[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pS[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pS[i, j, n_hyp] else 0)
update(Uv[, , ]) <- pU[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pU[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pU[i, j, n_hyp] else 0)
update(Av[, , ]) <- pA[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pA[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pA[i, j, n_hyp] else 0)
update(Dv[, , ]) <- pD[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pD[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pD[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) rcD_tot[i, j] else 0) else 0)
update(Trv[, , ]) <- pT[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pT[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pT[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) (1 - spc) * rcT_tot[i, j] else 0) else 0)
update(Trv_slow[, , ]) <- pTs[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pTs[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pTs[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) spc * rcT_tot[i, j] else 0) else 0)
update(Phv[, , , ]) <- pPh[i, j, k, l] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pPh[i, j, k - 1, l] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pPh[i, j, n_hyp, l] else 0)

## Immunity dynamics, in stock form (see the note on stocks, above maternal
## immunity). Each flow of people carries its cell's mean with it, which in
## stock terms is just the stock times the per-capita rate: ageing in
## r_age * J, ageing out and death re * J, batch clearance pbl * J, newborns
## nothing acquired, and the liver-stage clock after the day, as for the people.
## A bite moves the infected up a batch, so the stock carried up is
## sh_v * (inf_n / N) * J, where inf_n / N is a daily probability and cannot
## blow up.
##
## Boosting is by infection, once the refractory window since the last boost
## has passed: boost_immunity() fires when (timestep - last_boosted) >= u, so a
## boost blocks the next u_eff = ceiling(u) - 1 days. Who is still inside their
## window is carried as a stock of its own, R, per cell -- the people boosted in
## the last u_eff days -- and moved with the people like any other: aged, dying,
## losing batches, moved up a batch by a bite, taken by radical cure, and
## released when the window ends. A day's infections are then split by it: the
## share R / N of them fall inside a window and do not boost, the rest do.
##
## A renewal rate p / (p u + 1), the same in every cell, gets the share right for
## someone whose infection risk has been steady for the whole window, and wrong
## for vivax, because a bite infection boosts AND moves the person up a batch,
## into a cell whose relapse hazard is higher by f. A child with one batch
## relapses within IAA's 44-day window two times in three; those relapses fall
## inside the window of the bite that brought the batch, and do not boost. The
## renewal rate cannot see that the cell's newcomers are all freshly boosted, and
## put young children's IAA 10-11% above the IBM's. The stock sees it.
##
## The window is a chain of n stages, each passed on with probability n / u_eff
## a day, on a clock applied after the day's other moves (as the liver-stage one
## is), with the day's boosts joining its first stage after it: so a boost made
## today first counts tomorrow, and n = u_eff stages hold everyone exactly u_eff
## days, as the IBM does. Refractory people are infected at their cell's
## average rate, inf_n / N. A boosted person also skips the day's decay (their
## update replaces it), so the decay is taken back for the immunity they
## carried: boosts x (J / N) x (1 - e).
dim(RAv) <- c(n_age_v, n_het_v, n_hyp, n_ra)
dim(RCv) <- c(n_age_v, n_het_v, n_hyp, n_rc)
p_ra <- n_ra / max(ua_eff, 1)
p_rc <- n_rc / max(uc_eff, 1)
RA_tot[, , ] <- sum(RAv[i, j, k, ])
RC_tot[, , ] <- sum(RCv[i, j, k, ])
# the share of a cell's infections that boost; min(u_eff, 1) switches the stock
# off for a zero window
tA[, , ] <- 1 - min(ua_eff, 1) * min(RA_tot[i, j, k] / (Npop_v[i, j, k] + pop_floor), 1)
tC[, , ] <- 1 - min(uc_eff, 1) * min(RC_tot[i, j, k] / (Npop_v[i, j, k] + pop_floor), 1)
bstA_n[, , ] <- tA[i, j, k] * inf_n[i, j, k]
bstC_n[, , ] <- tC[i, j, k] * inf_n[i, j, k]
finf_v[, , ] <- inf_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)
fbA_v[, , ] <- bstA_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)   # boosts per head
fbC_v[, , ] <- bstC_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)
dim(RA_tot, RC_tot, tA, tC, bstA_n, bstC_n, finf_v, fbA_v, fbC_v) <- c(n_age_v, n_het_v, n_hyp)

bstA_in[, , ] <- bstA_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bstA_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
bstC_in[, , ] <- bstC_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bstC_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
upA_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * JAv[i, j, k]
upC_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * JCv[i, j, k]
dim(bstA_in, bstC_in, upA_v, upC_v) <- c(n_age_v, n_het_v, n_hyp)
dim(qJA, qJC, qKA, qKC, pJAv, pJCv, pKAv, pKCv) <- c(n_age_v, n_het_v, n_hyp)

## Radical cure takes its people's immunity to k_rc with them: what they carried
## from the cell they were infected in -- in place, or a batch down for a bite --
## plus the boosts that infection gave them. They are the day's clinical cases,
## and the clinical draw is taken node by node, so they carry the immunity of the
## nodes that fell ill rather than their cell's mean: the less immune fall ill
## more. cw_rc is each node's clinical cases and nA_rc, nC_rc its acquired
## immunity, m u^3. The node sums are taken relative to the quadrature's own, so
## that a draw taking every node alike carries exactly the cell's J / N and K / N
## per head, however far the Wilson-Hilferty nodes' moments are from the gamma's.
## fC_v is the clinical share of a cell, so as with finf_v nothing divides a flow
## by a population that could be empty.
fC_v[, , ] <- toC_n[i, j, k] / (Npop_v[i, j, k] + pop_floor)
cBA_v[, , ] <- tA[i, j, k] * toC_n[i, j, k]
cBC_v[, , ] <- tC[i, j, k] * toC_n[i, j, k]
cw_rc[, , , ] <- iS_v[i, j, k] * Sv[i, j, k] * w_lmc[i, j, k, l] + Uv[i, j, k] * w_iulmc[i, j, k, l] +
  iA_v[i, j, k] * Av[i, j, k] * w_c[i, j, k, l]
nA_rc[, , , ] <- mA[i, j, k] * max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0) *
  max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0) * max(whaA[i, j, k] + whbA[i, j, k] * qz[l], 0)
nC_rc[, , , ] <- mC[i, j, k] * max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0) *
  max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0) * max(whaC[i, j, k] + whbC[i, j, k] * qz[l], 0)
swA1[, , , ] <- cw_rc[i, j, k, l] * nA_rc[i, j, k, l]
swC1[, , , ] <- cw_rc[i, j, k, l] * nC_rc[i, j, k, l]
swA2[, , , ] <- swA1[i, j, k, l] * nA_rc[i, j, k, l]
swC2[, , , ] <- swC1[i, j, k, l] * nC_rc[i, j, k, l]
sqA1[, , , ] <- qw[l] * nA_rc[i, j, k, l]
sqC1[, , , ] <- qw[l] * nC_rc[i, j, k, l]
sqA2[, , , ] <- sqA1[i, j, k, l] * nA_rc[i, j, k, l]
sqC2[, , , ] <- sqC1[i, j, k, l] * nC_rc[i, j, k, l]
dim(cw_rc, nA_rc, nC_rc, swA1, swC1, swA2, swC2, sqA1, sqC1, sqA2, sqC2) <-
  c(n_age_rc, n_het_rc, n_hyp_rc, n_q)
clJA[, , ] <- mA[i, j, k] * sum(swA1[i, j, k, ]) / (sum(sqA1[i, j, k, ]) + pop_floor)
clJC[, , ] <- mC[i, j, k] * sum(swC1[i, j, k, ]) / (sum(sqC1[i, j, k, ]) + pop_floor)
clKA[, , ] <- KAv[i, j, k] / (Npop_v[i, j, k] + pop_floor) * sum(swA2[i, j, k, ]) /
  (sum(sqA2[i, j, k, ]) + pop_floor)
clKC[, , ] <- KCv[i, j, k] / (Npop_v[i, j, k] + pop_floor) * sum(swC2[i, j, k, ]) /
  (sum(sqC2[i, j, k, ]) + pop_floor)
cJA_n[, , ] <- clJA[i, j, k] + cBA_v[i, j, k]
cJC_n[, , ] <- clJC[i, j, k] + cBC_v[i, j, k]
dim(fC_v, cBA_v, cBC_v, clJA, clJC, clKA, clKC, cJA_n, cJC_n) <- c(n_age_rc, n_het_rc, n_hyp_rc)
rcJA_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cJA_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cJA_n[i, j, k - 1] else 0))
rcJC_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cJC_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cJC_n[i, j, k - 1] else 0))
dim(rcJA_v, rcJC_v) <- c(n_age_rc, n_het_rc, n_hyp_rc)
rcJA_tot[, ] <- sum(rcJA_v[i, j, ])
rcJC_tot[, ] <- sum(rcJC_v[i, j, ])
dim(rcJA_tot, rcJC_tot) <- c(n_age_rc, n_het_rc)

qJA[, , ] <- JAv[i, j, k] + bstA_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * JAv[i - 1, j, k] else 0) - re[i] * JAv[i, j, k] +
  (if (k > 1) upA_v[i, j, k - 1] else 0) - upA_v[i, j, k] -
  (if (rc_on == 1) rcJA_v[i, j, k] else 0) -
  (1 - fbA_v[i, j, k]) * JAv[i, j, k] * (1 - e1_iaa)
pJAv[, , ] <- qJA[i, j, k] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qJA[i, j, k + 1] else 0)
update(JAv[, , ]) <- pJAv[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pJAv[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pJAv[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) rcJA_tot[i, j] else 0) else 0)
qJC[, , ] <- JCv[i, j, k] + bstC_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * JCv[i - 1, j, k] else 0) - re[i] * JCv[i, j, k] +
  (if (k > 1) upC_v[i, j, k - 1] else 0) - upC_v[i, j, k] -
  (if (rc_on == 1) rcJC_v[i, j, k] else 0) -
  (1 - fbC_v[i, j, k]) * JCv[i, j, k] * (1 - e1_ica)
pJCv[, , ] <- qJC[i, j, k] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qJC[i, j, k + 1] else 0)
update(JCv[, , ]) <- pJCv[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pJCv[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pJCv[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) rcJC_tot[i, j] else 0) else 0)

## The refractory stocks, stage by stage. The day's moves first: ageing, death,
## a bite moving the refractory infected up a batch, and radical cure taking the
## refractory among the day's clinical away. Then the window's clock, the day's
## batch clearance and the liver-stage clock; then the radically cured arrive at
## k_rc, and the day's boosts join the first stage -- the radically cured
## among them at k_rc.
upRA[, , , ] <- sh_v[i, j, k] * finf_v[i, j, k] * RAv[i, j, k, l]
upRC[, , , ] <- sh_v[i, j, k] * finf_v[i, j, k] * RCv[i, j, k, l]
dim(upRA, mRA, cRA, qRA, dRA) <- c(n_age_v, n_het_v, n_hyp, n_ra)
dim(upRC, mRC, cRC, qRC, dRC) <- c(n_age_v, n_het_v, n_hyp, n_rc)
# radical cure's share of each cell's refractory people (the clinical share of
# the cell, as with fC_v), landing in place or a batch up as the people do, by
# stage; and of the day's boosts, which start their window at k_rc instead
cRA[, , , ] <- if (rc_on == 1) ft * hyp_mix * ((1 - sh_v[i, j, k]) * fC_v[i, j, k] * RAv[i, j, k, l] +
  (if (k > 1) sh_v[i, j, k - 1] * fC_v[i, j, k - 1] * RAv[i, j, k - 1, l] else 0)) else 0
cRC[, , , ] <- if (rc_on == 1) ft * hyp_mix * ((1 - sh_v[i, j, k]) * fC_v[i, j, k] * RCv[i, j, k, l] +
  (if (k > 1) sh_v[i, j, k - 1] * fC_v[i, j, k - 1] * RCv[i, j, k - 1, l] else 0)) else 0
rcBA[, , ] <- if (rc_on == 1) ft * hyp_mix * ((1 - sh_v[i, j, k]) * cBA_v[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cBA_v[i, j, k - 1] else 0)) else 0
rcBC[, , ] <- if (rc_on == 1) ft * hyp_mix * ((1 - sh_v[i, j, k]) * cBC_v[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cBC_v[i, j, k - 1] else 0)) else 0
dim(rcBA, rcBC) <- c(n_age_v, n_het_v, n_hyp)
cRA_tot[, , ] <- sum(cRA[i, j, , k])
cRC_tot[, , ] <- sum(cRC[i, j, , k])
dim(cRA_tot) <- c(n_age_v, n_het_v, n_ra)
dim(cRC_tot) <- c(n_age_v, n_het_v, n_rc)
rcBA_tot[, ] <- sum(rcBA[i, j, ])
rcBC_tot[, ] <- sum(rcBC[i, j, ])
dim(rcBA_tot, rcBC_tot) <- c(n_age_v, n_het_v)
mRA[, , , ] <- RAv[i, j, k, l] +
  (if (i > 1) r_age[i - 1] * RAv[i - 1, j, k, l] else 0) - re[i] * RAv[i, j, k, l] +
  (if (k > 1) upRA[i, j, k - 1, l] else 0) - upRA[i, j, k, l] - cRA[i, j, k, l]
mRC[, , , ] <- RCv[i, j, k, l] +
  (if (i > 1) r_age[i - 1] * RCv[i - 1, j, k, l] else 0) - re[i] * RCv[i, j, k, l] +
  (if (k > 1) upRC[i, j, k - 1, l] else 0) - upRC[i, j, k, l] - cRC[i, j, k, l]
qRA[, , , ] <- mRA[i, j, k, l] * (1 - p_ra) + (if (l > 1) p_ra * mRA[i, j, k, l - 1] else 0)
qRC[, , , ] <- mRC[i, j, k, l] * (1 - p_rc) + (if (l > 1) p_rc * mRC[i, j, k, l - 1] else 0)
dRA[, , , ] <- qRA[i, j, k, l] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qRA[i, j, k + 1, l] else 0)
dRC[, , , ] <- qRC[i, j, k, l] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qRC[i, j, k + 1, l] else 0)
# the radically cured refractory arrive at k_rc after the liver-stage clock, as
# the people do, their windows having run on through the day
wcRA[, , ] <- cRA_tot[i, j, k] * (1 - p_ra) + (if (k > 1) p_ra * cRA_tot[i, j, k - 1] else 0)
wcRC[, , ] <- cRC_tot[i, j, k] * (1 - p_rc) + (if (k > 1) p_rc * cRC_tot[i, j, k - 1] else 0)
dim(wcRA) <- c(n_age_v, n_het_v, n_ra)
dim(wcRC) <- c(n_age_v, n_het_v, n_rc)
update(RAv[, , , ]) <- dRA[i, j, k, l] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * dRA[i, j, k - 1, l] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * dRA[i, j, n_hyp, l] else 0) +
  (if (k == k_rc) wcRA[i, j, l] else 0) +
  (if (l == 1) min(ua_eff, 1) * (bstA_in[i, j, k] - rcBA[i, j, k] +
    (if (k == k_rc) rcBA_tot[i, j] else 0)) else 0)
update(RCv[, , , ]) <- dRC[i, j, k, l] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * dRC[i, j, k - 1, l] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * dRC[i, j, n_hyp, l] else 0) +
  (if (k == k_rc) wcRC[i, j, l] else 0) +
  (if (l == 1) min(uc_eff, 1) * (bstC_in[i, j, k] - rcBC[i, j, k] +
    (if (k == k_rc) rcBC_tot[i, j] else 0)) else 0)

## Second moments, K = sum of I^2 over a cell. Every flow of people carries its
## share of K exactly as it carries its share of J -- which is exact, not an
## approximation, because nothing that moves a person between cells or infects
## them depends on their own immunity once their cell is known. What differs is
## what a boost and a day of decay do to I^2: a boost takes I to I + 1 (from the
## undecayed value), adding 2 I + 1, and a day of decay multiplies I^2 by
## exp(-2/d). Boosted people are a random draw of their cell, so b boosts carry
## b J / N of immunity and add 2 b J / N to K, and skip the decay of b K / N.
##
## The rest of a boost's contribution is the spread it injects, and that is NOT
## one unit per boost. The refractory window makes a person's boosts a renewal
## process -- u_eff blocked days, then a geometric wait at daily probability p --
## far more regular than independent events: over a long stretch its count
## varies by the gap's squared CV times its mean, the Fano factor (1 - p) t^2,
## t being the share of infections that boost (1 / (p u + 1) at a steady p).
## For IAA's 44-day window that is ~0.06, so counting each boost as a unit of
## variance, as a Markov jump would, overstates the IBM's within-cell spread of
## IAA by 20-45%. The injection is therefore boosts x Fano, p (1 - p) t^3 per
## person. A boost also moves its person's mean, by the boost probability p t a
## day, which adds its square to I^2 over and above the variance: (p t)^2 per
## person. ICA carries both terms: with no window they give exactly one boost's
## worth per person-day, and simulated person by person against the IBM's rule
## they hold ICA's within-cell spread within 6% of the IBM's up to p = 0.2 a day.
## IAA's long window makes its boosts so regular that the second term
## overshoots, so IAA keeps the Fano term alone. That matches the IBM's spread
## at EIR 10 and below, and departs from it as p rises past about 0.05 a day
## (V.3).
vA_n[, , ] <- iS_v[i, j, k] * (1 - iS_v[i, j, k]) * Sv[i, j, k] +
  iU_v[i, j, k] * (1 - iU_v[i, j, k]) * Uv[i, j, k] +
  iA_v[i, j, k] * (1 - iA_v[i, j, k]) * Av[i, j, k] +
  iD_v[i, j, k] * (1 - iD_v[i, j, k]) * Dv[i, j, k]
vB_n[, , ] <- iS_v[i, j, k] * iS_v[i, j, k] * Sv[i, j, k] +
  iU_v[i, j, k] * iU_v[i, j, k] * Uv[i, j, k] +
  iA_v[i, j, k] * iA_v[i, j, k] * Av[i, j, k] +
  iD_v[i, j, k] * iD_v[i, j, k] * Dv[i, j, k]
bvA_n[, , ] <- tA[i, j, k] * tA[i, j, k] * tA[i, j, k] * vA_n[i, j, k]
bvC_n[, , ] <- tC[i, j, k] * tC[i, j, k] * (tC[i, j, k] * vA_n[i, j, k] + vB_n[i, j, k])
dim(vA_n, vB_n, bvA_n, bvC_n) <- c(n_age_v, n_het_v, n_hyp)
bKA_n[, , ] <- 2 * fbA_v[i, j, k] * JAv[i, j, k] + bvA_n[i, j, k]
bKC_n[, , ] <- 2 * fbC_v[i, j, k] * JCv[i, j, k] + bvC_n[i, j, k]
bKA_in[, , ] <- bKA_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bKA_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
bKC_in[, , ] <- bKC_n[i, j, k] * (1 - sh_v[i, j, k]) +
  (if (k > 1) bKC_n[i, j, k - 1] * sh_v[i, j, k - 1] else 0)
upKA_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * KAv[i, j, k]
upKC_v[, , ] <- sh_v[i, j, k] * finf_v[i, j, k] * KCv[i, j, k]
dim(bKA_n, bKC_n, bKA_in, bKC_in, upKA_v, upKC_v) <- c(n_age_v, n_het_v, n_hyp)
## radical cure moves the spread its people's boosts injected with them, by
## the same Fano rule, from each clinical source
cv_n[, , ] <- iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * (1 - iS_v[i, j, k]) +
  s_iulmc[i, j, k] * Uv[i, j, k] * (1 - iU_v[i, j, k]) +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * (1 - iA_v[i, j, k])
cvB_n[, , ] <- iS_v[i, j, k] * s_lmc[i, j, k] * Sv[i, j, k] * iS_v[i, j, k] +
  s_iulmc[i, j, k] * Uv[i, j, k] * iU_v[i, j, k] +
  iA_v[i, j, k] * s_c[i, j, k] * Av[i, j, k] * iA_v[i, j, k]
cvA_v[, , ] <- tA[i, j, k] * tA[i, j, k] * tA[i, j, k] * cv_n[i, j, k]
cvC_v[, , ] <- tC[i, j, k] * tC[i, j, k] * (tC[i, j, k] * cv_n[i, j, k] + cvB_n[i, j, k])
cKA_n[, , ] <- clKA[i, j, k] + 2 * tA[i, j, k] * clJA[i, j, k] + cvA_v[i, j, k]
cKC_n[, , ] <- clKC[i, j, k] + 2 * tC[i, j, k] * clJC[i, j, k] + cvC_v[i, j, k]
rcKA_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cKA_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cKA_n[i, j, k - 1] else 0))
rcKC_v[, , ] <- ft * hyp_mix * ((1 - sh_v[i, j, k]) * cKC_n[i, j, k] +
  (if (k > 1) sh_v[i, j, k - 1] * cKC_n[i, j, k - 1] else 0))
dim(cv_n, cvB_n, cvA_v, cvC_v, cKA_n, cKC_n, rcKA_v, rcKC_v) <- c(n_age_rc, n_het_rc, n_hyp_rc)
rcKA_tot[, ] <- sum(rcKA_v[i, j, ])
rcKC_tot[, ] <- sum(rcKC_v[i, j, ])
dim(rcKA_tot, rcKC_tot) <- c(n_age_rc, n_het_rc)

qKA[, , ] <- KAv[i, j, k] + bKA_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * KAv[i - 1, j, k] else 0) - re[i] * KAv[i, j, k] +
  (if (k > 1) upKA_v[i, j, k - 1] else 0) - upKA_v[i, j, k] -
  (if (rc_on == 1) rcKA_v[i, j, k] else 0) -
  (1 - fbA_v[i, j, k]) * KAv[i, j, k] * (1 - e2_iaa)
pKAv[, , ] <- qKA[i, j, k] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qKA[i, j, k + 1] else 0)
update(KAv[, , ]) <- pKAv[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pKAv[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pKAv[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) rcKA_tot[i, j] else 0) else 0)
qKC[, , ] <- KCv[i, j, k] + bKC_in[i, j, k] +
  (if (i > 1) r_age[i - 1] * KCv[i - 1, j, k] else 0) - re[i] * KCv[i, j, k] +
  (if (k > 1) upKC_v[i, j, k - 1] else 0) - upKC_v[i, j, k] -
  (if (rc_on == 1) rcKC_v[i, j, k] else 0) -
  (1 - fbC_v[i, j, k]) * KCv[i, j, k] * (1 - e2_ica)
pKCv[, , ] <- qKC[i, j, k] * (1 - pbl[k] * fdec_v[i, j, k]) +
  (if (k < n_hyp) pbl[k + 1] * fdec_v[i, j, k + 1] * qKC[i, j, k + 1] else 0)
update(KCv[, , ]) <- pKCv[i, j, k] * (1 - rls_k * lsmask[k]) +
  (if (k > n_bat + 1) rls_k * pKCv[i, j, k - 1] else 0) +
  (if (k == 1) rls_k * lsmask[n_hyp] * pKCv[i, j, n_hyp] else 0) +
  (if (rc_on == 1) (if (k == k_rc) rcKC_tot[i, j] else 0) else 0)

## vivax onward infectivity: A is LM-detectable by definition for vivax and
## carries the constant ca, where falciparum's asymptomatic infectivity depends
## on ID
inf_v[, , ] <- cD * Dv[i, j, k] * tbv_fD[i] + cA_v * Av[i, j, k] * tbv_fA[i] +
  cU * Uv[i, j, k] * tbv_fU[i] + cT * (Trv[i, j, k] + Trv_slow[i, j, k]) * tbv_fT[i]
infw_v[, , ] <- zeta[j] * psi[i] * inf_v[i, j, k]
dim(inf_v, infw_v) <- c(n_age_v, n_het_v, n_hyp)

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
zp_n[, ] <- zeta[j] * psi[i] * (Npop[i, j] + nv_ij[i, j])
dim(psi_n) <- n_age
dim(zp_n) <- c(n_age, n_het)
mpsi <- sum(psi_n) / sum(n_g_now)
mzp <- sum(zp_n) / sum(n_g_now)
inf_now <- (sum(infw) + sum(infw_v)) / mzp       # both blocks; one is empty
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
# The vivax block's share of the same outputs, and its own two. Vivax A is
# LM-detectable by definition, so LM counts D, Tr and all of A, as
# create_prevalence_renderer() does for vivax. Relapses are the relapse share of
# the day's infections, relapse_rates / infection_rates as in
# calculate_vivax_infections(). Severe disease does not exist for vivax.
hypmask[] <- if (i > 1 && i <= n_bat) 1 else 0     # levels holding at least one batch
dim(hypmask) <- n_hyp
detlm_v[, , ] <- Dv[i, j, k] + Trv[i, j, k] + Trv_slow[i, j, k] + Av[i, j, k]
detpcr_v[, , ] <- detlm_v[i, j, k] + Uv[i, j, k]
rel_v[, , ] <- inf_n[i, j, k] * kk[k] * ff / (r_tot_v[i, j, k] + pop_floor)
hyp_v[, , ] <- hypmask[k] * Npop_v[i, j, k]
dim(detlm_v, detpcr_v, rel_v, hyp_v) <- c(n_age_v, n_het_v, n_hyp)
dlm_va[] <- sum(detlm_v[i, , ])
dpcr_va[] <- sum(detpcr_v[i, , ])
clin_va[] <- sum(arr_C[i, , ])
inc_va[] <- sum(inf_n[i, , ])
rel_va[] <- sum(rel_v[i, , ])
hyp_va[] <- sum(hyp_v[i, , ])
S_va[] <- sum(Sv[i, , ])
D_va[] <- sum(Dv[i, , ])
A_va[] <- sum(Av[i, , ])
U_va[] <- sum(Uv[i, , ])
Tr_va[] <- sum(Trv[i, , ]) + sum(Trv_slow[i, , ])
Ph_va[] <- sum(Phv_tot[i, , ])
dim(dlm_va, dpcr_va, clin_va, inc_va, rel_va, hyp_va) <- n_age_v
dim(S_va, D_va, A_va, U_va, Tr_va, Ph_va) <- n_age_v
update(n_g[]) <- n_g_now[i]
update(det_lm_g[]) <- sum(detlm[i, ]) + (if (i <= n_age_v) dlm_va[i] else 0)
update(det_pcr_g[]) <- sum(detpcr[i, ]) + (if (i <= n_age_v) dpcr_va[i] else 0)
update(clin_g[]) <- sum(clin_n[i, ]) + (if (i <= n_age_v) clin_va[i] else 0)
update(sev_g[]) <- sum(sev_n[i, ])
update(inc_g[]) <- sum(inf_tot[i, ]) + (if (i <= n_age_v) inc_va[i] else 0)
update(relapse_g[]) <- if (i <= n_age_v) rel_va[i] else 0
update(hyp_g[]) <- if (i <= n_age_v) hyp_va[i] else 0
update(S_g[]) <- sum(S[i, ]) + (if (i <= n_age_v) S_va[i] else 0)
update(D_g[]) <- sum(D[i, ]) + (if (i <= n_age_v) D_va[i] else 0)
update(A_g[]) <- sum(A[i, ]) + (if (i <= n_age_v) A_va[i] else 0)
update(U_g[]) <- sum(U[i, ]) + (if (i <= n_age_v) U_va[i] else 0)
update(Tr_g[]) <- sum(Tr[i, ]) + sum(Tr_slow[i, ]) + sum(Trc_tot[i, ]) +
  (if (i <= n_age_v) Tr_va[i] else 0)
update(Ph_g[]) <- sum(Ph_tot[i, ]) + sum(Phc_tot[i, ]) + (if (i <= n_age_v) Ph_va[i] else 0)
update(EIR_yr) <- eir_lag * 365            # the EIR biting humans today, as EIR_<species>
update(FOIM) <- foim[1]
update(ft_out) <- ft
# Immunity summed over the people in each age group, for the IBM's
# population means (ica_mean and the rest) and their age bands. The falciparum
# block holds cell means, the vivax block stocks; both blocks add to the two
# they share, one of them empty.
imm_ica[, ] <- ICA[i, j] * Npop[i, j]
imm_icm[, ] <- ICM[i, j] * Npop[i, j]
imm_ib[, ] <- IB[i, j] * Npop[i, j]
imm_iva[, ] <- IVA[i, j] * Npop[i, j]
imm_ivm[, ] <- IVM[i, j] * Npop[i, j]
imm_id[, ] <- ID[i, j] * Npop[i, j]
dim(imm_ica, imm_icm, imm_ib, imm_iva, imm_ivm, imm_id) <- c(n_age, n_het)
imm_iamv[, ] <- IAMv[i, j] * Nv_ij[i, j]
imm_icmv[, ] <- ICMv[i, j] * Nv_ij[i, j]
dim(imm_iamv, imm_icmv) <- c(n_age_v, n_het_v)
hypk_v[, , ] <- kk[k] * Npop_v[i, j, k]
dim(hypk_v) <- c(n_age_v, n_het_v, n_hyp)
update(ica_g[]) <- sum(imm_ica[i, ]) + (if (i <= n_age_v) sum(JCv[i, , ]) else 0)
update(icm_g[]) <- sum(imm_icm[i, ]) + (if (i <= n_age_v) sum(imm_icmv[i, ]) else 0)
update(ib_g[]) <- sum(imm_ib[i, ])
update(iva_g[]) <- sum(imm_iva[i, ])
update(ivm_g[]) <- sum(imm_ivm[i, ])
update(id_g[]) <- sum(imm_id[i, ])
update(iaa_g[]) <- if (i <= n_age_v) sum(JAv[i, , ]) else 0
update(iam_g[]) <- if (i <= n_age_v) sum(imm_iamv[i, ]) else 0
update(hypk_g[]) <- if (i <= n_age_v) sum(hypk_v[i, , ]) else 0
dim(n_g, det_lm_g, det_pcr_g, clin_g, sev_g, inc_g, relapse_g, hyp_g) <- n_age
dim(S_g, D_g, A_g, U_g, Tr_g, Ph_g) <- n_age
dim(ica_g, icm_g, ib_g, iva_g, ivm_g, id_g, iaa_g, iam_g, hypk_g) <- n_age

## ---- initial conditions ---------------------------------------------------
S0 <- parameter(); D0 <- parameter(); A0 <- parameter()
U0 <- parameter(); Tr0 <- parameter(); Ph0 <- parameter(); Phc0 <- parameter()
IB_init <- parameter(); ICA_init <- parameter(); ID_init <- parameter(); IVA_init <- parameter()
dim(S0, D0, A0, U0, Tr0, IB_init, ICA_init, ID_init, IVA_init) <- c(n_age, n_het)
dim(Ph0) <- c(n_age, n_het, n_ph)
dim(Phc0) <- c(n_age, n_het, n_phc)
Sv0 <- parameter(); Dv0 <- parameter(); Av0 <- parameter(); Uv0 <- parameter()
Trv0 <- parameter(); JAv0 <- parameter(); JCv0 <- parameter()
KAv0 <- parameter(); KCv0 <- parameter()
dim(Sv0, Dv0, Av0, Uv0, Trv0, JAv0, JCv0, KAv0, KCv0) <- c(n_age_v, n_het_v, n_hyp)
Phv0 <- parameter(); dim(Phv0) <- c(n_age_v, n_het_v, n_hyp, n_phv)
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
initial(Sv[, , ]) <- Sv0[i, j, k]
initial(Dv[, , ]) <- Dv0[i, j, k]
initial(Av[, , ]) <- Av0[i, j, k]
initial(Uv[, , ]) <- Uv0[i, j, k]
initial(Trv[, , ]) <- (1 - spc0) * Trv0[i, j, k]
initial(Trv_slow[, , ]) <- spc0 * Trv0[i, j, k]
initial(Phv[, , , ]) <- Phv0[i, j, k, l]
initial(JAv[, , ]) <- JAv0[i, j, k]
initial(JCv[, , ]) <- JCv0[i, j, k]
initial(KAv[, , ]) <- KAv0[i, j, k]
initial(KCv[, , ]) <- KCv0[i, j, k]
# nobody is inside a refractory window at the start: the IBM's last_boosted_*
# start at -1
initial(RAv[, , , ]) <- 0
initial(RCv[, , , ]) <- 0
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
initial(relapse_g[]) <- 0
initial(hyp_g[]) <- 0
initial(S_g[]) <- 0
initial(D_g[]) <- 0
initial(A_g[]) <- 0
initial(U_g[]) <- 0
initial(Tr_g[]) <- 0
initial(Ph_g[]) <- 0
initial(EIR_yr) <- 0
initial(FOIM) <- 0
initial(ft_out) <- 0
initial(ica_g[]) <- 0
initial(icm_g[]) <- 0
initial(ib_g[]) <- 0
initial(iva_g[]) <- 0
initial(ivm_g[]) <- 0
initial(id_g[]) <- 0
initial(iaa_g[]) <- 0
initial(iam_g[]) <- 0
initial(hypk_g[]) <- 0
