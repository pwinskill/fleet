# Phase 2 intervention time-series builders: clinical treatment ft(t),
# antimalarial resistance (ETF/SPC), drug-linked treatment properties
# (efficacy, treated infectivity, prophylaxis duration), and the mean-field
# vector-control coefficients a(t)/mu(t) per species.
#
# All series start at their baseline value at t = 0 so the equilibrium seed
# (which has no interventions) is preserved; interventions act from their
# scheduled timesteps onward, exactly as malariasimulation applies them.

.weibull_mean <- function(shape, scale) scale * gamma(1 + 1 / shape)

# Coverage active at timestep t for one drug's (timesteps, coverages).
.cov_at <- function(ts, cov, t) {
  idx <- which(ts <= t)
  if (length(idx)) cov[[max(idx)]] else 0
}

# ---- clinical treatment ft(t): summed coverage across drugs ----------------
treatment_series <- function(p) {
  cov <- p$clinical_treatment_coverages
  if (is.null(cov) || length(cov) == 0) {
    return(list(times = 0, vals = 0))
  }
  change_ts <- sort(unique(c(0, unlist(p$clinical_treatment_timesteps))))
  change_ts <- change_ts[change_ts >= 0]
  vals <- vapply(change_ts, function(t) {
    total <- 0
    for (d in seq_along(cov)) {
      total <- total + .cov_at(p$clinical_treatment_timesteps[[d]], cov[[d]], max(t, 1))
    }
    min(total, 1)
  }, numeric(1))
  list(times = change_ts, vals = vals)
}

# ---- drug-linked treatment properties --------------------------------------
# Coverage-weighted (over each drug's peak scheduled coverage, so the mix is
# defined even for later-onset treatment): drug_eff (efficacy), cT (treated
# infectivity = cd*drug_rel_c), and rP (prophylaxis rate). rP subtracts the
# ~1/rT days already spent refractory in Tr from the Weibull mean protection.
drug_mix <- function(p, eqp) {
  drugs <- p$clinical_treatment_drugs
  if (is.null(drugs) || length(drugs) == 0) {
    return(list(drug_eff = 1, cT = eqp[["cT"]], rP = eqp[["rP"]]))
  }
  di <- w <- numeric(length(drugs))
  for (d in seq_along(drugs)) {
    di[d] <- p$clinical_treatment_drugs[[d]]
    w[d] <- max(p$clinical_treatment_coverages[[d]])   # peak coverage
  }
  if (sum(w) == 0) return(list(drug_eff = 1, cT = eqp[["cT"]], rP = eqp[["rP"]]))
  wn <- w / sum(w)
  drug_eff <- sum(wn * p$drug_efficacy[di])
  cT <- p$cd * sum(wn * p$drug_rel_c[di])
  mean_dur <- sum(wn * .weibull_mean(p$drug_prophylaxis_shape[di],
                                     p$drug_prophylaxis_scale[di]))
  dt <- 1 / eqp[["rT"]]                                # mean days in Tr (refractory)
  rP <- 1 / max(mean_dur - dt, 1e-6)
  list(drug_eff = drug_eff, cT = cT, rP = rP)
}

# ---- antimalarial resistance ETF(t)/SPC(t) and slow-clearance rate ----------
# Blended across every treatment drug that carries a resistance schedule, weighted
# by each drug's share of clinical treatment (the same coverage weighting as
# drug_mix). ETF/SPC affect only the fraction of treated who received a resistant
# drug, so the population early-treatment-failure / slow-clearance rates are the
# coverage-weighted sums over drugs.
resistance_series <- function(p, eqp) {
  rT_base <- eqp[["rT"]]
  none <- list(times = 0, etf = 0, spc = 0, rT_slow = rT_base)
  if (!isTRUE(p$antimalarial_resistance)) return(none)
  cdrugs <- vapply(p$clinical_treatment_drugs, function(x) x[[1]], numeric(1))
  cw <- vapply(p$clinical_treatment_coverages, function(x) max(x), numeric(1))  # peak coverage
  if (length(cdrugs) == 0 || sum(cw) == 0) return(none)
  wn <- cw / sum(cw)
  rdrug <- vapply(p$antimalarial_resistance_drug, function(x) x[[1]], numeric(1))
  all_ts <- sort(unique(c(0, unlist(p$antimalarial_resistance_timesteps))))
  step_at <- function(ts, vals, t) if (t < min(ts)) 0 else vals[max(which(ts <= t))]
  nT <- length(all_ts)
  etf <- spc <- slow_num <- numeric(nT)     # slow_num = sum wn*spc_c*dt_slow_c
  for (ci in seq_along(cdrugs)) {
    k <- which(rdrug == cdrugs[ci])          # resistance entry for this treatment drug
    if (!length(k)) next
    k <- k[1]
    ts   <- p$antimalarial_resistance_timesteps[[k]]
    art  <- p$artemisinin_resistance_proportion[[k]]
    etfp <- p$early_treatment_failure_probability[[k]]
    spcp <- p$slow_parasite_clearance_probability[[k]]
    dts  <- p$dt_slow_parasite_clearance[[k]]
    dts  <- if (length(dts)) dts[length(dts)] else 1 / rT_base
    for (i in seq_len(nT)) {
      t <- all_ts[i]; a <- step_at(ts, art, t)
      e <- a * step_at(ts, etfp, t); s <- a * step_at(ts, spcp, t)
      etf[i] <- etf[i] + wn[ci] * e
      # slow fraction of the treated cohort. The IBM applies SPC to post-ETF
      # survivors; here spc is the slow fraction of the (already ETF-net) Tr inflow
      # -> exact for a single drug (the common case), a small over-count only for
      # multiple drugs at high resistance.
      spc[i] <- spc[i] + wn[ci] * s
      slow_num[i] <- slow_num[i] + wn[ci] * s * dts
    }
  }
  # blended slow-clearance duration at PEAK resistance (single scalar, as modelled).
  # Peak, not the last timestep: a rise-then-fall schedule can end at spc = 0.
  nPeak <- which.max(spc)
  dt_slow_blend <- if (length(nPeak) && spc[nPeak] > 0) slow_num[nPeak] / spc[nPeak] else 1 / rT_base
  list(times = all_ts, etf = etf, spc = spc, rT_slow = 1 / dt_slow_blend)
}

# ---- mean-field vector control: per-species a(t), mu(t) ---------------------
# Reuses malariasimulation's bed-net / IRS decay + combination formulas,
# population-averaged over the net-using / sprayed fractions. Returns arrays in
# SPECIES-MAJOR order [n_spp, n_time] as required by odin2 array interpolate().
.spray_decay <- function(t, theta, gamma) 1 / (1 + exp(-(theta + gamma * t)))

# Fraction of a cohort still using a net `sn` days after distribution.
#   exponential (`bednet_retention` set): exp(-sn / retention).
#   logistic (`bednet_logistic_*`): survival of the malariasimulation logistic
#   net-retention time T = l*sqrt(a/(1+a)), a ~ Exp(k), l = half_life /
#   sqrt(1 - k/(k - log 0.5)). P(T > sn) = exp(-k r^2/(1-r^2)) for r = sn/l < 1,
#   else 0 (verified: S(half_life) = 0.5 exactly).
.net_survival <- function(sn, p) {
  if (!is.null(p$bednet_retention)) return(exp(-sn / p$bednet_retention))
  hl <- p$bednet_logistic_half_life; k <- p$bednet_logistic_k
  l <- hl / sqrt(1 - k / (k - log(0.5)))
  s <- numeric(length(sn)); below <- sn < l
  r <- sn[below] / l
  s[below] <- exp(-k * r^2 / (1 - r^2))
  s
}

vector_control_series <- function(p, timesteps) {
  n_spp <- length(p$species_proportions)
  nets <- isTRUE(p$bednets)
  spray <- isTRUE(p$spraying)
  grid <- sort(unique(c(seq(0, timesteps, by = 10), timesteps,
                        if (nets) p$bednet_timesteps,
                        if (spray) p$spraying_timesteps)))
  grid <- grid[grid >= 0 & grid <= timesteps]
  ng <- length(grid)
  # Net-usage weights (species-independent). At each grid time, every past
  # distribution d contributes weight w_d = coverage_d * P(no later re-net) *
  # retention-survival(t - t_d) -- i.e. the fraction of the population whose most
  # recent net came from d and is still in use. sum(w) is total usage; the net-age
  # mixture (via t - t_d) then decays each distribution's efficacy separately. This
  # reduces to the single-distribution formula and is exact for repeated nets.
  net_w <- vector("list", ng)
  if (nets) {
    bt <- p$bednet_timesteps; bc <- as.numeric(unlist(p$bednet_coverages))
    for (g in seq_len(ng)) {
      m <- which(bt <= grid[g])
      if (!length(m)) next
      cm <- bc[m]; sn <- grid[g] - bt[m]
      isuf <- rev(cumprod(rev(1 - cm)))          # inclusive suffix prod of (1-c)
      suff <- c(isuf[-1], 1)                      # prod_{j>d}(1 - c_j)
      net_w[[g]] <- list(idx = m, w = cm * suff * .net_survival(sn, p), sn = sn)
    }
  }
  a_out <- mum_out <- matrix(0, n_spp, ng)             # species-major
  for (s in seq_len(n_spp)) {
    Q0 <- p$Q0[[s]]; fmr <- p$blood_meal_rates[[s]]; ft_forage <- p$foraging_time[[s]]
    mum_base <- p$mum[[s]]; gono <- 1 / fmr - ft_forage
    phi_b <- if (nets) p$phi_bednets[[s]] else 0
    phi_i <- if (spray) p$phi_indoors[[s]] else 0
    for (g in seq_len(ng)) {
      t <- grid[g]
      sn_bar <- 1; rn_bar <- 0
      if (nets && !is.null(net_w[[g]])) {
        nw <- net_w[[g]]; idx <- nw$idx; w <- nw$w
        decay <- exp(-nw$sn / p$bednet_gamman[idx])          # net-age efficacy decay
        rnv <- (p$bednet_rn[idx, s] - p$bednet_rnm[idx, s]) * decay + p$bednet_rnm[idx, s]
        dnv <- p$bednet_dn0[idx, s] * decay
        snv <- 1 - rnv - dnv
        sn_bar <- sum(w * snv) + (1 - sum(w))                # non-users survive at 1
        rn_bar <- sum(w * rnv)
      }
      # IRS: rs (repel) and ss (survive) are the SAME spray outcome, so the
      # survive-and-not-repelled term is a JOINT per-individual mean. Sprayed
      # protection never expires in the IBM (spray_time is not reset), so the
      # population is a mixture over ALL past rounds weighted by most-recent-spray
      # probability cov_d * prod_{j>d}(1-cov_j) (no retention factor, unlike nets).
      ss_surv <- 1; rs_bar <- 0
      if (spray) {
        m <- which(p$spraying_timesteps <= t)
        if (length(m)) {
          k0 <- p$k0
          cm <- vapply(m, function(r) p$spraying_coverages[[r]], numeric(1))
          isuf <- rev(cumprod(rev(1 - cm))); suff <- c(isuf[-1], 1)
          w <- cm * suff
          ss_surv <- 0; rs_bar <- 0
          for (di in seq_along(m)) {
            r <- m[di]; ss <- t - p$spraying_timesteps[r]
            ls <- .spray_decay(ss, p$spraying_ls_theta[r, s], p$spraying_ls_gamma[r, s])
            ks <- k0 * .spray_decay(ss, p$spraying_ks_theta[r, s], p$spraying_ks_gamma[r, s])
            ms <- .spray_decay(ss, p$spraying_ms_theta[r, s], p$spraying_ms_gamma[r, s])
            js <- 1 - ls - ks; msc <- 1 - ms
            lsp <- ls * msc; ksp <- ks * msc; jsp <- js * msc + ms
            rsv <- (1 - ksp / k0) * (jsp / (lsp + jsp))
            ssv <- ksp / k0
            ss_surv <- ss_surv + w[di] * (1 - rsv) * ssv
            rs_bar  <- rs_bar  + w[di] * rsv
          }
          ss_surv <- ss_surv + (1 - sum(w))            # never-sprayed survive at 1
        }
      }
      # combination (net factor independent of spray; spray joint term = ss_surv)
      psurv <- (1 - phi_i) + phi_b * sn_bar * ss_surv + (phi_i - phi_b) * ss_surv
      prep <- phi_b * (1 - rs_bar) * rn_bar + phi_i * rs_bar
      W <- (1 - Q0) + Q0 * psurv
      Z <- Q0 * prep
      if (W <= 0 || Z >= 1) {
        stop("Vector-control parameters give a degenerate biting rate (W<=0 or ",
             "Z>=1); check dn0+rn<1 and coverage/efficacy inputs.", call. = FALSE)
      }
      f <- 1 / (ft_forage / (1 - Z) + gono)
      a <- (1 - (1 - Q0) / W) * f
      p1_0 <- exp(-mum_base * ft_forage); p2 <- exp(-mum_base * gono)
      p1 <- p1_0 * W / (1 - Z * p1_0)
      a_out[s, g] <- a
      mum_out[s, g] <- -f * log(p1 * p2)
    }
  }
  list(times = grid, a = a_out, mum = mum_out)
}

# ---- chemoprevention pulses (MDA / SMC / PMC) ------------------------------
# Each is a mass drug administration applied at scheduled timesteps to a target
# age band: a fraction (coverage x drug_efficacy) of the targeted population has
# infections cleared and receives prophylaxis. Applied as instantaneous state
# jumps between ODE integration segments. PMC (age-based delivery in the IBM) is
# approximated here as periodic pulses on the target age band.
# Early-treatment-failure fraction for a drug at time t (0 if no resistance on it).
# The IBM applies resistance to MDA/SMC/PMC as well as clinical treatment, so a
# resistant chemoprevention drug clears fewer infections.
.drug_etf_at <- function(p, drug, t) {
  if (!isTRUE(p$antimalarial_resistance)) return(0)
  rd <- vapply(p$antimalarial_resistance_drug, function(x) x[[1]], numeric(1))
  k <- which(rd == drug); if (!length(k)) return(0); k <- k[1]
  ts <- p$antimalarial_resistance_timesteps[[k]]
  if (t < min(ts)) return(0)
  i <- max(which(ts <= t))
  p$artemisinin_resistance_proportion[[k]][i] * p$early_treatment_failure_probability[[k]][i]
}

chemoprevention_events <- function(p, timesteps) {
  ev <- list()
  recycle <- function(x, n) if (length(x) == 1L) rep_len(x, n) else x
  add <- function(on, drug, ts, cov, lo, hi) {
    if (!isTRUE(on) || is.null(ts) || length(ts) == 0) return(invisible())
    n <- length(ts); cov <- recycle(cov, n); lo <- recycle(lo, n); hi <- recycle(hi, n)
    eff <- p$drug_efficacy[drug]
    for (k in seq_len(n)) {
      ev[[length(ev) + 1L]] <<- list(time = ts[k], lo = lo[k], hi = hi[k],
        frac = cov[k] * eff * (1 - .drug_etf_at(p, drug, ts[k])))   # resistance ETF
    }
  }
  add(p$mda, p$mda_drug, p$mda_timesteps, p$mda_coverages, p$mda_min_ages, p$mda_max_ages)
  add(p$smc, p$smc_drug, p$smc_timesteps, p$smc_coverages, p$smc_min_ages, p$smc_max_ages)
  # PMC delivery is age-triggered and continuous in the IBM; pmc_timesteps is a
  # COVERAGE schedule. Approximate as pulses at a ~monthly cadence over each dose
  # age band, using the coverage active at each pulse time (per-group overlap
  # weighting in the pulse keeps the effective dose ~coverage x efficacy).
  if (isTRUE(p$pmc)) {
    pts <- p$pmc_timesteps; pcv <- p$pmc_coverages; drug <- p$pmc_drug
    eff <- p$drug_efficacy[drug]; band <- 30
    grid <- seq(min(pts), timesteps, by = band)
    for (t in grid) {
      idx <- which(pts <= t)
      cv <- if (length(idx)) pcv[max(idx)] else 0
      if (cv <= 0) next
      fr <- cv * eff * (1 - .drug_etf_at(p, drug, t))
      for (a in p$pmc_ages)
        ev[[length(ev) + 1L]] <- list(time = t, frac = fr, lo = a, hi = a + band)
    }
  }
  ev
}

# Prophylaxis rate for the shared chemoprevention compartment Ph_c, from the
# MDA/SMC/PMC drug(s)' Weibull mean. Distinct from the clinical-treatment rP. When
# several chemoprevention types with different drugs are co-deployed, the single
# Ph_c decay rate is the coverage-weighted mean of their protection durations
# (better than an arbitrary priority ordering).
chemoprevention_prophylaxis_rate <- function(p, eqp) {
  tot_cov <- function(cov) if (is.null(cov)) 0 else sum(unlist(cov))
  active <- list()
  if (isTRUE(p$smc)) active[[length(active) + 1L]] <- c(p$smc_drug, tot_cov(p$smc_coverages))
  if (isTRUE(p$mda)) active[[length(active) + 1L]] <- c(p$mda_drug, tot_cov(p$mda_coverages))
  if (isTRUE(p$pmc)) active[[length(active) + 1L]] <- c(p$pmc_drug, tot_cov(p$pmc_coverages))
  if (!length(active)) return(eqp[["rP"]])
  w <- vapply(active, `[`, numeric(1), 2)
  dur <- vapply(active, function(a) .weibull_mean(p$drug_prophylaxis_shape[a[1]],
                                                  p$drug_prophylaxis_scale[a[1]]), numeric(1))
  mean_dur <- if (sum(w) > 0) sum(w * dur) / sum(w) else mean(dur)
  1 / mean_dur
}

# Apply a single chemoprevention pulse to the packed dust2 state in place: a
# fraction `frac` of the targeted age band is cleared of infection and moved to
# the chemoprevention prophylaxis compartment Ph_c (protected). (The brief
# treated-infectious phase of detectable cases is omitted; negligible for the
# fast-clearing MDA/SMC/PMC drugs.)
apply_chemoprevention_pulse <- function(sys, uidx, meta, event) {
  n_age <- meta$n_age; n_het <- meta$n_het
  # per-group fraction of the pulse: coverage x efficacy x (fraction of the group
  # that overlaps the target band). Weighting by overlap (not midpoint membership)
  # stops narrow bands from being silently dropped and coarse groups from being
  # fully treated when only partly targeted.
  ov <- pmax(0, pmin(event$hi, meta$age_hi) - pmax(event$lo, meta$age_lo)) /
        (meta$age_hi - meta$age_lo)
  arows <- which(ov > 0)
  if (length(arows) == 0) return(invisible())
  fr <- event$frac * ov[arows]                       # length(arows) vector
  sv <- dust2::dust_system_state(sys)
  gv <- function(nm) matrix(sv[uidx[[nm]]], n_age, n_het)
  S <- gv("S"); D <- gv("D"); A <- gv("A"); U <- gv("U"); Tr <- gv("Tr"); Phc <- gv("Ph_c")
  keep <- 1 - fr
  cleared <- fr * (S[arows, , drop = FALSE] + U[arows, , drop = FALSE] +
                   A[arows, , drop = FALSE] + D[arows, , drop = FALSE] +
                   Tr[arows, , drop = FALSE])
  S[arows, ]  <- S[arows, , drop = FALSE]  * keep
  U[arows, ]  <- U[arows, , drop = FALSE]  * keep
  A[arows, ]  <- A[arows, , drop = FALSE]  * keep
  D[arows, ]  <- D[arows, , drop = FALSE]  * keep
  Tr[arows, ] <- Tr[arows, , drop = FALSE] * keep
  Phc[arows, ] <- Phc[arows, , drop = FALSE] + cleared
  sv[uidx[["S"]]] <- S; sv[uidx[["U"]]] <- U; sv[uidx[["A"]]] <- A
  sv[uidx[["D"]]] <- D; sv[uidx[["Tr"]]] <- Tr; sv[uidx[["Ph_c"]]] <- Phc
  dust2::dust_system_set_state(sys, sv)
  invisible()
}

# ---- vaccines: PEV (pre-erythrocytic) and TBV (transmission-blocking) -------
# Mean-field: point-estimate antibody trajectories (means of the profile
# distributions), giving a per-(age,time) FOI reduction (PEV) or infectivity
# reduction (TBV). Factors are 1 with no vaccine, preserving the equilibrium.
.invlogit <- function(x) 1 / (1 + exp(-x))

# PEV efficacy as a function of days since the last primary dose.
pev_efficacy_curve <- function(profile, t) {
  cs <- exp(profile$cs[1]); rho <- .invlogit(profile$rho[1])
  ds <- exp(profile$ds[1]); dl <- exp(profile$dl[1])
  ab <- cs * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
  profile$vmax * (1 - 1 / (1 + (ab / profile$beta)^profile$alpha))
}

.combine <- function(a, b) 1 - (1 - a) * (1 - b)
.recycle <- function(x, n) if (length(x) == 1L) rep_len(x, n) else x

# Mean PEV efficacy among the vaccinated `tsince` days after primary completion,
# for an arbitrary sequence of boosters. `bspace[j]`/`bcov[j]`/`bprofs[[j]]` are
# booster j's spacing (from primary), conditional coverage (of those who reached
# the previous dose), and profile. The vaccinated are partitioned by the most
# recent booster each has reached; each stratum carries its own decayed efficacy.
pev_protection <- function(prof, bprofs, bspace, bcov, tsince) {
  nb <- 0L                                    # boosters actually reached (contiguous)
  for (j in seq_along(bspace)) {
    if (!is.na(bspace[j]) && bspace[j] <= tsince && !is.na(bcov[j]) && bcov[j] > 0) nb <- j
    else break
  }
  eff0 <- pev_efficacy_curve(prof, tsince)
  if (nb == 0L) return(eff0)
  cumc <- cumprod(bcov[seq_len(nb)])          # fraction reaching booster j
  prot <- (1 - bcov[1]) * eff0                # primary-only stratum
  if (nb >= 2L) for (j in seq_len(nb - 1L))   # reached j, not j+1
    prot <- prot + cumc[j] * (1 - bcov[j + 1L]) *
      pev_efficacy_curve(bprofs[[min(j, length(bprofs))]], tsince - bspace[j])
  prot + cumc[nb] *                            # most recent booster reached
    pev_efficacy_curve(bprofs[[min(nb, length(bprofs))]], tsince - bspace[nb])
}
# booster-coverage row for a cohort vaccinated at `vd` (matrix rows align with `ts`)
.booster_cov_row <- function(bcovm, ts, vd) {
  if (is.null(bcovm)) return(numeric(0))
  bcovm <- as.matrix(bcovm)
  r <- if (length(ts) > 1) { i <- which(ts <= vd); if (length(i)) max(i) else 1L } else 1L
  bcovm[min(r, nrow(bcovm)), ]
}

# Per-(age,time) PEV FOI multiplier [n_age, n_time] (time last for interpolate).
# Models the primary series AND a full booster sequence (EPI + mass).
pev_series <- function(p, age_mid, timesteps) {
  n_age <- length(age_mid)
  if (!isTRUE(p$pev)) return(list(times = c(0, timesteps), vals = matrix(1, n_age, 2)))
  last_dose <- if (length(p$pev_doses)) max(p$pev_doses) else 0
  onsets <- c(p$pev_epi_timesteps, p$mass_pev_timesteps) + last_dose
  grid <- sort(unique(c(seq(0, timesteps, by = 30), timesteps, timesteps + 365,
                        p$pev_epi_timesteps, p$mass_pev_timesteps, onsets,
                        onsets + 7, onsets + 14)))
  grid <- grid[grid >= 0]
  ng <- length(grid)
  red <- matrix(0, n_age, ng)
  or0 <- function(x) if (is.null(x)) 0 else x

  # --- EPI (age-based) primary + boosters ---
  if (!is.null(p$pev_epi_coverages)) {
    prof <- p$pev_profiles[[p$pev_epi_profile_indices[1]]]
    bidx <- p$pev_epi_profile_indices[-1]
    bprofs <- if (length(bidx)) lapply(bidx, function(ix) p$pev_profiles[[ix]]) else list(prof)
    bspace <- p$pev_epi_booster_spacing
    if (length(bspace)) bspace <- pmax(bspace, or0(p$pev_epi_min_wait))   # min_wait floor
    bcovm <- p$pev_epi_booster_coverage
    vax_complete <- p$pev_epi_age + last_dose
    start <- p$pev_epi_timesteps[1]
    for (i in seq_len(n_age)) {
      if (age_mid[i] < vax_complete) next
      tsince <- age_mid[i] - vax_complete
      for (g in seq_len(ng)) if (grid[g] >= start + tsince) {
        vd <- grid[g] - tsince                          # this cohort's vaccination date
        cov <- .cov_at(p$pev_epi_timesteps, p$pev_epi_coverages, vd)
        bcov <- .booster_cov_row(bcovm, p$pev_epi_timesteps, vd)
        red[i, g] <- .combine(red[i, g], cov * pev_protection(prof, bprofs, bspace, bcov, tsince))
      }
    }
  }

  # --- Mass campaigns: primary + boosters, over ALL age bands x ALL campaigns ---
  if (!is.null(p$mass_pev_timesteps)) {
    prof <- p$pev_profiles[[p$mass_pev_profile_indices[1]]]
    bidx <- p$mass_pev_profile_indices[-1]
    bprofs <- if (length(bidx)) lapply(bidx, function(ix) p$pev_profiles[[ix]]) else list(prof)
    bspace <- p$mass_pev_booster_spacing
    if (length(bspace)) bspace <- pmax(bspace, or0(p$mass_pev_min_wait))
    bcovm <- p$mass_pev_booster_coverage
    lo <- p$mass_pev_min_ages; hi <- p$mass_pev_max_ages   # bands applied at EVERY campaign
    nc <- length(p$mass_pev_timesteps)
    cov <- .recycle(p$mass_pev_coverages, nc)
    for (k in seq_len(nc)) {
      tc <- p$mass_pev_timesteps[k]
      bcov <- if (!is.null(bcovm)) as.matrix(bcovm)[min(k, nrow(as.matrix(bcovm))), ] else numeric(0)
      for (b in seq_along(lo)) {
        inband <- which(age_mid >= lo[b] & age_mid < hi[b])
        if (!length(inband)) next
        for (g in seq_len(ng)) if (grid[g] >= tc + last_dose) {
          tsince <- grid[g] - tc - last_dose
          red[inband, g] <- .combine(red[inband, g], cov[k] * pev_protection(prof, bprofs, bspace, bcov, tsince))
        }
      }
    }
  }
  list(times = grid, vals = 1 - red)
}

# TBV per-state infectivity multipliers [n_age, n_time] for U/A/D/Tr. The IBM
# maps transmission-reducing activity (TRA) to state-specific transmission-
# blocking activity (TBA) via per-state mx (tbv_mu/ma/md/mt) and shape tbv_k.
.calculate_TBA <- function(mx, k, tra) {
  offset <- (k / (k + mx))^k
  ((k / (k + mx * (1 - tra)))^k - offset) / (1 - offset)
}
tbv_series <- function(p, age_mid, timesteps) {
  n_age <- length(age_mid)
  ones <- matrix(1, n_age, 2)
  if (!isTRUE(p$tbv))
    return(list(times = c(0, timesteps), fU = ones, fA = ones, fD = ones, fT = ones))
  grid <- sort(unique(c(seq(0, timesteps, by = 30), timesteps, timesteps + 365,
                        p$tbv_timesteps)))
  grid <- grid[grid >= 0]
  ng <- length(grid)
  rU <- rA <- rD <- rT <- matrix(0, n_age, ng)
  ab_of <- function(t) p$tbv_tau * (p$tbv_rho * exp(-t * log(2) / p$tbv_ds) +
                                    (1 - p$tbv_rho) * exp(-t * log(2) / p$tbv_dl))
  tra_of <- function(ab) (ab / p$tbv_tra_mu)^p$tbv_gamma1 /
    ((ab / p$tbv_tra_mu)^p$tbv_gamma1 + p$tbv_gamma2)
  for (k in seq_along(p$tbv_timesteps)) {
    tc <- p$tbv_timesteps[k]; cov <- p$tbv_coverages[k]
    inband <- which(trunc(age_mid / 365) %in% p$tbv_ages)   # exact year set (IBM convention)
    for (g in seq_len(ng)) if (grid[g] >= tc) {
      tra <- tra_of(ab_of(grid[g] - tc))
      rU[inband, g] <- .combine(rU[inband, g], cov * .calculate_TBA(p$tbv_mu, p$tbv_k, tra))
      rA[inband, g] <- .combine(rA[inband, g], cov * .calculate_TBA(p$tbv_ma, p$tbv_k, tra))
      rD[inband, g] <- .combine(rD[inband, g], cov * .calculate_TBA(p$tbv_md, p$tbv_k, tra))
      rT[inband, g] <- .combine(rT[inband, g], cov * .calculate_TBA(p$tbv_mt, p$tbv_k, tra))
    }
  }
  list(times = grid, fU = 1 - rU, fA = 1 - rA, fD = 1 - rD, fT = 1 - rT)
}

# ---- seasonality / flexible carrying capacity ------------------------------
# Per-species larval carrying capacity K(t) = K0 * scaler(t) * rainfall(t)/R_bar,
# where rainfall is the truncated Fourier series (matches malariasimulation).
# Returns [n_spp, n_time] (time last for interpolate). Constant when aseasonal
# and no carrying-capacity schedule (preserving the equilibrium).
.rainfall <- function(t, g0, g, h, floor) {
  r <- g0
  for (i in seq_along(g)) {
    r <- r + g[i] * cos(2 * pi * t * i / 365) + h[i] * sin(2 * pi * t * i / 365)
  }
  pmax(r, floor)
}

carrying_capacity_series <- function(p, K0, timesteps) {
  n_spp <- length(K0)
  seasonal <- isTRUE(p$model_seasonality)
  flexcc <- isTRUE(p$carrying_capacity)
  if (!seasonal && !flexcc) {
    return(list(times = c(0, timesteps), vals = matrix(K0, n_spp, 2)))
  }
  # daily grid under seasonality (matches the IBM's per-day rainfall; avoids
  # clipping the seasonal peak), coarser when only a carrying-capacity schedule.
  by <- if (seasonal) 1 else 5
  grid <- sort(unique(c(seq(0, timesteps, by = by), timesteps, timesteps + 365,
                        if (flexcc) p$carrying_capacity_timesteps)))
  grid <- grid[grid >= 0]
  ng <- length(grid)
  R_bar <- if (seasonal) mean(.rainfall(1:365, p$g0, p$g, p$h, p$rainfall_floor)) else 1
  vals <- matrix(0, n_spp, ng)
  for (gi in seq_len(ng)) {
    t <- grid[gi]
    seas <- if (seasonal) .rainfall(t, p$g0, p$g, p$h, p$rainfall_floor) / R_bar else 1
    for (s in seq_len(n_spp)) {
      scal <- 1
      if (flexcc) {
        idx <- which(p$carrying_capacity_timesteps <= t)
        if (length(idx)) scal <- p$carrying_capacity_scalers[max(idx), s]
      }
      vals[s, gi] <- max(K0[s] * scal * seas, K0[s] * 1e-4)  # keep > 0
    }
  }
  list(times = grid, vals = vals)
}
