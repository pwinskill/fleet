# Prophylaxis chains, chemoprevention pulse renewal and the custom-demography
# mosquito-sizing convention.

## The protected days of a dose, computed here independently of the package: the
## IBM's protection on day n after the dose is W(n) = exp(-(n/scale)^shape), so
## the whole days protected T have P(T >= n) = W(n), mean sum(W) and variance
## sum((2n - 1) W) - mean^2.
days_moments <- function(shape, scale, w = 1) {
  n <- 1:2000; w <- w / sum(w)
  W <- Reduce(`+`, Map(function(sh, sc, wi) wi * exp(-(n / sc)^sh), shape, scale, w))
  list(n = n, W = W, M = sum(W), V = sum((2 * n - 1) * W) - sum(W)^2)
}

test_that("erlang_stages matches the protection's mean and variance on the daily clock", {
  # k geometric stages at p = k / M have mean M and variance M^2 / k - M, so the
  # matching count is M^2 / (V + M), held to k <= M (a stage lasts a day at least)
  k_of <- function(shape, scale) {
    d <- days_moments(shape, scale)
    erlang_stages(d$M, d$V)
  }
  expect_equal(k_of(4.3, 38.1), 10L)                  # SP-AQ
  expect_equal(k_of(4.4, 28.1), 9L)                   # DHA-PQP
  expect_equal(k_of(11.3, 10.6), 9L)                  # AL: nearly fixed ~9.6 days -> k <= M binds
  expect_equal(k_of(1, 30), 1L)                       # exponential Weibull -> one stage
  expect_equal(k_of(0.5, 30), 1L)                     # floor at 1
  d <- days_moments(2.1, 30)
  expect_equal(k_of(2.1, 30), as.integer(round(d$M^2 / (d$V + d$M))))
  # a protection of under a day cannot be split at all
  expect_equal(erlang_stages(0.37, 0.1), 1L)
  # the package's own day sums agree with the independent ones
  pd <- .protection_days(4.3, 38.1)
  expect_equal(pd$mean, days_moments(4.3, 38.1)$M, tolerance = 1e-12)
  expect_equal(pd$var, days_moments(4.3, 38.1)$V, tolerance = 1e-10)
})

## Split from the block above so the closed-form checks there keep running when
## malariasimulation (a Suggests-only GitHub Remote) is absent: everything from here
## down reads its Weibull parameters out of the IBM's own drug tables.
test_that("the post-treatment chain matches the variance of the whole Tr + chain sojourn", {
  skip_if_not_installed("malariasimulation")
  rT <- 1 - exp(-1 / 5)
  var_tr <- (1 - rT) / rT^2                             # Tr is left with probability rT a day
  al <- malariasimulation::AL_params; sp <- malariasimulation::SP_AQ_params
  mc <- function(pr) .chain_mean_after_tr(.protection_days(pr[3], pr[4]), rT)
  k <- function(pr) {
    d <- days_moments(pr[3], pr[4])
    erlang_stages(d$M, d$V, var_tr = var_tr, m = mc(pr))
  }
  # a protection less variable than Tr itself (AL: sd 1.1 d vs 5.0 d) gets ONE stage
  expect_equal(k(al), 1L)
  expect_equal(k(sp), 9L)
  # the chain carries the protected days left once out of Tr: nearly all of
  # SP-AQ's (the whole less ~1/rT, W being ~1 throughout Tr), less of AL's
  d_sp <- days_moments(sp[3], sp[4]); d_al <- days_moments(al[3], al[4])
  expect_equal(mc(sp), d_sp$M - 1 / rT, tolerance = 2e-3)
  expect_gt(mc(al), d_al$M - 1 / rT + 0.5)
  expect_equal(mc(sp), sum((1 - (1 - rT)^(d_sp$n - 1)) * d_sp$W), tolerance = 1e-12)
  # a 50/50 AL + SP-AQ mixture is bimodal: far fewer stages than either drug alone
  mix <- days_moments(c(al[3], sp[3]), c(al[4], sp[4]), c(0.5, 0.5))
  expect_lt(erlang_stages(mix$M, mix$V), 5L)
})

test_that("drug_mix / chemoprevention_prophylaxis derive chain length and rate from the drug", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::AL_params,
                                              malariasimulation::SP_AQ_params))
  eqp <- translate_parameters(p)
  for (nm in c("rA", "rD", "rU", "rT")) eqp[[nm]] <- 1 - exp(-eqp[[nm]])   # as build_inputs()
  expect_equal(drug_mix(p, eqp)$n_ph, 1L)                     # no clinical treatment -> 1 stage
  p1 <- malariasimulation::set_clinical_treatment(p, drug = 2, timesteps = 1, coverages = 0.5)
  expect_equal(drug_mix(p1, eqp)$n_ph, 9L)                    # SP-AQ first line
  p2 <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.5)
  expect_equal(drug_mix(p2, eqp)$n_ph, 1L)                    # AL first line
  expect_equal(chemoprevention_prophylaxis(p1, eqp)$n_stages, 1L)   # no chemoprevention
  ps <- malariasimulation::set_smc(p, drug = 2, timesteps = 100, coverages = 0.8,
                                   min_ages = 91, max_ages = 1825)
  chp <- chemoprevention_prophylaxis(ps, eqp)
  expect_equal(chp$n_stages, 10L)
  expect_equal(chp$rate, 1 / days_moments(4.3, 38.1)$M, tolerance = 1e-12)
  inp <- build_inputs(ps, 20)                                 # counts reach odin pars + meta
  expect_equal(inp$pars$n_phc, 10L); expect_equal(inp$meta$n_phc, 10L)
  expect_equal(inp$pars$n_ph, 1L)
  expect_equal(dim(inp$pars$Phc0), c(inp$meta$n_age, inp$meta$n_het, 10L))
  expect_equal(dim(inp$pars$Ph0), c(inp$meta$n_age, inp$meta$n_het, 1L))
  expect_equal(build_inputs(ps, 20, n_phc = 1)$pars$n_phc, 1L)
  # a count whose stages, with the day's ageing and deaths, would lose more than
  # they hold is refused
  expect_error(build_inputs(ps, 20, n_phc = 40), "more than it holds")
})

test_that("the default post-treatment count fits the shortest protection on the schedule", {
  skip_if_not_installed("malariasimulation")
  # SP-AQ first line at the seed (9 stages), switching to AL (a ~5-day chain) later:
  # 9 stages would each have to be left more than once a day after the switch
  p <- malariasimulation::set_drugs(malariasimulation::get_parameters(),
                                    list(malariasimulation::SP_AQ_params, malariasimulation::AL_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = c(1, 400), coverages = c(0.5, 0))
  p <- malariasimulation::set_clinical_treatment(p, drug = 2, timesteps = c(1, 400), coverages = c(0, 0.5))
  inp <- build_inputs(p, 20, timesteps = 800)
  expect_lte(inp$pars$n_ph * max(inp$pars$rP_vals), 1)
  o <- run_simulation_ode(800, eqm(p, 20))
  expect_true(all(is.finite(o$Ph_count)) && all(o$Ph_count >= 0))
})

test_that("the Ph_c chain follows the Weibull protection curve; n_phc = 1 is a single stage", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  # one MDA round over ALL ages at day 10 at coverage 1: 0.9 (efficacy) of everyone
  # enters Ph_c at the end of day 10, and there is no clinical treatment, so
  # Ph_count afterwards is exactly the protected cohort. Row 10 + s is the start
  # of day s after the dose.
  p <- malariasimulation::set_mda(p, drug = 1, timesteps = 10, coverages = 1,
                                  min_ages = 0, max_ages = 200 * 365)
  # at a vanishing EIR hardly anyone is LM-detectable, so the round's Tr phase
  # (the detectable pass through Tr before protection) carries nobody and the
  # whole covered population enters Ph_c
  hp <- p$human_population; eta <- 1 / p$average_age
  d <- days_moments(4.3, 38.1)
  W <- function(s) exp(-(s / 38.1)^4.3) * (1 - eta)^(s - 1)   # IBM protection x survival
  surv <- function(o, s) o$Ph_count[match(10 + s, o$timestep)] / (0.9 * hp)
  s <- 1:90
  o10 <- run_simulation_ode(100, eqm(p, 1e-3))         # n_phc auto = 10
  o1  <- run_simulation_ode(100, eqm(p, 1e-3), tuning = list(n_phc = 1))
  expect_equal(o10$Ph_count[o10$timestep == 10], 0)          # row 10 is pre-pulse
  expect_equal(surv(o10, 1), 1, tolerance = 1e-3)            # the first protected day
  # the chain is its closed form: 10 stages each left with probability 10 / M a
  # day, a sum of geometric sojourns
  nb <- function(n) 1 - stats::pnbinom(n - 10 - 1, 10, 10 / d$M)
  expect_lt(max(abs(surv(o10, s) - nb(s) * (1 - eta)^(s - 1))), 2e-3)
  # ... which tracks the Weibull to within a few points everywhere; one stage is
  # 28 points low at day 30 (0.42 vs 0.70) and lingers far too long at day 60
  expect_lt(max(abs(surv(o10, s) - W(s))), 0.06)
  expect_equal(surv(o10, 30), 0.666, tolerance = 0.01)
  expect_gt(max(abs(surv(o1, s) - W(s))), 0.3)
  expect_equal(surv(o1, 30), (1 - 1 / d$M)^29 * (1 - eta)^29, tolerance = 1e-3)
  expect_lt(surv(o10, 60), 0.02); expect_gt(surv(o1, 60), 0.15)
  # n_phc = 1 decays as one geometric stage (+ deaths): log-linear, exactly
  lp <- log(o1$Ph_count[o1$timestep %in% (15:70)])
  slope <- coef(lm(lp ~ seq_along(lp)))[[2]]
  expect_equal(slope, log(1 - 1 / d$M - eta), tolerance = 1e-3)
})

test_that("a chemoprevention pulse treats as the IBM does and renews protection", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params,
                                              malariasimulation::AL_params))
  # SP-AQ first line so the treatment chain has several stages, plus slow
  # parasite clearance from t = 0 so Tr_slow carries mass
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.5)
  p <- malariasimulation::set_antimalarial_resistance(p, drug = 1, timesteps = 0,
    artemisinin_resistance_proportion = 0.5, partner_drug_resistance_proportion = 0,
    slow_parasite_clearance_probability = 0.5, early_treatment_failure_probability = 0,
    late_clinical_failure_probability = 0, late_parasitological_failure_probability = 0,
    reinfection_during_prophylaxis_probability = 0, slow_parasite_clearance_time = 10)
  p <- malariasimulation::set_smc(p, drug = 1, timesteps = c(20, 50), coverages = c(0.6, 0.6),
                                  min_ages = c(0, 0), max_ages = rep(200 * 365, 2))
  inp <- build_inputs(p, 20, timesteps = 100)
  n_age <- inp$meta$n_age; n_het <- inp$meta$n_het; k <- inp$meta$n_ph; kc <- inp$meta$n_phc
  expect_gt(k, 1L)
  sys <- dust2::dust_system_create(malaria_daily, pars = inp$pars, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  uidx <- dust2::dust_unpack_index(sys)
  ev <- chemoprevention_events(p, 100)
  kt <- inp$meta$n_phct
  get <- function() {
    sv <- dust2::dust_system_state(sys)
    m <- function(nm) matrix(sv[uidx[[nm]]], n_age, n_het)
    list(S = m("S"), D = m("D"), A = m("A"), U = m("U"), Tr = m("Tr"), Trs = m("Tr_slow"),
         Trc = m("Tr_c"), Trcs = m("Tr_cs"), Jc = m("J_c"), Jcs = m("J_cs"), ID = m("ID"),
         Ph = array(sv[uidx$Ph], c(n_age, n_het, k)), Phc = array(sv[uidx$Ph_c], c(n_age, n_het, kc)),
         Phct = array(sv[uidx$Ph_ct], c(n_age, n_het, kt)))
  }
  tot <- function(x) x$S + x$D + x$A + x$U + x$Tr + x$Trs + x$Trc + x$Trcs +
    apply(x$Ph, c(1, 2), sum) + apply(x$Phc, c(1, 2), sum) + apply(x$Phct, c(1, 2), sum)
  # first round at day 20, then 30 days so Ph_c mass has spread down the chain
  dust2::dust_system_run_to_time(sys, 20); apply_chemoprevention_pulse(sys, uidx, inp$meta, ev[[1]])
  dust2::dust_system_run_to_time(sys, 50)
  before <- get()
  expect_gt(sum(before$Trs), 0); expect_gt(sum(before$Ph), 0); expect_gt(sum(before$Phc[, , -1]), 0)
  apply_chemoprevention_pulse(sys, uidx, inp$meta, ev[[2]])
  after <- get()
  fr <- ev[[2]]$frac
  expect_equal(fr, 0.6 * 0.9)                                 # coverage x SP-AQ efficacy (ETF = 0)
  expect_equal(ev[[2]]$spc, 0.5 * 0.5)                        # resistance: slow clearance
  # the covered fraction leaves EVERY state -- Tr_slow and every prophylaxis stage included
  for (nm in c("S", "D", "A", "U", "Tr", "Trs", "Ph"))
    expect_equal(after[[nm]], (1 - fr) * before[[nm]], tolerance = 1e-12)
  expect_equal(after$Phc[, , -1], (1 - fr) * before$Phc[, , -1], tolerance = 1e-12)
  # the clinical and the LM-detectable asymptomatic it treats go to Tr first,
  # infectious at their infectivity x drug_rel_c, the slow-clearing share to Tr_cs
  dp <- inp$meta$detect
  q <- dp$d1 + (1 - dp$d1) / (1 + (before$ID / dp$id0)^dp$kd * dp$fd)
  cA <- dp$cU + (dp$cD - dp$cU) * q^dp$g_inf
  to_tr <- fr * (before$D + q * before$A)
  j_in <- fr * (dp$cD * before$D + cA * q * before$A) * p$drug_rel_c[1]
  expect_equal(after$Trc, (1 - fr) * before$Trc + 0.75 * to_tr, tolerance = 1e-12)
  expect_equal(after$Trcs, (1 - fr) * before$Trcs + 0.25 * to_tr, tolerance = 1e-12)
  expect_equal(after$Jc, (1 - fr) * before$Jc + 0.75 * j_in, tolerance = 1e-12)
  # everyone else treated re-enters stage 1 of Ph_c with a fresh protection clock
  expect_equal(after$Phc[, , 1], (1 - fr) * before$Phc[, , 1] + fr * tot(before) - to_tr,
               tolerance = 1e-12)
  expect_equal(sum(tot(after)), sum(tot(before)), tolerance = 1e-12)   # population conserved
})

test_that("n_ph = 1 is the single-compartment seed; the chain seed conserves mass", {
  skip_if_not_installed("malariasimulation")
  old_sdb <- function(FOI, phi, prop, r, eta, rA, rD, rU, rT, rP, ft) {   # single-stage version
    n <- length(FOI); S <- Tc <- D <- A <- U <- P <- numeric(n)
    for (i in seq_len(n)) {
      re <- r[i] + eta; betaT <- rT + re; betaD <- rD + re
      betaA <- FOI[i] * phi[i] + rA + re; betaU <- FOI[i] + rU + re; betaP <- rP + re
      aT <- ft * phi[i] * FOI[i] / betaT; aP <- rT * aT / betaP
      aD <- (1 - ft) * phi[i] * FOI[i] / betaD
      if (i == 1) { bT <- 0; bD <- 0; bP <- 0 } else {
        bT <- r[i - 1] * Tc[i - 1] / betaT; bD <- r[i - 1] * D[i - 1] / betaD
        bP <- (rT * bT + r[i - 1] * P[i - 1]) / betaP
      }
      Y <- (prop[i] - (bT + bD + bP)) / (1 + aT + aD + aP)
      Tc[i] <- aT * Y + bT; D[i] <- aD * Y + bD; P[i] <- aP * Y + bP
      rA_in <- if (i == 1) 0 else r[i - 1] * A[i - 1]; rU_in <- if (i == 1) 0 else r[i - 1] * U[i - 1]
      A[i] <- (rA_in + (1 - phi[i]) * Y * FOI[i] + rD * D[i]) / (betaA + (1 - phi[i]) * FOI[i])
      U[i] <- (rU_in + rA * A[i]) / betaU; S[i] <- Y - A[i] - U[i]
    }
    list(S = S, T = Tc, D = D, A = A, U = U, P = P)
  }
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  eqp <- translate_parameters(p)
  for (nm in c("rA", "rD", "rU", "rT")) eqp[[nm]] <- 1 - exp(-eqp[[nm]])
  eqn <- malariaEquilibrium::human_equilibrium_no_het(EIR = 20, ft = 0.45, p = eqp, age = default_age_lower())
  args <- list(eqn[, "FOI"], eqn[, "phi"], eqn[, "prop"], eqn[, "r"], eqp$eta,
               eqp$rA, eqp$rD, eqp$rU, eqp$rT, 1 / 29, 0.45)
  old <- do.call(old_sdb, args)
  new1 <- do.call(solve_disease_block, c(args, n_ph = 1L))
  for (nm in c("S", "T", "D", "A", "U")) expect_equal(new1[[nm]], old[[nm]], tolerance = 1e-12)
  expect_equal(drop(new1$P), old$P, tolerance = 1e-12)
  new9 <- do.call(solve_disease_block, c(args, n_ph = 9L))
  expect_equal(dim(new9$P), c(nrow(eqn), 9L))
  expect_equal(with(new9, S + T + D + A + U + rowSums(P)), unname(eqn[, "prop"]), tolerance = 1e-12)
  expect_true(all(new9$P > 0))
  # every stage of the seed satisfies the chain's own balance, which is the same
  # on the daily clock as in continuous time (stage m>1 in the first age group:
  # n*rP*P[m-1] = (n*rP + r1 + eta) P[m])
  rPk <- 9 / 29; r1 <- eqn[1, "r"]
  expect_equal(rPk * new9$P[1, 1:8], (rPk + r1 + eqp$eta) * new9$P[1, 2:9], tolerance = 1e-12)
})

test_that("the chain seed adds no drift under treatment (Ph_count and EIR)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  # at efficacy 1, so the raw coverage the seed takes is also the coverage the
  # run clears (see test-equilibrium.R)
  spaq <- malariasimulation::SP_AQ_params; spaq[1] <- 1
  p <- malariasimulation::set_drugs(gp(), list(spaq))
  # treatment in force at timestep 0, so the humans are seeded treated (from
  # timestep 1 they start untreated, as the IBM's do)
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 0.5)
  p$bite_dedup <- 0; p$acquired_immunity_offset <- 0
  ok <- run_simulation_ode(730, eqm(p, 20))           # n_ph auto = 9
  o1 <- run_simulation_ode(730, eqm(p, 20), tuning = list(n_ph = 1))
  drift <- function(o, col) max(abs(o[[col]] - o[[col]][1])) / o[[col]][1]
  # the seed under treatment relaxes by a few percent in Ph_count and under 1% in
  # EIR (the seed is not the model's own fixed point, identically at n_ph = 1);
  # the chain must not add to it, and it must carry the same equilibrium
  # prophylaxis mass (the per-stage aging split moves it by ~0.2%)
  expect_lt(drift(ok, "Ph_count"), 0.05)
  expect_lt(drift(ok, "EIR"), 0.01)
  expect_lt(drift(ok, "Ph_count"), 1.1 * drift(o1, "Ph_count") + 1e-6)
  expect_lt(drift(ok, "EIR"), 1.1 * drift(o1, "EIR") + 1e-6)
  expect_equal(ok$Ph_count[1], o1$Ph_count[1], tolerance = 5e-3)
})

test_that("hold_init_EIR = TRUE realises init_EIR under custom demography", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  dr <- c(0.048, 0.007, 0.003, 0.004, 0.008, 0.020, 0.050, 0.120) / 365
  p <- malariasimulation::set_demography(gp(), agegroups = round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365),
                                         timesteps = 0, deathrates = matrix(dr, nrow = 1))
  p$hold_init_EIR <- TRUE; p$bite_dedup <- 0; p$acquired_immunity_offset <- 0
  inp <- build_inputs(p, 20)
  expect_equal(inp$meta$eir_seed, 20); expect_true(is.na(inp$meta$total_M_ibm))
  o <- run_simulation_ode(400, eqm(p, 20))
  expect_equal(o$EIR[1], 20, tolerance = 1e-6)
  expect_lt(max(abs(o$EIR - 20)) / 20, 1e-2)
  p0 <- p; p0$hold_init_EIR <- NULL                       # default: the IBM's (smaller) density
  expect_gt(inp$meta$total_M / build_inputs(p0, 20)$meta$total_M, 1.2)
  p$hold_init_EIR <- "yes"
  expect_error(build_inputs(p, 20), "hold_init_EIR")
})

test_that("chain stage and sub-step arguments are validated", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = 0)), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = 2.5)), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = "a")), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_phc = c(1, 2))), "n_phc")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_sub = 0)), "n_sub")
})
