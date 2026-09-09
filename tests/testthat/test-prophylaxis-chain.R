# Prophylaxis Erlang chains, chemoprevention pulse renewal and the custom-demography
# mosquito-sizing convention. (Drafted in the Sept 2026 parallel review.)

test_that("erlang_stages moment-matches the Weibull and is capped at 20", {
  cv2 <- function(s) gamma(1 + 2 / s) / gamma(1 + 1 / s)^2 - 1
  # chemoprevention chain (no Tr stage in front): k = 1/CV^2 (real scales, so the
  # per-stage rate cap of 4/day does not bind)
  expect_equal(erlang_stages(4.3, 38.1), 14L)                 # SP-AQ
  expect_equal(erlang_stages(4.4, 28.1), 15L)                 # DHA-PQP
  expect_equal(erlang_stages(11.3, 10.6), 20L)                # AL: 1/CV^2 = 87 -> cap
  expect_equal(erlang_stages(1, 30), 1L)                      # exponential Weibull -> one stage
  expect_equal(erlang_stages(0.5, 30), 1L)                    # floor at 1
  expect_equal(erlang_stages(2.1, 30), as.integer(round(1 / cv2(2.1))))
  # a one-day protection cannot take 14 stages at 14/day: the rate cap binds
  expect_equal(erlang_stages(4.3, 1), 3L)
})

## Split from the block above so the closed-form checks there keep running when
## malariasimulation (a Suggests-only GitHub Remote) is absent: everything from here
## down reads its Weibull parameters out of the IBM's own drug tables.
test_that("the post-treatment chain matches the variance of the whole Tr + chain sojourn", {
  skip_if_not_installed("malariasimulation")
  # a protection less variable than Tr itself (AL: sd 1.1 d vs 5.5 d) gets ONE stage
  rT <- 1 - exp(-1 / 5); tr <- 1 / rT
  mc <- function(pr) .chain_mean_after_tr(pr[3], pr[4], rT)
  al <- malariasimulation::AL_params; sp <- malariasimulation::SP_AQ_params
  expect_equal(erlang_stages(al[3], al[4], tr = tr, m_chain = mc(al)), 1L)
  expect_equal(erlang_stages(sp[3], sp[4], tr = tr, m_chain = mc(sp)), 16L)
  # the chain mean is the integrated protection left after Tr: ~ mean - 1/rT when
  # W ~ 1 throughout Tr (SP-AQ), well above it when W decays during Tr (AL)
  expect_equal(mc(sp), .weibull_mean(sp[3], sp[4]) - tr, tolerance = 2e-3)
  expect_gt(mc(al), .weibull_mean(al[3], al[4]) - tr + 0.5)
  # a 50/50 AL + SP-AQ mixture is bimodal: far fewer stages than either drug alone
  expect_lt(erlang_stages(c(al[3], sp[3]), c(al[4], sp[4]), c(0.5, 0.5)), 5L)
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
  expect_equal(drug_mix(p1, eqp)$n_ph, 16L)                   # SP-AQ first line
  p2 <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.5)
  expect_equal(drug_mix(p2, eqp)$n_ph, 1L)                    # AL first line
  expect_equal(chemoprevention_prophylaxis(p1, eqp)$n_stages, 1L)   # no chemoprevention
  ps <- malariasimulation::set_smc(p, drug = 2, timesteps = 100, coverages = 0.8,
                                   min_ages = 91, max_ages = 1825)
  chp <- chemoprevention_prophylaxis(ps, eqp)
  expect_equal(chp$n_stages, 14L)
  expect_equal(chp$rate, 1 / (38.1 * gamma(1 + 1 / 4.3)))
  inp <- build_inputs(ps, 20)                                 # counts reach odin pars + meta
  expect_equal(inp$pars$n_phc, 14L); expect_equal(inp$meta$n_phc, 14L)
  expect_equal(inp$pars$n_ph, 1L)
  expect_equal(dim(inp$pars$Phc0), c(inp$meta$n_age, inp$meta$n_het, 14L))
  expect_equal(dim(inp$pars$Ph0), c(inp$meta$n_age, inp$meta$n_het, 1L))
  expect_equal(build_inputs(ps, 20, n_phc = 1)$pars$n_phc, 1L)
})

test_that("the Ph_c chain follows the Weibull protection curve; n_phc = 1 is the old exponential", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  # one MDA round over ALL ages at day 10 at coverage 1: 0.9 (efficacy) of everyone
  # is moved into Ph_c at t = 10+, and there is no clinical treatment, so Ph_count
  # afterwards is exactly the protected cohort.
  p <- malariasimulation::set_mda(p, drug = 1, timesteps = 10, coverages = 1,
                                  min_ages = 0, max_ages = 200 * 365)
  hp <- p$human_population; eta <- 1 / p$average_age
  m <- 38.1 * gamma(1 + 1 / 4.3)                             # Weibull mean, 34.68 d
  W <- function(s) exp(-(s / 38.1)^4.3) * exp(-eta * s)      # IBM protection x survival
  E <- function(s) exp(-s / m) * exp(-eta * s)               # old single compartment
  surv <- function(o, s) o$Ph_count[match(10 + s, o$timestep)] / (0.9 * hp)
  s <- 1:90
  o14 <- run_simulation_ode(100, eqm(p, 20))           # n_phc auto = 14
  o1  <- run_simulation_ode(100, eqm(p, 20), tuning = list(n_phc = 1))
  expect_equal(o14$Ph_count[o14$timestep == 10], 0)          # row 10 is pre-pulse
  expect_equal(surv(o14, 1), 1 - pgamma(1, 14, rate = 14 / m), tolerance = 1e-3)
  # the chain tracks the Weibull to within a few points everywhere; the exponential
  # is 28 points low at day 30 (0.42 vs 0.70) and lingers far too long at day 60
  expect_lt(max(abs(surv(o14, s) - W(s))), 0.06)
  expect_equal(surv(o14, 30), 0.670 * exp(-30 * eta), tolerance = 0.01)
  expect_gt(max(abs(surv(o1, s) - W(s))), 0.3)
  expect_equal(surv(o1, 30), E(30), tolerance = 1e-3)
  expect_lt(surv(o14, 60), 0.02); expect_gt(surv(o1, 60), 0.15)
  # n_phc = 1 decays as one exponential at rP_c (+ deaths): log-linear
  lp <- log(o1$Ph_count[o1$timestep %in% (15:70)])
  slope <- coef(lm(lp ~ seq_along(lp)))[[2]]
  expect_equal(slope, -(1 / m + eta), tolerance = 1e-3)
})

test_that("a chemoprevention pulse renews Ph/Ph_c protection and clears Tr_slow", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params,
                                              malariasimulation::AL_params))
  # SP-AQ first line so the treatment chain has several stages (16), plus slow
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
  sys <- dust2::dust_system_create(malaria_ode, pars = inp$pars, n_particles = 1)
  dust2::dust_system_set_state_initial(sys)
  uidx <- dust2::dust_unpack_index(sys)
  ev <- chemoprevention_events(p, 100)
  get <- function() {
    sv <- dust2::dust_system_state(sys)
    list(S = sv[uidx$S], D = sv[uidx$D], A = sv[uidx$A], U = sv[uidx$U],
         Tr = sv[uidx$Tr], Trs = sv[uidx$Tr_slow],
         Ph = array(sv[uidx$Ph], c(n_age, n_het, k)), Phc = array(sv[uidx$Ph_c], c(n_age, n_het, kc)))
  }
  tot <- function(x) matrix(x$S + x$D + x$A + x$U + x$Tr + x$Trs, n_age, n_het) +
    apply(x$Ph, c(1, 2), sum) + apply(x$Phc, c(1, 2), sum)
  # first round at day 20, then 30 days so Ph_c mass has spread down the chain
  dust2::dust_system_run_to_time(sys, 20); apply_chemoprevention_pulse(sys, uidx, inp$meta, ev[[1]])
  dust2::dust_system_run_to_time(sys, 50)
  before <- get()
  expect_gt(sum(before$Trs), 0); expect_gt(sum(before$Ph), 0); expect_gt(sum(before$Phc[, , -1]), 0)
  apply_chemoprevention_pulse(sys, uidx, inp$meta, ev[[2]])
  after <- get()
  fr <- ev[[2]]$frac
  expect_equal(fr, 0.6 * 0.9)                                 # coverage x SP-AQ efficacy (ETF = 0)
  # the covered fraction leaves EVERY state -- Tr_slow and every prophylaxis stage included ...
  for (nm in c("S", "D", "A", "U", "Tr", "Trs", "Ph"))
    expect_equal(after[[nm]], (1 - fr) * before[[nm]], tolerance = 1e-12)
  expect_equal(after$Phc[, , -1], (1 - fr) * before$Phc[, , -1], tolerance = 1e-12)
  # ... and re-enters stage 1 of Ph_c with a fresh protection clock
  expect_equal(after$Phc[, , 1], (1 - fr) * before$Phc[, , 1] + fr * tot(before), tolerance = 1e-12)
  expect_equal(sum(tot(after)), sum(tot(before)), tolerance = 1e-12)   # population conserved
})

test_that("n_ph = 1 is the previous single-compartment seed; the chain seed conserves mass", {
  skip_if_not_installed("malariasimulation")
  old_sdb <- function(FOI, phi, prop, r, eta, rA, rD, rU, rT, rP, ft) {   # pre-chain version
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
  new14 <- do.call(solve_disease_block, c(args, n_ph = 14L))
  expect_equal(dim(new14$P), c(nrow(eqn), 14L))
  expect_equal(with(new14, S + T + D + A + U + rowSums(P)), unname(eqn[, "prop"]), tolerance = 1e-12)
  expect_true(all(new14$P > 0))
  # every stage of the seed satisfies the chain's own steady-state balance
  # (stage m>1 in the first age group: n*rP*P[m-1] = (n*rP + r1 + eta) P[m])
  rPk <- 14 / 29; r1 <- eqn[1, "r"]
  expect_equal(rPk * new14$P[1, 1:13], (rPk + r1 + eqp$eta) * new14$P[1, 2:14], tolerance = 1e-12)
})

test_that("the chain seed adds no drift under treatment (Ph_count and EIR)", {
  skip_if_not_installed("malariasimulation")
  gp <- malariasimulation::get_parameters
  p <- malariasimulation::set_drugs(gp(), list(malariasimulation::SP_AQ_params))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.5)
  p$bite_dedup <- 0; p$acquired_immunity_offset <- 0
  o16 <- run_simulation_ode(730, eqm(p, 20))          # n_ph auto = 16
  o1  <- run_simulation_ode(730, eqm(p, 20), tuning = list(n_ph = 1))
  drift <- function(o, col) max(abs(o[[col]] - o[[col]][1])) / o[[col]][1]
  # the seed under treatment relaxes by a few percent in Ph_count and ~0.5% in EIR
  # (pre-existing: the immunity-boost seed, identical at n_ph = 1); the chain must
  # not add to it, and it must carry the same equilibrium prophylaxis mass (the
  # per-stage aging split moves it by ~0.2%)
  expect_lt(drift(o16, "Ph_count"), 0.05)
  expect_lt(drift(o16, "EIR"), 0.01)
  expect_lt(drift(o16, "Ph_count"), 1.1 * drift(o1, "Ph_count") + 1e-6)
  expect_lt(drift(o16, "EIR"), 1.1 * drift(o1, "EIR") + 1e-6)
  expect_equal(o16$Ph_count[1], o1$Ph_count[1], tolerance = 5e-3)
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

test_that("chain stage arguments are validated", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters()
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = 0)), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = 2.5)), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_ph = "a")), "n_ph")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_phc = c(1, 2))), "n_phc")
  expect_error(run_simulation_ode(10, eqm(p, 20), tuning = list(n_eip = 0)), "n_eip")
})
