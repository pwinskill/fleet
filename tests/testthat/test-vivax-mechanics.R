# The P. vivax block's daily mechanisms, each pinned on its own so that a
# failure names the mechanism rather than only moving a reference value: the
# day's infection and relapse draw, relapse blocked by PEV, the hypnozoite
# ladder and the bite that overwrites its decay, immunity decay, the refractory
# windows and the boost that skips a day's decay, the liver-stage clock after
# radical cure, and the two vivax delay lines (de = 10, delay_gam = 0).

pv_bands <- function(overrides = list(), eir = 10) {
  p <- malariasimulation::get_parameters(utils::modifyList(list(
    human_population = 10000,
    clinical_incidence_rendering_min_ages = c(0, 0),
    clinical_incidence_rendering_max_ages = c(1824, 36499),
    incidence_rendering_min_ages = 0, incidence_rendering_max_ages = 36499), overrides),
    parasite = "vivax")
  eqm(p, eir)
}

# Step a system built from `pars` to day `t` and return its unpacked state.
state_at <- function(pars, t) {
  sys <- dust2::dust_system_create(get_generator(), pars = pars, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  if (t > 0) dust2::dust_system_run_to_time(sys, t)
  dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))
}
# population by hypnozoite level, summed over age, heterogeneity and state
by_level <- function(st, x) {
  d <- c(x$n_age_v, x$n_het_v, x$n_hyp)
  N <- Reduce(`+`, lapply(c("Sv", "Dv", "Av", "Uv", "Trv", "Trv_slow"),
                          function(nm) array(st[[nm]], d)))
  N <- N + apply(array(st$Phv, c(d, x$n_phv)), 1:3, sum)
  apply(N, 3, sum)
}

# Day 1's infections by cell, recomputed by hand from the seed of a build at
# `eir`. Every bite counts, so a person in S is infected with probability
# 1 - exp(-(b EPS + k f)); in A, D and U that hazard competes with the state's
# progression in one draw a day, U's at the immunity-dependent 1 / dpcr. At the
# seed nobody in a cell differs from its mean (K = J^2 / N), so the quadrature
# collapses onto the mean and the whole day can be written out.
seed_infections <- function(x, eir = 10) {
  d <- c(x$n_age, x$n_het, x$n_hyp)
  kk <- c(seq_len(x$n_bat) - 1, rep(0, x$n_hyp - x$n_bat))
  EPS <- eir / 365 * outer(x$psi, x$zeta) / sum(x$het_wt * x$zeta)
  rtot <- array(x$bv * EPS, d) + rep(kk * x$ff, each = x$n_age * x$n_het)
  iS <- 1 - exp(-rtot)
  hv <- rtot
  N <- x$Sv0 + x$Dv0 + x$Av0 + x$Uv0 + x$Trv0
  mA <- ifelse(N > 0, x$JAv0 / N, 0)
  mask <- x$mask20
  IAA20 <- apply(apply(x$JAv0, c(1, 2), sum) * mask, 2, sum) /
    apply(apply(N, c(1, 2), sum) * mask, 2, sum)
  IAM <- x$PM * outer(x$icm_factor, IAA20)
  xA <- mA + array(IAM, d)
  rU <- 1 / (x$dpcr_min + (x$dpcr_max - x$dpcr_min) / (1 + (xA / x$apcr50)^x$kpcr))
  comp <- function(h) (1 - exp(-(hv + h))) * hv / (hv + h)
  inf <- iS * x$Sv0 + comp(-log(1 - x$rA)) * x$Av0 + comp(-log(1 - x$rD)) * x$Dv0 +
    comp(rU) * x$Uv0
  list(inf = inf, rtot = rtot, kk = kk, N = N)
}

test_that("a vivax day infects and relapses at the IBM's probabilities", {
  skip_if_not_installed("malariasimulation")
  p <- pv_bands()
  inp <- build_inputs(p, 10, timesteps = 10)
  x <- inp$pars
  o <- run_simulation_ode(1, p)
  s <- seed_infections(x)
  inf <- s$inf
  hp <- p$human_population
  expect_equal(o$n_infections[1], sum(inf) * hp, tolerance = 1e-10)
  # relapses are the relapse share of each cell's infections
  rel <- inf * rep(s$kk * x$ff, each = x$n_age * x$n_het) / s$rtot
  expect_equal(o$n_relapses[1], sum(rel) * hp, tolerance = 1e-10)
  N <- s$N
  # and the force of infection on mosquitoes is today's infectivity (no lag for
  # vivax), weighted by zeta x psi over the live population
  infv <- x$cD * x$Dv0 + x$cA_v * x$Av0 + x$cU * x$Uv0 + x$cT_vals[1] * x$Trv0
  zp <- outer(x$psi, x$zeta)
  foim <- x$a_vals[1, 1] * sum(zp * apply(infv, c(1, 2), sum)) / sum(zp * apply(N, c(1, 2), sum))
  expect_equal(o$FOIM[1], foim, tolerance = 1e-10)
})

test_that("with no bites every infection is a relapse, and with neither none happen", {
  skip_if_not_installed("malariasimulation")
  p <- pv_bands()
  x <- build_inputs(p, 10, timesteps = 40)$pars
  x$bv <- 0
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  idx <- dust2::dust_unpack_index(sys)
  y <- dust2::dust_system_simulate(sys, 1:30)
  rel <- colSums(y[idx$relapse_g, , drop = FALSE])
  inc <- colSums(y[idx$inc_g, , drop = FALSE])
  expect_equal(rel, inc, tolerance = 1e-12)
  expect_true(all(rel > 0))
  x$ff <- 0
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, 1:30)
  expect_true(all(y[idx$inc_g, ] == 0))
})

test_that("PEV blocks relapses as well as bites", {
  skip_if_not_installed("malariasimulation")
  # calculate_vivax_infections() adds the relapse hazard to the bite hazard and
  # THEN applies the vaccine's efficacy, so a vaccinated carrier relapses less.
  # With bites switched off, only that route can move the relapse count.
  pev <- malariasimulation::set_mass_pev(
    pv_bands(), profile = malariasimulation::rtss_profile, timesteps = 50, coverages = 1,
    min_ages = 0, max_ages = 100 * 365, min_wait = 0, booster_spacing = 365,
    booster_coverage = matrix(0, 1, 1), booster_profile = list(malariasimulation::rtss_booster_profile))
  run_nobite <- function(p) {
    x <- build_inputs(p, 10, timesteps = 200)$pars
    x$bv <- 0
    sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
    dust2::dust_system_set_state_initial(sys)
    idx <- dust2::dust_unpack_index(sys)
    colSums(dust2::dust_system_simulate(sys, 1:200)[idx$relapse_g, , drop = FALSE])
  }
  base <- run_nobite(pv_bands())
  vac <- run_nobite(pev)
  expect_equal(vac[1:100], base[1:100])       # nobody protected until the last dose
  expect_lt(vac[180] / base[180], 0.8)
})

test_that("with no infection the hypnozoite ladder clears exactly as the IBM's does", {
  skip_if_not_installed("malariasimulation")
  # hypnozoite_batch_decay_process: a carrier of k batches loses one a day with
  # probability 1 - exp(-k gammal). With no bites and no relapses nothing else
  # moves a batch count. People die at the age-independent default rate and the
  # dead are replaced by newborns with none -- a reset the IBM queues after the
  # decay, so it overwrites it -- and the population by level follows
  #   N_k' = (1 - mu) (N_k (1 - p_k) + N_{k+1} p_{k+1}) + mu N delta_k0
  # exactly, whatever happens to their disease state.
  p <- pv_bands()
  x <- build_inputs(p, 10, timesteps = 400)$pars
  x$bv <- 0; x$ff <- 0
  mu <- unique(as.vector(x$mu_age_z))
  expect_length(mu, 1)
  pk <- 1 - exp(-x$gammal * (seq_len(x$n_bat) - 1))
  N <- by_level(state_at(x, 0), x)
  for (t in seq_len(365)) {
    N <- (1 - mu) * (N * (1 - pk) + c(N[-1] * pk[-1], 0)) + c(mu * sum(N), rep(0, x$n_bat - 1))
  }
  expect_equal(by_level(state_at(x, 365), x), N, tolerance = 1e-10)
  expect_gt(sum(N[-1]), 0.1)                # a real ladder left to clear
})

test_that("a batch-forming bite overwrites the day's batch decay", {
  skip_if_not_installed("malariasimulation")
  # The decay is queued first and a bite's pmin(k + 1, kmax), taken from the
  # start-of-day k, later, so it wins: the day's batch-formers move up, or stay
  # at the top batch, without losing a batch, and everyone else decays. With no
  # relapse every infection is a bite; with no ageing, death or treatment the
  # population by level after day 1 is, for the start-of-day N_k and the day's
  # bite-infected B_k (all of whom form a batch),
  #   N_k' = N_k - B_k [k < kmax] + B_{k-1} - p_k (N_k - B_k) + p_{k+1} (N_{k+1} - B_{k+1})
  x <- build_inputs(pv_bands(), 10, timesteps = 40)$pars
  x$ff <- 0; x$r_age[] <- 0; x$mu_age_z[] <- 0
  B <- apply(seed_infections(x)$inf, 3, sum)
  N <- by_level(state_at(x, 0), x)
  nb <- x$n_bat
  pk <- 1 - exp(-x$gammal * (seq_len(nb) - 1))
  move <- c(B[-nb], 0)                       # the top batch's bitten stay put
  pre <- N - move + c(0, move[-nb])
  dec <- pk * (N - B)
  expect_equal(by_level(state_at(x, 1), x), pre - dec + c(dec[-1], 0), tolerance = 1e-12)
  # decaying the bitten too, as a separate additive flow would, is measurably different
  expect_gt(max(abs(pk * B)) / max(abs(dec)), 0.01)
})

test_that("vivax immunity decays by exp(-1/d) a day, and its spread by exp(-2/d)", {
  skip_if_not_installed("malariasimulation")
  p <- pv_bands()
  x <- build_inputs(p, 10, timesteps = 40)$pars
  x$bv <- 0; x$ff <- 0
  s0 <- state_at(x, 0); s1 <- state_at(x, 1)
  mu <- unique(as.vector(x$mu_age_z))
  # summed over every cell, ageing moves the stocks and death removes them at
  # the cell mean; with no boost the day's decay is all that is left
  expect_equal(sum(s1$JAv), sum(s0$JAv) * (exp(-1 / x$d_iaa) - mu), tolerance = 1e-12)
  expect_equal(sum(s1$JCv), sum(s0$JCv) * (exp(-1 / x$d_ica) - mu), tolerance = 1e-12)
  expect_equal(sum(s1$KAv), sum(s0$KAv) * (exp(-2 / x$d_iaa) - mu), tolerance = 1e-12)
})

test_that("a refractory window is a stock that withholds the boost, and the boost skips decay", {
  skip_if_not_installed("malariasimulation")
  # boost_immunity() boosts only those whose last boost is at least u days old,
  # and a boost replaces that day's decay. fleet carries the people inside each
  # window as a stock per cell, RAv for IAA and RCv for ICA; a cell's infections
  # boost in the share 1 - R / N.
  x <- build_inputs(pv_bands(), 10, timesteps = 40)$pars
  s <- seed_infections(x)
  inf <- s$inf; N <- s$N
  mu <- unique(as.vector(x$mu_age_z))
  s1 <- state_at(x, 1)
  # the IBM starts everyone's last boost at -1, so nobody is refractory on day 1:
  # every infection boosts ICA, and the boosted skip the decay of what they carry
  e <- exp(-1 / x$d_ica)
  fr <- ifelse(N > 0, inf / N, 0)
  expect_equal(sum(s1$JCv), sum(x$JCv0 * (e - mu) + inf + fr * x$JCv0 * (1 - e)),
               tolerance = 1e-12)
  # decay-then-boost (everyone decays, the boosted then gain 1) is measurably different
  expect_gt(abs(sum(s1$JCv) - sum(x$JCv0 * (e - mu) + inf)) / sum(inf), 1e-3)
  # and each of the day's infected opens a window, in the first stage of both stocks
  RA <- array(s1$RAv, c(dim(N), x$n_ra)); RC <- array(s1$RCv, c(dim(N), x$n_rc))
  expect_equal(sum(RA[, , , 1]), sum(inf), tolerance = 1e-12)
  expect_equal(sum(RC[, , , 1]), sum(inf), tolerance = 1e-12)
  expect_equal(sum(RA[, , , -1]), 0)
  # A population entirely inside IAA's window takes no IAA boost, however many
  # are infected, while ICA, whose window is its own, is boosted as before.
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  idx <- dust2::dust_unpack_index(sys)
  st <- dust2::dust_system_state(sys)
  RA0 <- array(0, c(dim(N), x$n_ra)); RA0[, , , 1] <- N
  st[idx$RAv] <- RA0
  dust2::dust_system_set_state(sys, st)
  dust2::dust_system_run_to_time(sys, 1)
  r1 <- dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))
  expect_equal(sum(r1$JAv), sum(x$JAv0) * (exp(-1 / x$d_iaa) - mu), tolerance = 1e-12)
  expect_equal(sum(r1$JCv), sum(s1$JCv), tolerance = 1e-12)
  # And in between, each cell's infections boost in the share 1 - R / N of that
  # cell, whatever the others hold: a share that differs from cell to cell
  # pins the split cell by cell, not just its two ends.
  f <- array((seq_len(length(N)) %% 5) / 5, dim(N))
  RA0 <- array(0, c(dim(N), x$n_ra)); RA0[, , , 1] <- f * N
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  st <- dust2::dust_system_state(sys)
  st[idx$RAv] <- RA0
  dust2::dust_system_set_state(sys, st)
  dust2::dust_system_run_to_time(sys, 1)
  h1 <- dust2::dust_unpack_state(sys, dust2::dust_system_state(sys))
  e <- exp(-1 / x$d_iaa)
  frA <- ifelse(N > 0, (1 - f) * inf / N, 0)
  expect_equal(sum(h1$JAv), sum(x$JAv0 * (e - mu) + (1 - f) * inf + frA * x$JAv0 * (1 - e)),
               tolerance = 1e-12)
})

test_that("a refractory window holds its people ceiling(u) - 1 days, as a chain of stages", {
  skip_if_not_installed("malariasimulation")
  # n stages, each passed on with probability n / u_eff a day, so someone
  # boosted d days ago is still inside with probability P(Binomial(d, n / u_eff)
  # < n): u_eff days on average, as the IBM's (timestep - last_boosted) >= u
  # blocks the u_eff = ceiling(u) - 1 days after a boost. Nothing else moves
  # here: no infection, ageing, death or batch loss.
  p <- pv_bands()
  x <- build_inputs(p, 10, timesteps = 400)$pars
  expect_equal(x$ua_eff, ceiling(p$ua) - 1)
  expect_equal(x$n_ra, min(x$ua_eff, 4))
  x$bv <- 0; x$ff <- 0; x$gammal <- 0; x$r_age[] <- 0; x$mu_age_z[] <- 0
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  idx <- dust2::dust_unpack_index(sys)
  st <- dust2::dust_system_state(sys)
  N <- x$Sv0 + x$Dv0 + x$Av0 + x$Uv0 + x$Trv0
  RA0 <- array(0, c(dim(N), x$n_ra)); RA0[, , , 1] <- 0.5 * N
  st[idx$RAv] <- RA0
  dust2::dust_system_set_state(sys, st)
  days <- 1:300
  y <- dust2::dust_system_simulate(sys, days)
  held <- colSums(y[idx$RAv, , drop = FALSE]) / (0.5 * sum(N))
  pr <- x$n_ra / x$ua_eff
  expect_equal(held, stats::pbinom(x$n_ra - 1, days, pr), tolerance = 1e-12)
  expect_equal(1 + sum(held), x$ua_eff, tolerance = 1e-6)   # day 0 counts once
  # ICA's window the same way; its few days fit whole stages, so it holds its
  # people for exactly uc_eff days
  sys <- dust2::dust_system_create(get_generator(), pars = x, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  st <- dust2::dust_system_state(sys)
  RC0 <- array(0, c(dim(N), x$n_rc)); RC0[, , , 1] <- 0.5 * N
  st[idx$RCv] <- RC0
  dust2::dust_system_set_state(sys, st)
  y <- dust2::dust_system_simulate(sys, days)
  heldC <- colSums(y[idx$RCv, , drop = FALSE]) / (0.5 * sum(N))
  expect_equal(heldC, stats::pbinom(x$n_rc - 1, days, x$n_rc / x$uc_eff), tolerance = 1e-12)
  expect_equal(1 + sum(heldC), x$uc_eff, tolerance = 1e-9)
})

test_that("radical cure delivers its people to the liver-stage clock, which starts the next day", {
  skip_if_not_installed("malariasimulation")
  p <- malariasimulation::get_parameters(list(human_population = 10000), parasite = "vivax")
  p <- malariasimulation::set_drugs(p, list(malariasimulation::CQ_PQ_params_vivax))
  p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 10, coverages = 0.8)
  x <- build_inputs(eqm(p, 10), 10, timesteps = 40)$pars
  n_ls <- x$n_hyp - x$n_bat
  # primaquine's Weibull(10, 5.5): 4.7 protected days, four stages, each passed
  # on with probability 0.85 a day
  expect_equal(n_ls, 4)
  lv <- function(t) by_level(state_at(x, t), x)
  ls <- x$n_bat + seq_len(n_ls)
  # coverage scheduled for timestep 10 treats day 10's cases, as the IBM reads a
  # setting on its own timestep
  expect_true(all(lv(9)[ls] == 0))
  l10 <- lv(10)
  expect_gt(l10[ls[1]], 0)                  # the day's radically cured arrive
  expect_equal(l10[ls[-1]], rep(0, n_ls - 1))   # and none has left stage 1 yet
  expect_gt(lv(11)[ls[2]], 0)
  # the clock acts on what the day leaves, so no cell goes negative however fast
  # a stage empties beside a Tr stay (0.85 a day on top of Tr's 0.63)
  st <- state_at(x, 40)
  for (nm in c("Sv", "Dv", "Av", "Uv", "Trv", "Trv_slow", "Phv", "JAv", "JCv", "KAv", "KCv",
               "RAv", "RCv"))
    expect_gte(min(st[[nm]]), -1e-15, label = nm)
  # and a refractory stock never holds more people than its cell
  d <- c(x$n_age, x$n_het, x$n_hyp)
  Ncell <- array(st$Sv + st$Dv + st$Av + st$Uv + st$Trv + st$Trv_slow, d) +
    apply(array(st$Phv, c(d, x$n_phv)), 1:3, sum)
  for (nm in c("RAv", "RCv")) {
    R <- apply(array(st[[nm]], c(d, dim(st[[nm]])[4])), 1:3, sum)
    expect_lte(max(R - Ncell), 1e-15, label = nm)
  }
})

test_that("the vivax delay lines are the IBM's: EIR lagged de = 10 days, infectivity not at all", {
  skip_if_not_installed("malariasimulation")
  day <- 100
  nets <- malariasimulation::set_bednets(pv_bands(), timesteps = day, coverages = 0.8,
    retention = 5 * 365, dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
    rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)
  tbv <- malariasimulation::set_tbv(pv_bands(), timesteps = day, coverages = 0.9, ages = 5:40)
  base <- run_simulation_ode(140, pv_bands())
  o <- run_simulation_ode(140, nets)
  first <- function(a, b) which(a != b)[1]
  expect_equal(first(o$FOIM, base$FOIM), day + 1)
  expect_equal(first(o$EIR, base$EIR), day + 1 + 10)
  expect_equal(first(o$n_infections, base$n_infections), day + 1 + 10)
  # TBV cuts infectivity from day 101, and with delay_gam = 0 mosquitoes feel it
  # the same day; the fewer infected mosquitoes leave incubation ceiling(dem) - 1
  # = 9 days later, and their bites reach humans de = 10 days after that
  o <- run_simulation_ode(140, tbv)
  expect_equal(first(o$FOIM, base$FOIM), day + 1)
  expect_equal(first(o$EIR, base$EIR), day + 1 + 9 + 10)
  expect_equal(first(o$n_infections, base$n_infections), day + 1 + 9 + 10)
})
