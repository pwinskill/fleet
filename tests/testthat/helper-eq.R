## The target EIR reaches fleet on the parameter list, exactly as it reaches the
## IBM -- run_simulation_ode() has no init_EIR argument. This alias keeps a test
## call on one line; it is nothing but set_equilibrium(), and tests that care
## about the seeding convention itself should spell that out instead.
eqm <- function(p, eir) malariasimulation::set_equilibrium(p, init_EIR = eir)

## LM prevalence in a band, as postie computes it: the detected count over the
## band's population. (p_detect_lm_* is a count too, malariasimulation's expected
## number detected, not a proportion.)
pfpr <- function(o, tag = "730_3650") {
  o[[paste0("n_detect_lm_", tag)]] / o[[paste0("n_age_", tag)]]
}

## The population, summed over the infection states. A drug-protected person is
## uninfected and counted in S_count, as the IBM counts them; Ph_count is how many
## of S_count that is, so it is not added again.
pop_total <- function(o) with(o, S_count + A_count + D_count + U_count + Tr_count)

## Every human compartment of each block, each indexed by age group first.
PF_COMPARTMENTS <- c("S", "D", "A", "U", "Tr", "Tr_slow", "Tr_c", "Tr_cs", "Ph", "Ph_c", "Ph_ct")
PV_COMPARTMENTS <- c("Sv", "Dv", "Av", "Uv", "Trv", "Trv_slow", "Phv")

## A run's state at the start of every day, stepped as run_simulation_ode() steps
## it (chemoprevention rounds included) but scanned here in R, independently of the
## model's own positivity check: the lowest any human compartment falls to, and
## the first day one falls below -NEG_TOL, with its age group and value (the
## youngest group, if several hold the day's lowest), which is what the check
## reports.
scan_positivity <- function(timesteps, p, tuning = ode_tuning()) {
  tn <- as_ode_tuning(tuning)
  inp <- build_inputs(p, p$init_EIR, tn$age_lower, n_ph = tn$n_ph, n_phc = tn$n_phc,
                      n_sub = tn$n_sub, timesteps = timesteps)
  sys <- dust2::dust_system_create(get_generator(), pars = inp$pars, n_particles = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  uidx <- dust2::dust_unpack_index(sys)
  hum <- c(PF_COMPARTMENTS, PV_COMPARTMENTS)
  idx <- unlist(uidx[hum], use.names = FALSE)
  n_age <- c(rep(inp$pars$n_age, length(PF_COMPARTMENTS)),
             rep(inp$pars$n_age_v, length(PV_COMPARTMENTS)))
  age <- unlist(Map(function(nm, n) rep_len(seq_len(n), length(uidx[[nm]])), hum, n_age),
                use.names = FALSE)
  events <- chemoprevention_events(inp$meta$parameters, timesteps)
  ev_day <- vapply(events, function(e) e$time, numeric(1))
  out <- list(lowest = Inf, day = 0, age = NA, value = NA)
  for (t in seq_len(timesteps)) {
    s <- dust2::dust_system_state(sys, index_state = idx)
    out$lowest <- min(out$lowest, s)
    if (out$day == 0 && min(s) < -NEG_TOL)
      out[c("day", "age", "value")] <- list(t, min(age[s == min(s)]), min(s))
    dust2::dust_system_run_to_time(sys, t)
    for (e in events[ev_day == t]) apply_chemoprevention_pulse(sys, uidx, inp$meta, e)
  }
  out
}

## A parameter list rendering incidence as well as prevalence: malariasimulation's
## defaults render the 2-10 prevalence band and no incidence at all, and fleet
## renders exactly what the list asks for.
gp_bands <- function(overrides = list()) {
  malariasimulation::get_parameters(utils::modifyList(list(
    clinical_incidence_rendering_min_ages = c(0, 730, 0),
    clinical_incidence_rendering_max_ages = c(1825, 3650, 36500),
    severe_incidence_rendering_min_ages = c(730, 0),
    severe_incidence_rendering_max_ages = c(3650, 36500),
    incidence_rendering_min_ages = c(730, 0),
    incidence_rendering_max_ages = c(3650, 36500)), overrides))
}
