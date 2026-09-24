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
