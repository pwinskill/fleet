## The target EIR reaches blink on the parameter list, exactly as it reaches the
## IBM -- run_simulation_ode() has no init_EIR argument. This alias keeps a test
## call on one line; it is nothing but set_equilibrium(), and tests that care
## about the seeding convention itself should spell that out instead.
eqm <- function(p, eir) malariasimulation::set_equilibrium(p, init_EIR = eir)
