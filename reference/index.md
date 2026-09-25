# Package index

## Run the model

Simulate malaria transmission dynamics with the deterministic mean-field
model, one day at a time.

- [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
  : Run the mean-field malaria model.

## Discretisation settings

The age grid, the prophylaxis chains’ stage counts and the mosquito
sub-steps per day. Every default is the validated choice: change one
only when you have a reason to.

- [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  : Discretisation settings.

## Age grid

The default graded age discretisation used throughout the model.

- [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md)
  : Default (graded) age grid

## Package overview

Background, scope, and pointers to the mean-field approximations.

- [`fleet`](https://pwinskill.github.io/fleet/reference/fleet-package.md)
  [`fleet-package`](https://pwinskill.github.io/fleet/reference/fleet-package.md)
  : fleet: a mean-field twin of malariasimulation
