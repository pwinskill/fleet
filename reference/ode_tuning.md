# Discretisation settings

Numerical settings for
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md),
passed as its `tuning` argument. Every default is the validated choice,
and a run with the defaults is the reference configuration.
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
also accepts a plain named list of the fields to change: fields left out
keep their defaults, and an unknown field is an error.

## Usage

``` r
ode_tuning(
  age_lower = default_age_lower(),
  n_ph = NULL,
  n_phc = NULL,
  n_sub = 32L,
  odin_file = NULL,
  ...,
  batch_bins = TRUE
)
```

## Arguments

- age_lower:

  age-group lower edges in years. The default is
  [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md),
  a 118-group grid log-spaced between pinned reporting ages and weighted
  towards the under-fives; its documentation gives the grid's accuracy
  and cost. Must be finite, start at 0 and increase strictly, with at
  least two edges, and must be in years: a top edge above 1000 is
  rejected as a grid given in days.

  Each day every group loses its ageing and mortality fraction out of
  the same stock as the day's other exits, so the narrowest group must
  be wide enough to hold them together. Progression and clearance do not
  depend on transmission, and a grid with a group too narrow for them, a
  little over a day, is an error before the run. Infection is checked as
  the run goes: a run in which any compartment falls below zero stops
  with an error naming the day and the age group. That takes several
  consecutive groups a day or two wide at an EIR in the hundreds; the
  default grid's narrowest groups are 16.6 days wide.

- n_ph, n_phc:

  stage counts of the post-treatment (`Ph`) and chemoprevention (`Ph_c`)
  prophylaxis chains. `NULL` (default) matches each chain to the drug's
  Weibull protection curve. A stage is left with a fixed probability
  each day, so it lasts a geometric number of days, and a chain matches
  the protection's mean `M` and variance `V` with `M^2 / (V + M)`
  stages: 10 for SP-AQ as chemoprevention, 9 for DHA-PQP. `Ph` follows
  the treated stage `Tr`, so its count matches the variance of the whole
  `Tr + Ph` sojourn and its mean is the protection left after `Tr`: 9
  stages for SP-AQ, 10 for DHA-PQP, and 1 for AL, whose 10-day
  protection is already less variable than `Tr` itself. A drug mixture
  is moment-matched as a mixture.

  A stage lasts at least a day and also loses its ageing and mortality
  fraction each day, so the default is limited to what fits the
  shortest-protecting drug on the schedule and capped at 20, and an
  explicit count that does not fit is an error. The count is fixed for
  the run, sized at the seed's drug mix, so a first-line switch moves
  the chain's mean, not its shape.

- n_sub:

  sub-steps per day of the mosquito model, integrated by an exponential
  midpoint rule (second order, and stable however stiff the larval
  equations become). Once a run is on its cycle the error falls fourfold
  with each doubling. At the default of 32, a seasonal run's daily EIR
  stays within 2.4e-4 of the converged solution at the seasonal trough
  (5.5e-6 of its peak), comparable to the per-step tolerance to which
  `malariasimulation` integrates its own mosquito model
  (`r_tol = a_tol = 1e-4`); 8 leaves it 3.6e-3 off, a slight shift in
  the phase of the seasonal cycle. Annual means agree to about 4e-6 at
  32 and 5e-5 at 8, and a run with no seasonality and no interventions
  is exact to 1e-11 at any count. The mosquito model is a small part of
  the day, so 32 costs 2 to 13% more than 8.

- odin_file:

  optional path to an odin source to compile instead of the generator
  built into the package: a discrete-time model with the built-in one's
  parameters and states. For development; needs 'odin2' and a C++
  toolchain. The compiled model is cached by path for the session, so an
  edit to the file needs a new session.

- ...:

  the ODE solver's controls (`atol`, `rtol`, `step_size_max`) and the
  Erlang lag stage counts (`n_eir`, `n_foim`, `n_eip`) are accepted with
  a warning and ignored: the model advances one day at a time, with the
  IBM's own delay lines, so there is no solver to control and no lag to
  discretise. Any other name is an error.

- batch_bins:

  under *P. vivax*, whether hypnozoite batch counts share levels.
  `malariasimulation` gives each person a count of hypnozoite batches, 0
  to `kmax`, and `fleet` holds every vivax compartment and immunity
  stock at each level of that ladder, so its length accounts for most of
  a vivax run's cost. With `TRUE` (default), counts 0 to 3 and `kmax`
  have a level each and the counts between share levels in pairs, the
  people within a pair spread as a Poisson fitted to the levels either
  side. That gives eight levels for the default `kmax = 10`, which take
  a quarter less time than eleven and move no output of the `fleetcheck`
  vivax suite by more than 0.4%. `FALSE` gives every count its own
  level. Ignored under *P. falciparum*.

## Value

A list of class `fleet_ode_tuning`.

## See also

[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md),
[`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md)

## Examples

``` r
ode_tuning()$n_sub
#> [1] 32
# a finer age grid, to check how much of a severe-incidence level is the grid's
length(ode_tuning(age_lower = default_age_lower(n_group = 417))$age_lower)
#> [1] 417
```
