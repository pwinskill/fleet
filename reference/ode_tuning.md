# ODE solver and discretisation settings.

Numerical settings for
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md).
Every default is the validated choice; a run with `ode_tuning()`
untouched is the reference configuration. Pass either this object or a
plain named list of the fields you want to change; anything left out
keeps its default.

## Usage

``` r
ode_tuning(
  age_lower = default_age_lower(),
  n_eir = 10L,
  n_foim = 10L,
  n_eip = 20L,
  n_ph = NULL,
  n_phc = NULL,
  atol = 1e-08,
  rtol = 1e-08,
  step_size_max = 1,
  odin_file = NULL
)
```

## Arguments

- age_lower:

  age-group lower edges in **years**. The default is
  [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md),
  a 53-group grid log-spaced between pinned reporting ages. Must start
  at 0, increase strictly, and stay in years: a top edge above 1000 is
  rejected as a grid supplied in days.

  The default is **not** converged for severe disease or for adult
  clinical incidence, and both carry a discretisation bias of the order
  of a few per cent on it. That bias is removable: `fleet`'s own profile
  converges at first order, so `default_age_lower(n_group = 105)`
  roughly halves it and `n_group = 209` halves it again, at a
  proportionate cost in run time. A result that rests on the *level* of
  severe incidence or on adult bands is worth re-running on a finer grid
  to see how much of it is the grid.

  What a finer grid will **not** do is close the remaining gap to the
  IBM. That is the mean field itself – one immunity value per stratum
  against a spread of individual infection histories at the same age –
  and it does not shrink with group width. See the
  `age-profile-clinical` claim in `fleetcheck` for the measurement.

  The default is left where it is because every published comparison is
  stated on it; changing it moves every number.

- n_eir, n_foim, n_eip:

  Erlang-chain stage counts for the EIR lag, FOIM lag and mosquito EIP.
  Larger values sharpen the (otherwise gamma-shaped) lags toward the
  IBM's fixed delays; equilibrium is exact for any value.

- n_ph, n_phc:

  Erlang-chain stage counts for the post-treatment (`Ph`) and
  chemoprevention (`Ph_c`) prophylaxis compartments. `NULL` (default)
  matches the chain's variance to the drug's Weibull protection curve,
  capped at 20. For `Ph_c` that is `1/CV²` of the Weibull: 14 for SP-AQ,
  15 for DHA-PQP. `Ph` follows the exponential treated stage `Tr`, so
  its count matches the variance of the whole `Tr + Ph` sojourn and its
  mean is the integrated protection left after `Tr`: 16 stages for
  SP-AQ, 20 for DHA-PQP, and 1 for AL, whose 10-day protection is
  already less variable than `Tr` itself. A drug mixture is
  moment-matched as a mixture. `1` is a single exponential stage, which
  for `Ph_c` leaks protection early between monthly SMC rounds. The
  count is fixed at the seed's drug mix: a first-line switch moves the
  chain's mean, not its shape.

- atol, rtol, step_size_max:

  dust2 ODE-solver controls. The defaults (`1e-8`, `1e-8`, `1`) preserve
  the flat equilibrium exactly. For long dynamic projections
  `rtol = 1e-6` runs about 1.4x faster on seasonal ones (7.2 s to 5.0 s
  over 30 years) with a maximum deviation of 4e-7 in daily clinical
  incidence – five orders of magnitude inside the IBM's replicate band –
  and not at all faster on aseasonal ones.

  That last point is a hard floor, and it is worth knowing where it
  comes from. On an equilibrium run every accepted step is 0.62 days,
  and that step does not move with `atol`, `rtol`, `step_size_max`, the
  output grid, or the delay-chain lengths below 4 per day. A step that
  ignores both tolerances is a stability limit, not an accuracy one: the
  stiffest eigenvalue in the system is the late-larval density-dependent
  mortality, `ml * gamma * (E + 2L) / K`, which at the seed is 5.3 per
  day (gamma is 13.25 and the early-larval stock sits at about ten times
  K), and Dormand-Prince's stability boundary of ~3.3 on that mode gives
  3.3 / 5.31 = 0.621 days. It is intrinsic to malariasimulation's larval
  model at its default parameters, so the 1.6 steps per day it forces is
  what an explicit stepper costs here, and only an implicit one (dust2
  has none) could take longer steps at equilibrium. The daily output
  grid adds about 2,700 rejected trial steps on top over 30 years – the
  controller re-grows the step after each forced stop – which is why
  output every 30 days is ~25% faster on aseasonal runs and no more.

  **`step_size_max` is a safety rail, not a speed control**, despite
  travelling with `rtol` in the preset above. On 30-year runs, seasonal
  and not, the solver takes the same number of steps at `1`, `5`, `30`
  and `Inf` to within four in nineteen thousand, and the outputs agree
  to 1e-10. What the cap is for is stopping a trial step overshooting
  the end of an interpolation grid. The intervention series are
  interpolated on grids that extend just past `timesteps`, and near
  equilibrium the stepper will propose a step of tens of days; uncapped,
  such a step could land beyond the end of a coarse grid and abort the
  run ("Tried to interpolate at time = ..., which is ... after the last
  time"). None of the runs measured here did, but the cap costs nothing
  measurable, so leave it at `1` unless you have a specific reason.

  Keep `atol` at `1e-8`: the individual prophylaxis chain stages hold
  occupancies of order `1e-6`, which a looser absolute tolerance lets
  dip below zero.

- odin_file:

  optional path to an odin source to compile instead of the generator
  built into the package (development use; needs 'odin2' and a C++
  toolchain).

## Value

a `fleet_ode_tuning` list.

## See also

[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)

## Examples

``` r
ode_tuning()$atol
#> [1] 1e-08
# a faster long projection: the tolerance is the lever, not the step cap
ode_tuning(rtol = 1e-6)$rtol
#> [1] 1e-06
```
