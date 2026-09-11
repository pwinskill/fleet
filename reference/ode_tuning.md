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

  age-group lower edges in **years** (default graded grid: monthly to 1
  year, quarterly to 5, yearly to 15, then 5-yearly to an absorbing top
  group). Must start at 0, increase strictly, and stay in years: a top
  edge above 1000 is rejected as a grid supplied in days.

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
  the flat equilibrium exactly; for long dynamic projections
  `rtol = 1e-6` with `step_size_max = 10` runs about 1.4–1.7x faster on
  seasonal projections, and barely faster on aseasonal ones (there the
  daily output grid, not the tolerance, sets the step count), with
  negligible effect on aggregate outputs. Keep `atol` at `1e-8`: the
  individual prophylaxis chain stages hold occupancies of order `1e-6`,
  which a looser absolute tolerance lets dip below zero.

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
# a faster long projection
ode_tuning(rtol = 1e-6, step_size_max = 10)$step_size_max
#> [1] 10
```
