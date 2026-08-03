# Post-process an ODE run into postie-format rates and prevalence.

Returns a list with a `rates` table (long: one row per timestep x age
band, with `clinical`, `severe`, `mortality`, `yld`, `yll`, `dalys`,
`person_days`) and a `prevalence` table (wide: one
`<diagnostic>_prevalence_<lo>_<hi>` column per age band), matching
[`postie::get_rates()`](https://rdrr.io/pkg/postie/man/get_rates.html) /
[`postie::get_prevalence()`](https://rdrr.io/pkg/postie/man/get_prevalence.html).

## Usage

``` r
get_epi_outputs(
  x,
  diagnostic = "lm",
  rates_args = list(),
  prevalence_args = list(),
  ...
)
```

## Arguments

- x:

  output of
  [`run_simulation_ode()`](https://pwinskill.github.io/blink/reference/run_simulation_ode.md).

- diagnostic:

  prevalence diagnostic, "lm" or "pcr".

- rates_args:

  named list of extra arguments for
  [`postie::get_rates()`](https://rdrr.io/pkg/postie/man/get_rates.html)
  (e.g. `scaler`, `treatment_scaler`, `life_expectancy`, `infer_ft`).

- prevalence_args:

  named list of extra arguments for
  [`postie::get_prevalence()`](https://rdrr.io/pkg/postie/man/get_prevalence.html).

- ...:

  shared arguments accepted by both postie functions (e.g.
  `baseline_year`, `ages_as_years`).

## Value

list(rates, prevalence).

## Examples

``` r
if (FALSE) { # \dontrun{
p <- malariasimulation::get_parameters()
out <- run_simulation_ode(timesteps = 3650, parameters = p, init_EIR = 20)
epi <- get_epi_outputs(out, diagnostic = "lm")
head(epi$rates)
epi$prevalence$lm_prevalence_2_10
} # }
```
