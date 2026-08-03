# Run the mean-field (ODE) malaria model.

Run the mean-field (ODE) malaria model.

## Usage

``` r
run_simulation_ode(
  timesteps,
  parameters = NULL,
  correlations = NULL,
  init_EIR = NULL,
  age_lower = default_age_lower(),
  n_eir = 10L,
  n_foim = 10L,
  n_eip = 20L,
  atol = 1e-08,
  rtol = 1e-08,
  step_size_max = 1,
  odin_file = NULL
)
```

## Arguments

- timesteps:

  number of days to simulate. Output has `timesteps + 1` rows (day 0 is
  the seeded equilibrium; days 1..timesteps are integrated).

- parameters:

  a malariasimulation::get_parameters() list (falciparum), optionally
  with interventions layered on via the `set_*` builders. The following
  are applied automatically from the list (no extra arguments): clinical
  treatment (time-varying), antimalarial resistance, bed nets, IRS,
  MDA/SMC/PMC, PEV (EPI + mass) and TBV vaccines, seasonality and
  flexible carrying capacity, and custom demography (`set_demography`:
  age-specific mortality and the resulting equilibrium age
  structure). P. vivax is rejected; the model is always compartmental
  (the individual-mosquito path does not apply). See the README for the
  mean-field approximations used by each module.

- correlations:

  accepted so the first three arguments mirror
  `malariasimulation::run_simulation(timesteps, parameters, correlations)`
  exactly (drop-in call compatibility). Intervention correlation is an
  individual-level feature with no mean-field analogue, so a non-NULL
  value is ignored with a warning (the ODE assumes independence between
  interventions).

- init_EIR:

  target adult EIR (bites/adult/year). If NULL, taken from
  `parameters$init_EIR` (set by malariasimulation::set_equilibrium()).

- age_lower:

  age-group lower edges in **years** (default graded grid).

- n_eir, n_foim, n_eip:

  Erlang-chain stage counts for the EIR lag, FOIM lag and mosquito EIP.
  Larger values sharpen the (otherwise gamma-shaped) lags toward the
  IBM's fixed delays; equilibrium is exact for any value.

- atol, rtol, step_size_max:

  dust2 ODE-solver controls. The defaults (`1e-8`, `1e-8`, `1`) preserve
  the flat equilibrium exactly; for long dynamic projections a looser
  tolerance and larger step cap (e.g. `1e-6`, `1e-6`, `10`) run several
  times faster with negligible effect on aggregate outputs.

- odin_file:

  optional path to the odin source (development use).

## Value

a wide, malariasimulation-style daily count table, one row per output
day. Columns:

- `timestep` — output day (0..timesteps), the row key.

- per age band (tags in **days**, `<lo>_<hi>`): `n_age_*` (population),
  `n_detect_lm_*` and `p_detect_lm_*` (LM-positive count and PfPR
  proportion), `n_detect_pcr_*`, `n_inc_clinical_*`, `n_inc_severe_*`
  and `n_inc_*` (all-infection incidence) — all per-day counts.

- `ft` (treated fraction), `EIR` (per adult per year), `FOIM`.

- population-total infection-state counts `S_count`, `D_count`,
  `A_count`, `U_count`, `Tr_count`, `Ph_count` (diagnostic).

Incidence columns are per-day counts; do not thin rows before postie.
Can be passed directly to postie::get_rates() /
postie::get_prevalence().

## Details

Seasonal runs oscillate around a limit cycle rather than holding flat;
the state is seeded at the annual-mean (aseasonal) equilibrium, so the
first ~10 years are a transient onto the cycle and the seasonal
annual-mean EIR sits a few percent below the aseasonal `init_EIR` target
(nonlinear averaging). Use a burned-in cycle for calibration/comparison.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- malariasimulation::get_parameters()

# daily wide count table, init_EIR supplied directly
out <- run_simulation_ode(timesteps = 3650, parameters = p, init_EIR = 20)
head(out[, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650", "EIR")])

# or let set_equilibrium() seed init_EIR into the parameter list
p2 <- malariasimulation::set_equilibrium(
  malariasimulation::get_parameters(), init_EIR = 5)
out2 <- run_simulation_ode(timesteps = 3650, parameters = p2)
} # }
```
