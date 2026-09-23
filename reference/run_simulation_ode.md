# Run the mean-field (ODE) malaria model.

Run the mean-field (ODE) malaria model.

## Usage

``` r
run_simulation_ode(
  timesteps,
  parameters = NULL,
  correlations = NULL,
  tuning = ode_tuning(),
  ...
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
  (the individual-mosquito path does not apply).
  [`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
  indexes the mean-field approximations and when each one matters.

  **The target EIR is read from this list**, as `parameters$init_EIR`,
  which is where `malariasimulation` puts it and reads it from too. Seed
  it the same way you would for the IBM, before calling this function:
  `parameters <- malariasimulation::set_equilibrium(parameters, init_EIR = 20)`.
  That call also stores `eq_params`, which `fleet` honours if present.
  Setting `parameters$init_EIR` by hand works but skips that.

  **Call `set_equilibrium()` last.** It freezes all 56 equilibrium
  parameters, 51 of them translated from the IBM list, at the values
  they held when it ran, and `fleet` merges that stored set over the
  live list. So an edit made *after* it – another `set_*()` builder, or
  `parameters$du <- 10` by hand – is a silent no-op for anything in that
  set, and the run is bit-identical to one without the edit. The IBM
  reads those fields directly and would honour the edit, so this is a
  real difference between the two models rather than a detail of this
  one. `fleet` warns when a live value disagrees with the stored one. If
  a calibration loop varies any of them, re-call `set_equilibrium()`
  each iteration.

  Three further optional fields tune how faithfully the mean field
  mirrors the IBM. The first two are per-individual arithmetic, the
  third is a seeding convention; all three default to the validated
  choice, so you normally leave them alone:

  - `bite_dedup` (default `1`): reproduce the IBM's per-timestep bite
    deduplication. malariasimulation collects the day's bitten
    individuals in a bitset, so a person bitten repeatedly in one
    timestep is infected at most once: the daily infection probability
    is `(1 - exp(-EPS)) * b`, which saturates. Setting `0` uses the
    unbounded `b * EPS` instead, the linear form `malariaEquilibrium`
    assumes, which over-predicts infection where exposure approaches one
    bite/person/day (seasonal peaks, high-`zeta` strata). It does
    **not** make the `malariaEquilibrium` seed an exact fixed point, and
    it does not make a run flat: deduplication is only the largest of
    four departures from that solution, which
    [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
    lists in its initial-conditions section. What it does do is roughly
    halve the largest excursion from the seed: over an undisturbed
    15-year run at EIR 20, PfPR(2-10) departs from its seeded value by
    at most 0.14% with `0`, against 0.28% at the default. That is the
    metric the package's own equilibrium tests assert, which is why they
    set it. Leave it at `1` for any scientific run.

  - `acquired_immunity_offset` (default `0`): the IBM adds `+0.5` to
    positive acquired immunity inside the `b`/`phi`/`theta` Hill
    functions. That is a per-individual detail which does not carry over
    to a stratum mean; an A/B against the IBM ensemble mean favours `0`.
    Set `0.5` to reproduce the IBM's literal Hill calls.

  - `hold_init_EIR` (default `FALSE`): only matters with
    `set_demography()`. malariasimulation's `set_equilibrium()` sizes
    the mosquito population from the equilibrium under its *default*
    exponential age structure, so under a custom demography the IBM
    drifts to whatever transmission that density supports. By default
    fleet replicates that: it takes the IBM's mosquito density and seeds
    at the EIR its own equilibrium under the custom age structure then
    supports (a fixed point, so no burn-in), which is generally *not*
    `init_EIR`. Set `TRUE` to seed at `init_EIR` exactly instead.

- correlations:

  accepted so the first three arguments mirror
  `malariasimulation::run_simulation(timesteps, parameters, correlations)`
  exactly (drop-in call compatibility). Intervention correlation is an
  individual-level feature with no mean-field analogue, so a non-NULL
  value is ignored with a warning (the ODE assumes independence between
  interventions). `tuning` is the *fourth* argument: an
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  object (or a list of its fields) passed positionally here would land
  in `correlations` and be thrown away, so that is an error rather than
  a silently defaulted run.

- tuning:

  ODE solver and discretisation settings, from
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md).
  A plain named list of the fields you want to change is also accepted;
  every field left out keeps its default. These are
  numerical-approximation knobs, not model parameters: nothing
  epidemiological lives here.

- ...:

  not used. Named arguments that were top-level before the signature was
  cut back to
  [`malariasimulation::run_simulation()`](https://rdrr.io/pkg/malariasimulation/man/run_simulation.html)'s
  are caught here and reported with the call that replaces them.

## Value

a wide, malariasimulation-style daily count table, one row per output
day. The column names do not change with the parameter set, so runs with
and without interventions are directly comparable and `rbind`-able.
Columns:

- `timestep`: output day (0..timesteps), the row key.

- per age band (tags in **days**, `<lo>_<hi>`): `n_age_*` (population),
  `n_detect_lm_*` and `p_detect_lm_*` (LM-positive count and PfPR
  proportion), `n_detect_pcr_*`, `n_inc_clinical_*`, `n_inc_severe_*`
  and `n_inc_*` (all-infection incidence), all per-day counts. A band is
  aggregated by exact overlap: each model age group contributes the
  fraction of its own width that lies inside the band, so a band need
  not align with the `age_lower` grid, and bands that partition a span
  of the age axis partition the population living in that span exactly.
  The one age group that is not a finite interval is the absorbing top
  group, which holds everyone from its lower edge upward: a band is
  given all of it if it extends *above* that edge and none of it if it
  stops exactly on it. So a band `[0, 36500)` renders the whole
  population on the default grid (top group from 29200 days) but leaves
  the top group out on a grid built with `max_age = 100` – 1.3% of the
  population, at the default demography. That is not silent: the run
  warns, names the share of the population no `n_age_*` band renders,
  and names the age groups it sits in. A band that overlaps no age group
  at all (it lies above the open-ended top group) warns too, and its
  `p_detect_lm_*` is `NA` rather than 0.

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

## See also

[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
for the solver and discretisation settings.

## Examples

``` r
if (FALSE) { # \dontrun{
# exactly the call you would make to malariasimulation::run_simulation()
p <- malariasimulation::set_equilibrium(
  malariasimulation::get_parameters(), init_EIR = 20)
out <- run_simulation_ode(timesteps = 3650, parameters = p)
head(out[, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650", "EIR")])

# a long projection, trading a little accuracy for speed
fast <- run_simulation_ode(3650, p, tuning = list(rtol = 1e-6, step_size_max = 10))
} # }
```
