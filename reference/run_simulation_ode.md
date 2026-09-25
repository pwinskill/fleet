# Run the mean-field malaria model.

The population is advanced one day at a time, in the order
malariasimulation resolves a day: immunity decay, biting and infection,
one competing draw per person between infection and progression, a day
of the mosquito model, and deaths replaced by births. See
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md).

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

  number of days to simulate. Output has one row per day,
  `1..timesteps`, as malariasimulation's does: row `t` holds the state
  at the start of day `t` and the incidence during it, so row 1 is the
  seed.

- parameters:

  a malariasimulation::get_parameters() list, optionally with
  interventions layered on via the `set_*` builders. The following are
  applied automatically from the list (no extra arguments): clinical
  treatment (time-varying), antimalarial resistance, bed nets, IRS,
  MDA/SMC/PMC, PEV (EPI + mass) and TBV vaccines, seasonality and
  flexible carrying capacity, and custom demography (`set_demography`:
  age-specific mortality and the resulting equilibrium age structure).
  The parasite is `parameters$parasite`, as in malariasimulation: a list
  from `get_parameters(parasite = "vivax")` runs the P. vivax model –
  hypnozoite batches and relapse, radical cure with liver-stage
  protection (`CQ_PQ_params_vivax`, `CQ_TQ_params_vivax`), no severe
  disease. MDA, SMC and PMC are refused under vivax because
  malariasimulation itself fails on them. The model is always
  compartmental (the individual-mosquito path does not apply).
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

  Four further optional fields tune how faithfully the mean field
  mirrors the IBM: two per-individual details of its arithmetic, the
  vivax immunity closure, and a seeding convention. All four default to
  the validated choice, so you normally leave them alone:

  - `bite_dedup` (default `1`; P. falciparum only, since vivax counts
    every bite, as the IBM does): reproduce the IBM's per-timestep bite
    deduplication. malariasimulation collects the day's bitten
    individuals in a bitset, so a person bitten repeatedly in one
    timestep is infected at most once: the daily infection probability
    is `(1 - exp(-EPS)) * b`, which saturates. Setting `0` lets every
    bite infect independently instead, `1 - exp(-b * EPS)`, which
    over-predicts infection where exposure approaches one
    bite/person/day (seasonal peaks, high-`zeta` strata). Neither makes
    the `malariaEquilibrium` seed an exact fixed point:
    [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
    lists the departures from that solution in its initial-conditions
    section. Leave it at `1` for any scientific run.

  - `acquired_immunity_offset` (default `0` for P. falciparum, `0.5`
    for P. vivax): the IBM adds `+0.5` to positive acquired immunity
    inside the `b`/`phi`/`theta` Hill functions. That is a
    per-individual detail. The falciparum model reads its curves at a
    stratum mean, where it does not carry over and an A/B against the
    IBM ensemble mean favours `0`; the vivax model reads them at
    quadrature nodes that stand for individuals (see `immunity_spread`),
    where it does, so it takes the IBM's `0.5` – and `0` with
    `immunity_spread = FALSE`, back at the cell mean.

  - `immunity_spread` (default `TRUE`; P. vivax only): model the spread
    of immunity among people who share an age, heterogeneity group and
    hypnozoite batch count. In the IBM they differ by their batch
    *history*, and the vivax immunity curves are steep enough that the
    mean probability of clinical disease can be several times the
    probability at the mean immunity. fleet carries each cell's second
    moment of immunity and averages the curves over a gamma with that
    mean and variance. `FALSE` evaluates them at the cell mean, which
    under-predicts clinical incidence in school-age children and young
    adults by up to 40% once the IBM has built up its spread.

  - `hold_init_EIR` (default `FALSE`): only matters with
    `set_demography()`. malariasimulation's `set_equilibrium()` sizes
    the mosquito population from the equilibrium under its *default*
    exponential age structure, so under a custom demography the IBM
    drifts to whatever transmission that density supports. By default
    fleet replicates that: it takes the IBM's mosquito density and seeds
    at the EIR its own equilibrium under the custom age structure then
    supports, which is generally *not* `init_EIR`. Set `TRUE` to seed at
    `init_EIR` exactly instead. Either way the seed is the
    `malariaEquilibrium` solution at that EIR, and it relaxes slightly,
    as every seed does.

- correlations:

  accepted so the first three arguments mirror
  `malariasimulation::run_simulation(timesteps, parameters, correlations)`
  exactly (drop-in call compatibility). Intervention correlation is an
  individual-level feature with no mean-field analogue, so a non-NULL
  value is ignored with a warning (the model assumes independence
  between interventions). `tuning` is the *fourth* argument: an
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  object (or a list of its fields) passed positionally here would land
  in `correlations` and be thrown away, so that is an error rather than
  a silently defaulted run.

- tuning:

  discretisation settings, from
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

a wide daily table in malariasimulation's layout, one row per day. The
columns follow malariasimulation's rendering rules: a parameter list
gives the columns the IBM would give it, under the same names and with
the same meanings, as far as a mean field carries them. Columns:

- `timestep`: output day (1..timesteps), the row key.

- age-banded counts, one set per band in the parameter list's rendering
  lists and none for an empty list: `n_detect_lm_*`, `p_detect_lm_*`
  (not for vivax) and `n_detect_pcr_*` over `prevalence_rendering_*`;
  `n_inc_*` and `p_inc_*` over `incidence_rendering_*`;
  `n_inc_clinical_*` and `p_inc_clinical_*` over
  `clinical_incidence_rendering_*`; `n_inc_severe_*` and
  `p_inc_severe_*` over `severe_incidence_rendering_*`; and `n_age_*`
  (population) over all of those and `age_group_rendering_*` together.
  The tag is `<lower>_<upper>` in days, formatted as malariasimulation
  formats it. A band `[lower, upper]` holds the whole-day ages `lower`
  to `upper`, both included, as in the IBM, so `0_1824` and `1825_5474`
  partition the age axis while `0_1825` and `1825_5475` share one day's
  cohort.

        In the IBM an `n_*` column is a sampled count and its `p_*` partner the
        expected count behind it: `p_detect_lm_*` is the clinically detected
        count plus each asymptomatic person's probability of detection, a count
        and not a prevalence. A mean field carries expectations, so here the two
        are the same number. The exception is `p_inc_*`, which is 0 because
        malariasimulation 3.0.0 renders it from `incidence_min_ages`, a field no
        parameter list sets, and so returns 0 for it too. Prevalence is
        `n_detect_lm_*` / `n_age_*`, as postie computes it.

        A band is aggregated by exact overlap with the model's age groups, so it
        need not align with the `age_lower` grid. The one group that is not a
        finite interval is the absorbing top group, which holds everyone from its
        lower edge upward: a band is given all of it if it reaches above that edge
        and none of it otherwise. A band ending inside it (60-85 years on the
        default grid, whose top group starts at 80) is therefore given all of it,
        and a band ending exactly on its lower edge none of it, which warns, as
        does a band that selects no age group at all.
      \item `n_infections` (new infections at all ages) and, when clinical
        treatment is deployed, `ft` (the fraction of clinical cases treated), as
        in the IBM.
      \item P. vivax: `n_relapses` (the day's relapse infections) and
        `n_with_hypnozoites` (people carrying at least one batch), and by age
        `n_inc_relapse_*` over `incidence_relapse_rendering_*` and
        `n_with_hypnozoites_*`, with the band's population as `n_*`, over
        `n_with_hypnozoites_rendering_*`. A vivax `A` is LM-detectable by
        definition, so `n_detect_lm_*` counts all of it, and the IBM renders no
        `p_detect_lm_*` for vivax; a severe band on a vivax list renders, as
        the IBM's does, as 0.
      \item population totals by infection state, `S_count`, `A_count`,
        `D_count`, `U_count` and `Tr_count`, counted as the IBM counts them:
        someone a drug is protecting is uninfected, and so in `S_count`.
      \item the immunity means over the population: `ica_mean`, `icm_mean`,
        `ib_mean`, `iva_mean`, `ivm_mean` and `id_mean` for falciparum,
        `ica_mean`, `icm_mean`, `iaa_mean`, `iam_mean` and `hypnozoites_mean`
        (the mean batch count) for vivax; and each by age, `<name>_mean_*` with
        the band's population as `n_*`, over its `<name>_rendering_*`. The IBM
        seeds falciparum's maternal immunity from the equilibrium bin 0.1 years
        older than each person, so for the first months its `icm_mean` and
        `ivm_mean` sit below `fleet`'s, whose births take their mothers'.
      \item three columns malariasimulation does not have: `EIR` (infectious
        bites per adult per year), `FOIM` (the force of infection on the first
        mosquito species) and `Ph_count` (how many of `S_count` a drug is
        protecting).

malariasimulation columns this table does not carry: `n_bitten`,
`infectivity`, `EIR_<species>`, `FOIM_<species>`, `mu_<species>`, the
mosquito compartment counts and `total_M_<species>`, `natural_deaths`,
and the treatment and intervention tallies (`n_treated`, `n_smc_treated`
and the rest). The IBM renders `n_inc_relapse_*` only on days with a
relapse in the band, and NA otherwise; a mean field's expected count is
never zero.

Incidence columns are per-day counts; do not thin rows before postie.
The table can be passed directly to postie::get_rates() /
postie::get_prevalence().

## Details

The seed is where malariasimulation's own seed is: the
`malariaEquilibrium` solution for `init_EIR`, the humans under the
treatment coverage in force at timestep 0 and the mosquitoes sized under
the coverage at timestep 1, as `set_equilibrium()` and the IBM's
initialisation do. It is close to the daily model's fixed point but not
on it: prevalence relaxes by under half a per cent, incidence by a few
per cent over the first years, more at high transmission. Seasonal runs
oscillate around a limit cycle rather than holding flat; the state is
seeded at the annual-mean (aseasonal) equilibrium, so the first ~10
years are a transient onto the cycle and the seasonal annual-mean EIR
sits about 7% below the aseasonal `init_EIR` target (nonlinear
averaging). Use a burned-in window for calibration and comparison.

## See also

[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
for the discretisation settings.

## Examples

``` r
if (FALSE) { # \dontrun{
# exactly the call you would make to malariasimulation::run_simulation()
p <- malariasimulation::set_equilibrium(
  malariasimulation::get_parameters(), init_EIR = 20)
out <- run_simulation_ode(timesteps = 3650, parameters = p)
head(out[, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650", "EIR")])

# the same run on a coarser age grid: a quarter of the groups, four times as fast
coarse <- run_simulation_ode(3650, p,
  tuning = list(age_lower = default_age_lower(n_group = 53)))
} # }
```
