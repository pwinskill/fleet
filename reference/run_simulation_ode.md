# Run the mean-field model

Simulates a `malariasimulation` parameter list one day at a time,
resolving each day in the order `malariasimulation` does: immunity
decay, biting and infection, one competing draw per person between
infection and progression, a day of the mosquito model, and deaths
replaced by births. The first three arguments match
[`malariasimulation::run_simulation()`](https://rdrr.io/pkg/malariasimulation/man/run_simulation.html).
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
specifies the model.

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

  number of days to simulate, a single whole number of at least 1. The
  output has one row per day, `1..timesteps`, as `malariasimulation`'s
  does: row `t` holds the state at the start of day `t` and the
  incidence during it, so row 1 is the seed.

- parameters:

  a
  [`malariasimulation::get_parameters()`](https://rdrr.io/pkg/malariasimulation/man/get_parameters.html)
  list, with any interventions added by the `set_*()` functions.
  Required. Everything is read from the list, with no extra arguments:
  clinical treatment (time-varying), antimalarial resistance, bed nets,
  IRS, MDA, SMC and PMC, PEV (EPI and mass) and TBV vaccines,
  seasonality and flexible carrying capacity, and custom demography
  (`set_demography()`: age-specific mortality and the resulting
  equilibrium age structure).

  The parasite is `parameters$parasite`, as in `malariasimulation`: a
  list from `get_parameters(parasite = "vivax")` runs the *P. vivax*
  model, with hypnozoite batches and relapse, radical cure with
  liver-stage protection (`CQ_PQ_params_vivax`, `CQ_TQ_params_vivax`)
  and no severe disease. MDA, SMC and PMC are refused under *P. vivax*,
  because `malariasimulation` itself fails on them. The model is always
  compartmental, so the individual-mosquito option does not apply.
  [`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
  lists the mean-field approximations and when each one matters.

  The target EIR is `parameters$init_EIR`, where `malariasimulation`
  stores and reads it; a list without it is an error. Set it as for the
  IBM, before the run:
  `parameters <- malariasimulation::set_equilibrium(parameters, init_EIR = 20)`.
  That call also stores `eq_params`, which `fleet` uses if present.
  Setting `parameters$init_EIR` by hand works but skips that step.

  Call `set_equilibrium()` last. It stores all 56 equilibrium
  parameters, 51 of them translated from the IBM list, at the values
  they held when it ran, and `fleet` merges that stored set over the
  live list, taking only `iv0` from the live list. An edit made
  afterwards to anything in the set, by another `set_*()` function or by
  hand (`parameters$du <- 10`), has no effect, and the run is
  bit-identical to one without it. The IBM reads those fields directly
  and would apply the edit, so this is a real difference between the two
  models. `fleet` warns when a live value disagrees with the stored one.
  If a calibration loop varies any of them, call `set_equilibrium()`
  again on each iteration.

  Four optional fields, which are not `malariasimulation` parameters,
  control how the mean field mirrors the IBM. Each defaults to the
  validated choice, and is normally left unset:

  - `bite_dedup` (default `1`; *P. falciparum* only, since the *P.
    vivax* model counts every bite, as the IBM does): reproduce the
    IBM's per-timestep bite deduplication. `malariasimulation` collects
    the day's bitten individuals in a bitset, so a person bitten several
    times in one timestep is infected at most once, and the daily
    infection probability is `(1 - exp(-EPS)) * b`, which saturates. `0`
    lets every bite infect independently, `1 - exp(-b * EPS)`, which
    over-predicts infection where exposure approaches one bite per
    person per day (seasonal peaks, high-`zeta` strata). Neither makes
    the `malariaEquilibrium` seed an exact fixed point; the
    initial-conditions section of
    [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
    lists the departures from that solution. Leave it at `1` for any
    scientific run.

  - `acquired_immunity_offset` (default `0` for *P. falciparum*, `0.5`
    for *P. vivax*): the IBM adds 0.5 to positive acquired immunity
    inside the `b`, `phi` and `theta` Hill functions, a per-individual
    detail. The falciparum model reads its curves at a stratum mean,
    where the offset does not carry over and an A/B comparison against
    the IBM ensemble mean favours `0`. The vivax model reads them at
    quadrature nodes that stand for individuals (see `immunity_spread`),
    where it does, so it takes the IBM's `0.5`, or `0` with
    `immunity_spread = FALSE`, back at the cell mean.

  - `immunity_spread` (default `TRUE`; *P. vivax* only): model the
    spread of immunity among people who share an age, heterogeneity
    group and hypnozoite batch count. In the IBM they differ by batch
    history, and the vivax immunity curves are steep enough that the
    mean probability of clinical disease can be several times the
    probability at the mean immunity. `fleet` carries each cell's second
    moment of immunity and averages the curves over a gamma distribution
    with that mean and variance. `FALSE` evaluates them at the cell
    mean, which under-predicts clinical incidence in school-age children
    and young adults by up to 40% once the IBM has built up its spread.

  - `hold_init_EIR` (default `FALSE`): matters only with
    `set_demography()`.
    [`malariasimulation::set_equilibrium()`](https://rdrr.io/pkg/malariasimulation/man/set_equilibrium.html)
    sizes the mosquito population from the equilibrium under the default
    exponential age structure, so under a custom demography the IBM
    drifts to whatever transmission that density supports. By default
    `fleet` does the same: it takes the IBM's mosquito density and seeds
    at the EIR that its own equilibrium under the custom age structure
    supports, which is generally not `init_EIR`. `TRUE` seeds at
    `init_EIR` exactly. Either way the seed is the `malariaEquilibrium`
    solution at that EIR, and it relaxes slightly, as every seed does.

- correlations:

  accepted so that the first three arguments match
  `malariasimulation::run_simulation(timesteps, parameters, correlations)`.
  Intervention correlation is an individual-level feature with no
  mean-field analogue, so a non-`NULL` value is ignored with a warning,
  and interventions are assumed independent. Because `tuning` is the
  fourth argument, an
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  object, or a list of only its fields, passed here is an error rather
  than being silently discarded.

- tuning:

  discretisation settings: an
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  object, or a named list of the fields to change, the rest keeping
  their defaults. These are numerical settings only; model parameters
  stay on `parameters`.

- ...:

  not used. Arguments that were top-level before the signature was
  reduced to
  [`malariasimulation::run_simulation()`](https://rdrr.io/pkg/malariasimulation/man/run_simulation.html)'s
  (`init_EIR`, the fields now in `tuning`, and the retired ODE-solver
  settings) are caught and reported with the call that replaces them.
  Any other argument is an error.

## Value

A data frame in `malariasimulation`'s layout, one row per day. A
parameter list gives the columns the IBM would give it, with the same
names and meanings, as far as a mean field carries them:

- `timestep`: the output day, `1..timesteps`.

- Age-banded counts, one set per band in the parameter list's rendering
  lists and none for an empty list: `n_detect_lm_*`, `p_detect_lm_*`
  (not for *P. vivax*) and `n_detect_pcr_*` over
  `prevalence_rendering_*`; `n_inc_*` and `p_inc_*` over
  `incidence_rendering_*`; `n_inc_clinical_*` and `p_inc_clinical_*`
  over `clinical_incidence_rendering_*`; `n_inc_severe_*` and
  `p_inc_severe_*` over `severe_incidence_rendering_*`; and `n_age_*`
  (population) over all of those and `age_group_rendering_*` together.
  The tag is `<lower>_<upper>` in days, formatted as `malariasimulation`
  formats it. A band `[lower, upper]` holds the whole-day ages `lower`
  to `upper`, both included, as in the IBM, so `0_1824` and `1825_5474`
  partition the age axis while `0_1825` and `1825_5475` share one day's
  cohort.

        In the IBM an `n_*` column is a sampled count and its `p_*` partner the
        expected count behind it: `p_detect_lm_*` is the clinically detected
        count plus each asymptomatic person's probability of detection, a count
        and not a prevalence. A mean field carries expectations, so here the two
        are the same number. The exception is `p_inc_*`, which is 0 because
        `malariasimulation` 3.0.0 renders it from `incidence_min_ages`, a field
        no parameter list sets, and so returns 0 for it too. Prevalence is
        `n_detect_lm_*` / `n_age_*`, as `postie` computes it.

        A band is aggregated by exact overlap with the model's age groups, so it
        need not align with the `age_lower` grid. The exception is the absorbing
        top group, which holds everyone from its lower edge upward: a band is
        given all of it if the band reaches above that edge and none of it
        otherwise. A band ending inside it (60-85 years on the default grid,
        whose top group starts at 80) is therefore given all of it, and a band
        ending exactly on its lower edge none of it, which warns, as does a band
        that selects no age group at all.
      \item `n_infections` (new infections at all ages) and, when clinical
        treatment is deployed, `ft` (the fraction of clinical cases treated), as
        in the IBM.
      \item *P. vivax* only: `n_relapses` (the day's relapse infections) and
        `n_with_hypnozoites` (people carrying at least one batch), and by age
        `n_inc_relapse_*` over `incidence_relapse_rendering_*` and
        `n_with_hypnozoites_*`, with the band's population as `n_*`, over
        `n_with_hypnozoites_rendering_*`. A *P. vivax* `A` is LM-detectable by
        definition, so `n_detect_lm_*` counts all of it, and the IBM renders no
        `p_detect_lm_*` for *P. vivax*. A severe band on a *P. vivax* list
        renders as 0, as the IBM's does.
      \item Population totals by infection state, `S_count`, `A_count`,
        `D_count`, `U_count` and `Tr_count`, counted as the IBM counts them:
        someone a drug is protecting is uninfected, and so in `S_count`.
      \item Immunity means over the population: `ica_mean`, `icm_mean`,
        `ib_mean`, `iva_mean`, `ivm_mean` and `id_mean` for *P. falciparum*;
        `ica_mean`, `icm_mean`, `iaa_mean`, `iam_mean` and `hypnozoites_mean`
        (the mean batch count) for *P. vivax*; and each by age, `<name>_mean_*`
        with the band's population as `n_*`, over its `<name>_rendering_*`. The
        IBM seeds falciparum maternal immunity from the equilibrium bin 0.1
        years older than each person, so for the first months its `icm_mean` and
        `ivm_mean` are below `fleet`'s, whose births take their mothers'.
      \item Three columns `malariasimulation` does not have: `EIR` (infectious
        bites per adult per year), `FOIM` (the force of infection on the first
        mosquito species) and `Ph_count` (how many of `S_count` a drug is
        protecting).

`malariasimulation` columns this table does not carry: `n_bitten`,
`infectivity`, `EIR_<species>`, `FOIM_<species>`, `mu_<species>`, the
mosquito compartment counts and `total_M_<species>`, `natural_deaths`,
and the treatment and intervention tallies (`n_treated`, `n_smc_treated`
and the rest). The IBM renders `n_inc_relapse_*` only on days with a
relapse in the band, and `NA` otherwise; a mean field's expected count
is never zero.

Incidence columns are daily counts, so do not thin rows before passing
the table to
[`postie::get_rates()`](https://rdrr.io/pkg/postie/man/get_rates.html)
or
[`postie::get_prevalence()`](https://rdrr.io/pkg/postie/man/get_prevalence.html),
which accept it directly.

## Details

The seed is where `malariasimulation` seeds: the `malariaEquilibrium`
solution for `init_EIR` (`malariaEquilibriumVivax` for *P. vivax*), with
the humans under the treatment coverage in force at timestep 0 and the
mosquitoes sized under the coverage at timestep 1, as
`set_equilibrium()` and the IBM's initialisation do. It is close to the
daily model's fixed point but not on it, so outputs relax over the first
years. A seasonal run is seeded at the annual-mean (aseasonal)
equilibrium and oscillates around a limit cycle rather than holding
flat. Use a burned-in window for calibration and comparison;
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
gives the relaxation and settling times.

A run in which any compartment falls below zero stops with an error
naming the day and the age group. Each day an age group loses its ageing
and mortality fraction out of the same stock as its infections, and
several consecutive groups a day or two wide can lose more than they
hold at an EIR in the hundreds; widen them in `age_lower`
([`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)).

## See also

[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
for the discretisation settings.

## Examples

``` r
if (FALSE) { # \dontrun{
# the same call as malariasimulation::run_simulation()
p <- malariasimulation::set_equilibrium(
  malariasimulation::get_parameters(), init_EIR = 20)
out <- run_simulation_ode(timesteps = 3650, parameters = p)
head(out[, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650", "EIR")])

# the same run on a 53-group age grid, over twice as fast
coarse <- run_simulation_ode(3650, p,
  tuning = list(age_lower = default_age_lower(n_group = 53)))
} # }
```
