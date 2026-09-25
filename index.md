# fleet

[![R-CMD-check](https://github.com/pwinskill/fleet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pwinskill/fleet/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/pwinskill/fleet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/pwinskill/fleet/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Status: work in
progress](https://img.shields.io/badge/status-work%20in%20progress-orange.svg)](#not-ready-for-real-use)
[![Not for production
use](https://img.shields.io/badge/not%20for-production%20use-red.svg)](#not-ready-for-real-use)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/pwinskill/fleet/blob/main/LICENSE.md)

> ## Not ready for real use
>
> ⚠️ **`fleet` is a work in progress and is not validated for research,
> policy or operational use.** It is published early so the approach and
> the comparison against `malariasimulation` can be looked at and argued
> with, not so that anyone can rely on its numbers.
>
> Concretely, and honestly:
>
> - **The API is unstable.**
>   [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)’s
>   signature has already changed once and may change again without
>   deprecation. **Known discrepancies against the IBM are open, not
>   resolved.** Across the 63-country site-file comparison `fleet` runs
>   above the IBM on falciparum clinical and severe incidence, an excess
>   that is **not** explained. On *P. vivax* sites it runs further above
>   the IBM, and that is not explained either; and at school age its
>   vivax LM prevalence and clinical incidence run a few per cent high,
>   a limit of the mean field. The current verdict on every claim, with
>   the numbers behind it, is in
>   *[fleetcheck](https://pwinskill.github.io/fleetcheck/)* — which is
>   kept current, unlike any figure quoted here would be.
> - **Severe incidence and anything derived from it (including DALYs)
>   should be treated as indicative only.**
> - Nothing here has been peer reviewed, and there is no versioned
>   release.
>
> If you need results you can defend today, use
> [malariasimulation](https://github.com/mrc-ide/malariasimulation).

> A fast, deterministic **mean-field twin** of the
> [malariasimulation](https://github.com/mrc-ide/malariasimulation)
> individual-based model of *Plasmodium falciparum* and *P. vivax*
> malaria, run on its daily clock: same inputs, same outputs,
> population-independent.

## What it is

`fleet` reproduces the Griffin-style human model, structured by age and
biting heterogeneity, coupled to the compartmental mosquito model, and
for *P. vivax* the White-style model with its hypnozoite batches and
relapse. Like `malariasimulation`, it advances one day at a time: every
day it applies immunity decay, biting, one competing draw per person
between infection and progression, a day of the mosquito model, and
deaths and births, in the IBM’s order. It is written in
[odin2](https://github.com/mrc-ide/odin2) /
[dust2](https://github.com/mrc-ide/dust2).

|  |  |
|----|----|
| **Human states** | `S / D / A / U / Tr`, plus the post-treatment and chemoprevention prophylaxis chains and the chemoprevention treated phase; for vivax each over 0 to 10 hypnozoite batches, with liver-stage protection after radical cure |
| **Immunity** | falciparum: four acquired states `IB / ICA / ID / IVA` and two maternal terms `ICM / IVM`, algebraic rather than state variables; vivax: `IAA / ICA` with their within-cell spread and refractory windows, and maternal `IAM / ICM` |
| **Mosquito, per species** | `E / L / P / Sm / Em / Im`, with the IBM’s own delay lines for incubation, the EIR and human infectivity |

It exists to give the malariasimulation ecosystem a **deterministic,
Monte-Carlo-free companion** that:

- **Takes the same inputs as the individual-based model (IBM).**
  [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
  accepts a
  [`malariasimulation::get_parameters()`](https://rdrr.io/pkg/malariasimulation/man/get_parameters.html)
  list, with the usual `set_*` intervention builders layered on,
  unchanged, for either parasite: `get_parameters(parasite = "vivax")`
  runs the vivax model.
- **Is seeded where the IBM is seeded.** Initial conditions come from
  `malariaEquilibrium`, or `malariaEquilibriumVivax` for vivax, under
  the same treatment coverage the IBM seeds with. An undisturbed
  falciparum run’s prevalence relaxes off that seed by under half a
  percent and then holds, while incidence relaxes by a few per cent over
  the first years, more at high transmission. A vivax run moves further,
  as the IBM’s does: its clinical incidence settles about 20% above the
  seed and its realised EIR 4 to 17% above `init_EIR`, taking up to two
  decades at the lowest EIRs. So burn in before comparing levels, as
  with the IBM.
- **Returns malariasimulation’s output table.** A parameter list gives
  the columns the IBM would give it, under the same names and with the
  same meanings, so a post-processing pipeline written for the IBM,
  postie included, works on a `fleet` run unchanged.
- **Is fast and population-independent.** All compartments are
  per-capita densities, so a run costs the same whether you model a
  thousand people or ten million: a 30-year falciparum run takes about 3
  s on the default 209-group age grid (0.7 s on 53 groups), a vivax run,
  with its hypnozoite dimension, about a minute.

Reach for the IBM instead when you need stochastic variation or
individual heterogeneity beyond the mean field.

## Install

`fleet` compiles C++ at install time, as do several of its dependencies,
so you need a working C++ toolchain first: Rtools on Windows, the Xcode
command line tools on macOS, the usual build tools (`r-base-dev` or
equivalent) on Linux.

``` r

# install.packages("remotes")
remotes::install_github("pwinskill/fleet")
```

That also installs the GitHub-only hard dependencies (`dust2`, `monty`,
`malariaEquilibrium` and `malariaEquilibriumVivax`), which `DESCRIPTION`
`Remotes` points at. It installs nothing from `Suggests`; the examples
need two of those:

``` r

remotes::install_github(c("mrc-ide/malariasimulation", "mrc-ide/postie"))
```

`malariasimulation` builds the parameter list and `postie`
post-processes the output.

## Quick start

``` r

library(fleet)

# ask for incidence by age as you would of the IBM: fleet renders exactly the
# bands the list asks for (the defaults give 2-10 prevalence and no incidence)
bands <- list(min = c(0, 5, 15) * 365, max = c(5, 15, 100) * 365 - 1)
p <- malariasimulation::get_parameters(list(
  clinical_incidence_rendering_min_ages = bands$min,
  clinical_incidence_rendering_max_ages = bands$max,
  severe_incidence_rendering_min_ages = bands$min,
  severe_incidence_rendering_max_ages = bands$max))

# 10-year daily table, in malariasimulation's columns
out <- run_simulation_ode(timesteps = 3650, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))

# postie-format rates and prevalence, exactly as for an IBM run
postie::get_prevalence(out, diagnostic = "lm")$lm_prevalence_2_10
postie::get_rates(out)[, c("time", "age_lower", "age_upper", "clinical", "severe", "dalys")]
```

Layer interventions with the ordinary malariasimulation builders and
re-run; everything is applied automatically from the parameter list,
with no extra arguments. A vivax list runs the same way:

``` r

p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
p <- malariasimulation::set_bednets(
  p, timesteps = 365, coverages = 0.6, retention = 3 * 365,
  dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
  rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)

out <- run_simulation_ode(timesteps = 3650, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))

# P. vivax, with radical cure: relapses and hypnozoite carriage come back as columns.
# Radical cure starts on day 1 here, from a seed without it, so these rows are
# still settling; burn in first for a level to compare.
pv <- malariasimulation::get_parameters(parasite = "vivax")
pv <- malariasimulation::set_drugs(pv, list(malariasimulation::CQ_PQ_params_vivax))
pv <- malariasimulation::set_clinical_treatment(pv, drug = 1, timesteps = 1, coverages = 0.4)
out_pv <- run_simulation_ode(timesteps = 3650, parameters = malariasimulation::set_equilibrium(pv, init_EIR = 3))
tail(out_pv[, c("n_relapses", "n_with_hypnozoites", "iaa_mean")])
```

Three exported functions:
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
runs the model,
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
holds the discretisation settings, and
[`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md)
gives the default graded age grid. Everything else comes off the
parameter list. (The `ode` in two of the names is kept for
compatibility: the model is a daily update, not an ODE.)

## How well does it match the IBM?

That question has its own project:
**[fleetcheck](https://pwinskill.github.io/fleetcheck/)**, a register of
claims about how closely `fleet` reproduces `malariasimulation`, each
with the criterion that decides it, the value measured against it, and a
verdict.

It lives outside this repository on purpose. Evidence kept beside the
code it vouches for is evidence the code’s author can quietly leave out
of date, so this README quotes no comparison figures of its own.

`fleetcheck` runs its own CI, fails when a verdict and the register
disagree, and refuses to tolerate a failing claim that has no written
reason.

## Documentation

|  |  |
|----|----|
| **[Get started](https://pwinskill.github.io/fleet/articles/fleet.html)** | A worked tour: run the model, read the outputs with `postie`, layer on each intervention, handle seasonality and burn-in. |
| **[Using fleet well](https://pwinskill.github.io/fleet/articles/using.html)** | Where the mean field departs from the IBM and what to do about it: things to do, things to leave alone, results to treat with caution, and what a run costs. |
| **[fleetcheck](https://pwinskill.github.io/fleetcheck/)** | The evidence, as a separate project: a register of claims about agreement with the IBM — core relationships, age structure, demography, interventions, and 63 country site files — each with its criterion, measurement and verdict. |
| **[Model specification](https://pwinskill.github.io/fleet/articles/model.html)** | The formal version: scope, what the state space does and does not carry, and the full daily update. |
| **[Parameter reference](https://pwinskill.github.io/fleet/articles/parameters.html)** | Every `malariasimulation` `set_*()` function argument by argument: what `fleet` reproduces exactly, what it approximates, and what it rejects. |

Function reference:
[`?run_simulation_ode`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.html),
[`?ode_tuning`](https://pwinskill.github.io/fleet/reference/ode_tuning.html),
[`?default_age_lower`](https://pwinskill.github.io/fleet/reference/default_age_lower.html).
Contributing:
[CONTRIBUTING.md](https://github.com/pwinskill/fleet/blob/main/CONTRIBUTING.md).

## License

MIT. Copyright (c) 2026 Peter Winskill. Full text:
[LICENSE.md](https://github.com/pwinskill/fleet/blob/main/LICENSE.md).
