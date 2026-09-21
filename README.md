# fleet <img src="man/figures/logo.png" align="right" width="30%" alt="fleet hex logo: a mosquito trailing light over a network of connected points" />

[![R-CMD-check](https://github.com/pwinskill/fleet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pwinskill/fleet/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/pwinskill/fleet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/pwinskill/fleet/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Status: work in progress](https://img.shields.io/badge/status-work%20in%20progress-orange.svg)](#not-ready-for-real-use)
[![Not for production use](https://img.shields.io/badge/not%20for-production%20use-red.svg)](#not-ready-for-real-use)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/pwinskill/fleet/blob/main/LICENSE.md)

> ## Not ready for real use
>
> ⚠️ **`fleet` is a work in progress and is not validated for research, policy or
> operational use.** It is published early so the approach and the comparison
> against `malariasimulation` can be looked at and argued with, not so that
> anyone can rely on its numbers.
>
> Concretely, and honestly:
>
> - **The API is unstable.** `run_simulation_ode()`'s signature has already
>   changed once and may change again without deprecation.
> - **Known discrepancies against the IBM are open, not resolved.** The age
>   distribution of clinical incidence does not match in every age band, and
>   across the 63-country site-file comparison `fleet` runs roughly 9% above the
>   IBM on clinical and severe incidence, an excess that is **not** explained.
>   The current verdict on every claim, with the numbers behind it, is in
>   *[fleetcheck](https://pwinskill.github.io/fleetcheck/)* — which is kept
>   current, unlike any figure quoted here would be.
> - **Severe incidence and anything derived from it (including DALYs) should be
>   treated as indicative only.**
> - Nothing here has been peer reviewed, and there is no versioned release.
>
> If you need results you can defend today, use
> [malariasimulation](https://github.com/mrc-ide/malariasimulation).

> A fast, deterministic **mean-field (ODE) twin** of the [malariasimulation](https://github.com/mrc-ide/malariasimulation) individual-based model of *Plasmodium falciparum* malaria: same inputs, seconds per run, population-independent.

## What it is

`fleet` reproduces the Griffin-style human model, structured by age and biting heterogeneity, coupled to the compartmental mosquito model. It is written in [odin2](https://github.com/mrc-ide/odin2) / [dust2](https://github.com/mrc-ide/dust2).

| | |
| --- | --- |
| **Human states** | `S / D / A / U / Tr`, plus `Ph` (post-treatment) and `Ph_c` (chemoprevention) prophylaxis |
| **Immunity** | four acquired states `IB / ICA / ID / IVA`; two maternal terms `ICM / IVM`, algebraic rather than state variables |
| **Mosquito, per species** | `E / L / P / Sm / EIP-chain / Im` |

It exists to give the malariasimulation ecosystem a **deterministic, Monte-Carlo-free companion** that:

- **Takes the same inputs as the individual-based model (IBM).** `run_simulation_ode()` accepts a `malariasimulation::get_parameters()` list, with the usual `set_*` intervention builders layered on, unchanged. *P. falciparum only.*
- **Is seeded at, and checked against, equilibrium.** Initial conditions come from `malariaEquilibrium`; with no interventions an undisturbed run relaxes off that seed by under half a percent over the first years and then holds, as the IBM does from the same seed.
- **Produces postie-compatible outputs.** The returned wide count table is malariasimulation-shaped, so a post-processing pipeline written for the IBM works on a `fleet` run unchanged, with no wrapper in between.
- **Is fast and population-independent.** All compartments are per-capita densities, so a 30-year daily run takes a few seconds whether you model a thousand people or ten million.

Reach for the IBM instead when you need stochastic variation, individual heterogeneity beyond the mean field, or *P. vivax*.

## Install

`fleet` compiles C++ at install time, as do several of its dependencies, so you need a working C++ toolchain first: Rtools on Windows, the Xcode command line tools on macOS, the usual build tools (`r-base-dev` or equivalent) on Linux.

```r
# install.packages("remotes")
remotes::install_github("pwinskill/fleet")
```

That also installs the GitHub-only hard dependencies (`dust2`, `monty` and `malariaEquilibrium`), which `DESCRIPTION` `Remotes` points at. It installs nothing from `Suggests`; the examples need two of those:

```r
remotes::install_github(c("mrc-ide/malariasimulation", "mrc-ide/postie"))
```

`malariasimulation` builds the parameter list and `postie` post-processes the output.

## Quick start

```r
library(fleet)

p <- malariasimulation::get_parameters()

# 10-year daily wide count table (malariasimulation-style columns)
out <- run_simulation_ode(timesteps = 3650, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))

# postie-format rates and prevalence, exactly as for an IBM run
postie::get_prevalence(out, diagnostic = "lm")$lm_prevalence_2_10
postie::get_rates(out)[, c("time", "age_lower", "age_upper", "clinical", "severe", "dalys")]
```

Layer interventions with the ordinary malariasimulation builders and re-run; everything is applied automatically from the parameter list, with no extra arguments:

```r
p <- malariasimulation::set_drugs(p, list(malariasimulation::AL_params))
p <- malariasimulation::set_clinical_treatment(p, drug = 1, timesteps = 1, coverages = 0.4)
p <- malariasimulation::set_bednets(
  p, timesteps = 365, coverages = 0.6, retention = 3 * 365,
  dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
  rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)

out <- run_simulation_ode(timesteps = 3650, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))
```

Three exported functions: `run_simulation_ode()` runs the model, `ode_tuning()` holds the solver and discretisation settings, and `default_age_lower()` gives the default graded age grid. Everything else comes off the parameter list.

## How well does it match the IBM?

That question has its own project: **[fleetcheck](https://pwinskill.github.io/fleetcheck/)**, a register of claims about how closely `fleet` reproduces `malariasimulation`, each with the criterion that decides it, the value measured against it, and a verdict.

It lives outside this repository on purpose. Evidence kept beside the code it vouches for is evidence the code's author can quietly leave out of date — and that is exactly what happened here: this README used to carry its own comparison figures, drawn from a copy of the scenario data that went eleven days and two model changes stale without anything noticing, still quoting bands built from ten IBM replicates after the count had been raised to twenty for being too narrow.

`fleetcheck` runs its own CI, fails when a verdict and the register disagree, and refuses to tolerate a failing claim that has no written reason.

## Documentation

| | |
| --- | --- |
| **[Get started](https://pwinskill.github.io/fleet/articles/fleet.html)** | A worked tour: run the model, read the outputs with `postie`, layer on each intervention, handle seasonality and burn-in. |
| **[Using fleet well](https://pwinskill.github.io/fleet/articles/using.html)** | Where the mean field departs from the IBM and what to do about it: things to do, things to leave alone, results to treat with caution, and what a run costs. |
| **[fleetcheck](https://pwinskill.github.io/fleetcheck/)** | The evidence, as a separate project: a register of claims about agreement with the IBM — core relationships, age structure, demography, interventions, and 63 country site files — each with its criterion, measurement and verdict. |
| **[Model specification](https://pwinskill.github.io/fleet/articles/model.html)** | The formal version: scope, what the state space does and does not carry, and the full ODE system. |
| **[Parameter reference](https://pwinskill.github.io/fleet/articles/parameters.html)** | Every `malariasimulation` `set_*()` function argument by argument: what `fleet` reproduces exactly, what it approximates, and what it rejects. |

Function reference: [`?run_simulation_ode`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.html), [`?ode_tuning`](https://pwinskill.github.io/fleet/reference/ode_tuning.html), [`?default_age_lower`](https://pwinskill.github.io/fleet/reference/default_age_lower.html). Contributing: [CONTRIBUTING.md](https://github.com/pwinskill/fleet/blob/main/CONTRIBUTING.md).

## License

MIT. Copyright (c) 2026 Peter Winskill. Full text: [LICENSE.md](https://github.com/pwinskill/fleet/blob/main/LICENSE.md).
