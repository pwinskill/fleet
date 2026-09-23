# Using fleet well: approximations and what to do differently

> **⚠️ Work in progress: not ready for real use.** `fleet` is published
> early so the approach, and its comparison against `malariasimulation`,
> can be examined and argued with, not so that anyone can rely on its
> numbers: the API is unstable and nothing here has been peer reviewed.
> The known discrepancies against the IBM are open rather than resolved,
> the largest being all-age severe incidence, which runs about 4 to 6%
> below the IBM, alongside a roughly 9% excess on clinical and severe
> incidence across the 63-country site files that is not yet explained.
> The full statement is on the [package front
> page](https://pwinskill.github.io/fleet/); if you need results you can
> defend today, use
> [malariasimulation](https://github.com/mrc-ide/malariasimulation).

`fleet` is a mean field. It carries the same mechanisms as the
individual-based model and, for most outcomes, lands inside its
replicate spread, but a stratum mean is not a population of individuals
and some things do not survive the averaging. What follows is where that
bites: what to do, what to leave alone, and which results to treat with
caution.

It assumes you can already run the model; if not, start with
[`vignette("fleet")`](https://pwinskill.github.io/fleet/articles/fleet.md).
Section references of the form §B.4 are to
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md),
which carries the formal system; function names link into
[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md);
the measurements behind the claims here are in
[fleetcheck](https://pwinskill.github.io/fleetcheck/), which is kept
current as the model changes, and this article is not.

## Things to do

**Call `set_equilibrium()` last.** It freezes all 56 equilibrium
parameters, 51 of them translated from the IBM list, at the values they
held when it ran, and `fleet` merges that stored set over the live list.
An edit made *after* it — another `set_*()` builder, or
`parameters$du <- 10` by hand — is a silent no-op, and the run comes out
bit-identical to one without the edit. The IBM reads those fields
directly and would honour the edit, so this is a real difference between
the two models, not a quirk of this one. `fleet` warns when a live value
disagrees with the stored one. If a calibration loop varies any of them,
re-call `set_equilibrium()` each iteration.

**Burn in before calibrating.** The seed is a close approximation, not
the exact fixed point (§F), and a seasonal run needs ~10 years to settle
onto its limit cycle. The seasonal annual-mean EIR sits a few percent
below the aseasonal `init_EIR` target through nonlinear averaging.
Discard the transient.

**Match the age grid to your demography.** If the oldest model age group
runs past the top `deathrate_agegroups`, `fleet` applies the top death
rate to those ages while the IBM removes them. Fix it either way round
(widen the grid, or extend the death rates):

``` r

run_simulation_ode(t, p, tuning = list(age_lower = default_age_lower(max_age = 100)))
```

A warning fires when they diverge; don’t ignore it, because it changes
the equilibrium age structure the whole run is seeded on
([`set_demography()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_demography)).

**Make your output bands a partition.** `postie` treats every column as
an independent age stratum, so overlapping bands are double-counted by
any person-day-weighted aggregate (§G). Set each family’s own fields:

``` r

p$prevalence_rendering_min_ages <- c(0, 2, 10) * 365
p$prevalence_rendering_max_ages <- c(2, 10, 100) * 365
```

**Loosen the solver for long projections.** The equilibrium-preserving
defaults are slow and only needed near equilibrium:

``` r

run_simulation_ode(t, p, tuning = list(rtol = 1e-6, step_size_max = 10))
```

**Convert EIR before comparing it.** `fleet` reports `EIR` per adult per
year. The IBM reports a per-species daily total, so the like-for-like
comparison is `EIR_<species> / human_population * 365`.

**Sanity-check conservation.**
`S_count + D_count + A_count + U_count + Tr_count + Ph_count` should
equal `human_population` to numerical tolerance.

## Things to leave alone

**`acquired_immunity_offset`.** Leave it at 0. Setting 0.5 reproduces
the IBM’s literal per-individual offset but is a *worse* match to the
IBM ensemble mean (§B.3).

**`bite_dedup`.** Leave it at 1 for any scientific run. Set 0 only for
an equilibrium test, where it cuts the drift off the seed to well under
1% (§B.2).

**Erlang stage counts (`n_eir`, `n_foim`, `n_eip`).** The defaults are
fine. Raising them sharpens the gamma-shaped lags toward the IBM’s fixed
delays at linear cost in state size; equilibrium is exact for any value.
The prophylaxis chains (`n_ph`, `n_phc`) are sized from the drug’s
Weibull by default (§B.4) and capped at 20; each stage adds $`n_a n_z`$
states and the solver slows steeply beyond that.

**`parameterise_total_M()` / `parameterise_mosquito_equilibrium()`.**
They have no effect; `fleet` re-derives `total_M` itself
([`set_equilibrium()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_equilibrium):
from `init_EIR` under the default demography, as the IBM’s own value
under a custom one).

## Results to treat with caution

### Severe incidence

`fleet` runs below the IBM here, and DALYs derived from it inherit that.
The size of the gap is not quoted in this article on purpose: it moves
whenever the model does, and a number written here would go stale
silently. The measurement that is kept current is the
`severe-allage-eir` claim in
[fleetcheck](https://pwinskill.github.io/fleetcheck/).

The obvious explanation is the wrong one: $`\theta`$ is the *shallowest*
of the three acquired-immunity Hill functions, not the steepest. What
singles it out is position. $`I_{V0}`$ is sixteen times smaller than
$`I_{C0}`$ while the two acquired immunities take nearly the same values
at every age, so by age 3 $`\theta`$ is already out on its convex tail
while $`\phi`$ is still nearly linear. `fleet` evaluates $`\theta`$ once
per stratum, at that stratum’s mean immunity; the IBM averages it over
individuals whose bite histories differ *within* the stratum. Averaging
a convex function the second way gives the larger answer, so the mean
field returns the smaller one and severe incidence comes out low. §B.3
derives this and gives the numbers.

Two things make severe incidence uniquely exposed. $`\theta`$ is a
**post-hoc multiplier** — severe counts are $`\theta`$ times an
infection count the state equations have already determined (§G), so
nothing downstream absorbs the error, where $`\phi`$ splits the
infection inflow and feeds back through infectivity, treatment and
immunity. And **the bias needs the tail**: at EIR 1 $`\theta`$ sits near
its plateau and all-age severe incidence is within 1% of the IBM, the
gap opening only over the EIR 20–120 range where the tail is where the
population is.

None of this is a correction you can apply. It sets the sign and the
order of magnitude, and it says which comparisons to distrust first:
severe incidence in narrow age bands, and any scenario difference taken
on severe incidence.

### Vector-control eras

Population-averaged coverage loses the correlation of protection within
individuals across bites, so `fleet` tends to sit above the IBM while
nets or IRS are active (§2.2).

### Repeated mass campaigns

Overlapping mass PEV or TBV campaigns combine as independent protections
rather than most-recent-receipt, and `min_wait` re-vaccination exclusion
is not applied
([`set_mass_pev()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mass_pev),
[`set_tbv()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_tbv)).
A single campaign is fine.

### Stochastic scatter

`fleet` returns the ensemble mean, so it has no replicate spread of its
own. A single IBM replicate is not a fair yardstick for it: compare
against the IBM’s replicate median and read its 10–90% band as the noise
`fleet` is not carrying.

## Cost model

Seconds for one complete
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
call: building the inputs, seeding the equilibrium, integrating and
rendering the outputs, which is what you actually pay. Single core, at
`tuning = list(rtol = 1e-6)`.

| Scenario                     | 5 years | 10 years | 30 years |
|------------------------------|---------|----------|----------|
| No interventions             | 0.5 s   | 1.0 s    | 3.1 s    |
| Seasonal                     | 0.9 s   | 1.7 s    | 5.5 s    |
| Seasonal + treatment (AL)    | 1.0 s   | 1.9 s    | 5.8 s    |
| Seasonal + nets + IRS        | 1.6 s   | 1.8 s    | 4.8 s    |
| Seasonal + SMC (4 rounds/yr) | 1.4 s   | 2.5 s    | 7.8 s    |
| All of the above + RTS,S     | 1.5 s   | 3.0 s    | 11.2 s   |

Three things to read off it.

**Seasonality adds about 80%** (1.8x, 1.7x and 1.8x on the three
horizons), because rainfall drives the larval carrying capacity on a
daily interpolation grid, and the oscillating solution costs the
adaptive stepper 4.7 steps per output day against 2.0 on an aseasonal
run.

**Chemoprevention is the expensive intervention**, and not because of
the compartments it adds: each round splits the integration into a fresh
segment (§E), so a schedule with many rounds pays for many solver
restarts. A 20-year monthly PMC schedule means ~240 of them. Prefer a
coarser cadence when the extra resolution buys nothing.

**Cost is close to proportional to the horizon.** Fit the
no-intervention row and you get 0.102 s per simulated year on an
intercept of -0.01 s, so the equilibrium solve is not a meaningful fixed
charge. Three rows are dearer per simulated year at 5 years than at 30
(seasonal + nets + IRS by a wide margin, 0.32 s/y against 0.16 s/y; then
seasonal + SMC and seasonal + treatment by a little). What is fixed in
those rows is building the intervention time series, not the seed.

Population size is irrelevant, as it should be: the compartments are
per-capita densities, and `human_population` only rescales the count
columns on the way out.

| Population           | 1,000 | 10,000 | 100,000 | 1,000,000 | 10,000,000 |
|----------------------|-------|--------|---------|-----------|------------|
| 30-year seasonal run | 5.6 s | 5.8 s  | 5.6 s   | 5.6 s     | 5.6 s      |

For long projections the solver settings are worth knowing. The same
30-year seasonal run takes **7.9 s** at the
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
defaults and **5.5 s** at `rtol = 1e-6`. It is the relative tolerance
that buys this, by roughly halving the step count; `step_size_max` is a
safety rail rather than a speed control, and changes neither the step
count nor the outputs. Loosening `atol` as well saves almost nothing
(5.4 s) and is not worth having: the prophylaxis chain stages hold
occupancies of order `1e-6`, which a looser absolute tolerance lets go
slightly negative.

Measured as the minimum of five repeats, since contention can only add
time, on R 4.5.2, `aarch64-w64-mingw32`. Treat them as indicative: one
machine, one core, and a laptop under load will be slower. For how this
compares with the IBM’s cost, see the `speed` claim in
[fleetcheck](https://pwinskill.github.io/fleetcheck/).

## Index of approximations

Every 🟡 in
[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md),
and when it matters. The first column links to the argument-by-argument
account; § refers to
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md).

| Mechanism | Where | Matters when |
|----|----|----|
| Ages above the top `deathrate_agegroups` take the top death rate | [`set_demography()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_demography) | A custom demography whose age grid stops short of the model’s oldest group; the IBM removes those people instead. Warned |
| Weibull prophylaxis → Erlang chain (moment-matched) | [`set_drugs()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_drugs), §B.4 | The far tail of chemoprevention protection; for post-treatment protection the exponential $`T`$ stage in front of the chain spreads AL’s sharp 10-day protection, reproducing its integral rather than its shape |
| Multi-drug blend rather than parallel sub-populations | [`set_clinical_treatment()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_clinical_treatment) | A genuinely mixed first line (a *switch* is modelled fine) |
| Single scalar slow-clearance rate blended across drugs | [`set_antimalarial_resistance()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_antimalarial_resistance), §B.4 | High artemisinin resistance with several treatment drugs |
| Population-averaged nets | [`set_bednets()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_bednets) | Any nets era, the largest routine gap |
| Population-averaged IRS | [`set_spraying()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_spraying) | Any IRS era |
| Carrying-capacity floor at $`K_0\times10^{-4}`$ | [`set_carrying_capacity()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_carrying_capacity) | Modelling complete vector elimination |
| Chemoprevention as pulses; half-open age bands | [`set_mda()`, `set_smc()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mda), §E | Very narrow target bands |
| PMC age trigger → monthly pulses | [`set_pmc()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_pmc) | PMC dose timing specifically |
| Seasonal PEV boosters → fixed delay | [`set_pev_epi()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_pev_epi) | Seasonally-timed booster campaigns (warned) |
| Repeated campaigns: no `min_wait`, independent combination, no ageing out of band | [`set_mass_pev()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mass_pev) | Repeated mass PEV |
| TBV campaign combination | [`set_tbv()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_tbv) | Repeated TBV campaigns |
