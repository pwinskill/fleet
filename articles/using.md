# Using fleet well: approximations and what to do differently

> **⚠️ Work in progress: not ready for real use.** `fleet` is published
> early so the approach, and its comparison against `malariasimulation`,
> can be examined and argued with, not so that anyone can rely on its
> numbers: the API is unstable and nothing here has been peer reviewed.
> The known discrepancies against the IBM are open rather than resolved:
> an excess on falciparum clinical and severe incidence across the
> 63-country site files that is not yet explained; a larger excess on
> *P. vivax* sites, not yet explained either; and, at school age, vivax
> LM prevalence and clinical incidence a few per cent high. The full
> statement is on the [package front
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
the exact fixed point (§F): prevalence relaxes by under half a per cent,
but incidence by a few per cent over the first years, more at high
transmission (+6% for under-5 clinical incidence at EIR 120). A vivax
run moves much further, as the IBM’s does, because both build up a
within-cell spread of immunity the seed does not have: its clinical
incidence settles about 20% above the seed and its realised EIR 4 to 17%
above `init_EIR` (more at lower EIR), within about eight years at EIR 3
but over two decades at EIR 0.3, so `init_EIR` is not the EIR a vivax
run realises. A treated run whose treatment starts at timestep 1 starts
untreated, as the IBM’s does, and settles onto its treated state over
the first months. A seasonal run needs ~10 years to settle onto its
limit cycle, and its annual-mean EIR sits about 7% below the aseasonal
`init_EIR` target through nonlinear averaging. Discard the transient.

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
any person-day-weighted aggregate (§G). A band `[lower, upper]` holds
both its end days, as in the IBM, so a partition ends each band the day
before the next one starts. Set each family’s own fields:

``` r

p$prevalence_rendering_min_ages <- c(0, 2, 10) * 365
p$prevalence_rendering_max_ages <- c(2, 10, 100) * 365 - 1
```

A family whose list is empty renders nothing, in `fleet` as in the IBM:
the defaults give the 2-10 prevalence band and no incidence at all.

**Compare rows by timestep.** Row $`t`$ of a run is day $`t`$ — the
state at its start and the incidence during it — exactly as
`malariasimulation` renders its timestep $`t`$, so the two tables line
up row for row.

**Read prevalence as a count over a population.** `p_detect_lm_*` is the
IBM’s expected number LM-detected, a count like `n_detect_lm_*`, not a
prevalence. Prevalence is `n_detect_lm_*` / `n_age_*`, as `postie`
computes it.

**Convert EIR before comparing it.** `fleet` reports `EIR` per adult per
year, the value biting humans today. The IBM reports the same lagged
quantity as a per-species daily total, so the like-for-like comparison
is `EIR_<species> / human_population * 365`.

**Sanity-check conservation.**
`S_count + D_count + A_count + U_count + Tr_count` should equal
`human_population` to numerical tolerance. `Ph_count` is the
drug-protected part of `S_count`, so it is not added again.

## Things to leave alone

**`acquired_immunity_offset`.** Leave it unset. It defaults to 0 for
falciparum, whose curves are read at a stratum mean, where the IBM’s
per-individual 0.5 is a *worse* match to the IBM ensemble mean (§B.3);
to 0.5 for vivax, whose curves are read at quadrature nodes that stand
for people (§V.3); and to 0 for vivax with `immunity_spread = FALSE`,
back at the cell mean.

**`bite_dedup`.** Leave it at 1 for any scientific falciparum run.
Setting 0 lets every bite infect independently, which the IBM does not
do for falciparum (§B.2). Vivax counts every bite, as the IBM does, and
ignores it.

**Mosquito sub-steps (`n_sub`).** The default of 32 integrates each day
about as accurately as the IBM’s own mosquito solver does, and the lags
the model carries are the IBM’s own delay lines, so there is nothing to
sharpen. The prophylaxis chains (`n_ph`, `n_phc`) are sized from the
drug’s Weibull by default (§B.4) and capped at 20; each stage adds
$`n_a n_z`$ states, and no chain can have more stages than its mean has
days.

**`parameterise_total_M()` / `parameterise_mosquito_equilibrium()`.**
They have no effect; `fleet` re-derives `total_M` itself
([`set_equilibrium()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_equilibrium):
from `init_EIR` under the default demography, as the IBM’s own value
under a custom one).

## Results to treat with caution

### Severe incidence, and the age grid

Severe incidence is where the age grid shows most. $`\theta`$ sits out
on the convex tail of its Hill function from early childhood (§B.3 has
the numbers), and a group’s severe incidence is evaluated at its mean
immunity, so wherever immunity changes across a group, as it does
fastest in the young groups where severe disease concentrates, the group
returns too small a number. The error is first order in the group width.
On the default 118 groups severe incidence at EIR 120 is 2% below the
IBM median; on 53 groups, under half the cost, it runs 7% below, and
clinical incidence in young children a few per cent low. So a result
that rests on the level of severe incidence, or of clinical incidence in
under-5s at high transmission, wants the default grid or a finer one,
and a coarser grid is for everything else:

``` r

run_simulation_ode(t, p, tuning = list(age_lower = default_age_lower(n_group = 53)))
```

The same convexity inside an (age, $`\zeta`$) stratum, where the IBM
holds a spread of infection histories and the mean field one value, is a
much smaller effect for falciparum. The measurements, kept current as
the model changes, are the `severe-allage-eir` and
`age-profile-clinical` claims in
[fleetcheck](https://pwinskill.github.io/fleetcheck/).

$`\theta`$ is also a **post-hoc multiplier** — severe counts are
$`\theta`$ times an infection count the state equations have already
determined (§G) — so nothing downstream absorbs an error in it, where
$`\phi`$ splits the infection inflow and feeds back through infectivity,
treatment and immunity.

### *P. vivax* where biteless days are common

`malariasimulation` 3.0.0 evaluates vivax relapses only on days when at
least one person is bitten, so wherever days without a bite are common
the IBM loses relapses that `fleet` keeps, and runs below it, sometimes
to elimination. That is a defect in the IBM, which `fleet` does not
copy; comparing the two there compares `fleet` with the defect. The
chance of a biteless day in the IBM is about
$`e^{-\text{EIR}_{\text{day}}\, \bar\psi\, N}`$ for a daily EIR per
adult, mean relative biting rate $`\bar\psi`$ and population $`N`$: in a
10,000-person run it passes 5% only below an EIR of about 0.14, and in
50,000 below about 0.03. So it matters in small populations, at the
lowest transmission, and where vector control drives the day’s bites
down. It does not explain `fleet`’s excess over the IBM on the vivax
site files, whose 50,000-person runs almost never go a day without a
bite; that excess is not yet explained.

### *P. vivax* at school age

The vivax immunity curves are steep, and the IBM’s immunity-dependent
transitions sort people within an age, heterogeneity and batch cell:
those with more anti-parasite immunity clear sub-patent infections
sooner, so the uninfected hold more of the immune, and the clinical draw
takes the least immune. A mean field with one immunity distribution per
cell does not carry that sorting (§V.7), and at EIR 10 its LM prevalence
and clinical incidence in school-age children run 1–3% above the IBM’s.

### Vector-control eras

Population-averaged coverage loses the correlation of protection within
individuals across bites, so `fleet` tends to sit above the IBM while
nets or IRS are active (§2.2).

### Repeated mass campaigns

A campaign’s protection ages with the cohort it vaccinated, and repeated
rounds combine as the most recent dose, as in the IBM. What the mean
field cannot carry is who was reached before: each round reaches a
person independently of the last, and `min_wait` re-vaccination
exclusion is not applied
([`set_mass_pev()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mass_pev),
[`set_tbv()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_tbv)).
A single campaign is unaffected.

### Stochastic scatter

`fleet` returns the ensemble mean, so it has no replicate spread of its
own. A single IBM replicate is not a fair yardstick for it: compare
against the IBM’s replicate median, and read the replicate band (median
± 1.28 SD, `fleetcheck`’s definition) as the noise `fleet` is not
carrying.

## Cost model

Seconds for one complete
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
call: building the inputs, seeding the equilibrium, stepping the days
and rendering the outputs, which is what you actually pay. Single core,
at the
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
defaults, on the default 118-group age grid. Falciparum at EIR 20, the
SMC scenarios at 15:

| Scenario                     | 5 years | 10 years | 30 years |
|------------------------------|---------|----------|----------|
| No interventions             | 0.30 s  | 0.59 s   | 1.65 s   |
| Seasonal                     | 0.31 s  | 0.60 s   | 1.67 s   |
| Seasonal + treatment (AL)    | 0.33 s  | 0.62 s   | 1.75 s   |
| Seasonal + nets + IRS        | 0.31 s  | 0.61 s   | 1.79 s   |
| Seasonal + SMC (4 rounds/yr) | 0.42 s  | 0.75 s   | 2.15 s   |
| All of the above + RTS,S     | 0.43 s  | 0.81 s   | 2.34 s   |

Three things to read off it.

**Cost is proportional to the horizon.** Every day costs the same, so
the no-intervention row is 0.05 s per simulated year on an intercept of
a few hundredths of a second: building the inputs and solving the seed
are a small fixed charge next to the days.

**Seasonality costs almost nothing**, 1% at 30 years. The update does
the same work each day however fast the solution is moving; rainfall
only changes the day’s carrying capacity.

**Chemoprevention is the dearer intervention**, SMC 29% over a seasonal
run, and everything together, vaccination included, 40%. What SMC adds
is state – the two chemoprevention protection chains and the treated
phase, about six hundred states per stage on the default grid – and each
round pauses the run for its pulse.

A vivax run costs about sixteen times a falciparum one: every
compartment carries the hypnozoite-batch dimension, eleven levels, and
each cell carries its immunity’s spread. At EIR 3, with 40% of clinical
cases treated where a drug is named:

| Scenario                        | 5 years | 10 years | 30 years |
|---------------------------------|---------|----------|----------|
| No interventions                | 4.7 s   | 8.8 s    | 25.8 s   |
| Treatment (chloroquine)         | 5.5 s   | 9.9 s    | 27.2 s   |
| Radical cure (CQ + primaquine)  | 10.9 s  | 20.4 s   | 58.5 s   |
| Radical cure (CQ + tafenoquine) | 19.8 s  | 38.4 s   | 107 s    |

Radical cure is the dear vivax intervention: it adds the
liver-stage-protected levels to the batch dimension, four for primaquine
and eleven for tafenoquine, and the day’s radical-cure transfers, 2.1
and 3.9 times the cost of treatment alone at 30 years.

Population size is irrelevant, as it should be: the compartments are
per-capita densities, and `human_population` only rescales the count
columns on the way out.

| Population           | 1,000  | 10,000 | 100,000 | 1,000,000 | 10,000,000 |
|----------------------|--------|--------|---------|-----------|------------|
| 30-year seasonal run | 1.74 s | 1.73 s | 1.70 s  | 1.72 s    | 1.70 s     |

Of the discretisation settings, the mosquito sub-steps barely move the
cost, since the mosquito model is a small part of the day: a 30-year
seasonal run takes 1.70 s at `n_sub = 8`, 1.72 s at the default 32 and
1.75 s at 64. The age grid is the setting that costs, a little more than
in proportion to its number of groups: on the 53 groups of
`default_age_lower(n_group = 53)` the same run takes 0.70 s.

Measured as the minimum of five repeats, two for vivax, since contention
can only add time, on R 4.5.2, `aarch64-w64-mingw32`. Treat them as
indicative: one machine, one core, and a laptop under load will be
slower. The 30-year seasonal run appears three times above, at 1.67 to
1.75 s; that spread is the machine’s, not the model’s. For how this
compares with the IBM’s cost, see the `speed` claims in
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
| Weibull prophylaxis → stage chain (moment-matched), protecting fully while it lasts | [`set_drugs()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_drugs), §B.4 | SMC and MDA at high transmission, where the chain’s all-or-nothing protection shields dosed children a little more than the IBM’s partial, fading protection does; the far tail of chemoprevention protection; for post-treatment protection the $`T`$ stage in front of the chain spreads AL’s sharp 10-day protection, reproducing its integral rather than its shape |
| Multi-drug blend rather than parallel sub-populations, weighted by coverage where the IBM’s treated are weighted by success | [`set_clinical_treatment()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_clinical_treatment) | A genuinely mixed first line, especially with resistance on one drug (a *switch* is modelled fine) |
| Single scalar slow-clearance rate blended across drugs | [`set_antimalarial_resistance()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_antimalarial_resistance), §B.4 | High artemisinin resistance with several treatment drugs |
| Population-averaged nets | [`set_bednets()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_bednets) | Any nets era, the largest routine gap |
| Population-averaged IRS | [`set_spraying()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_spraying) | Any IRS era |
| Carrying-capacity floor at $`K_0\times10^{-4}`$ | [`set_carrying_capacity()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_carrying_capacity) | Modelling complete vector elimination |
| Chemoprevention as pulses; target bands half-open, where the IBM’s hold both end days | [`set_mda()`, `set_smc()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mda), §E | Very narrow target bands |
| PMC age trigger → monthly pulses | [`set_pmc()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_pmc) | PMC dose timing specifically |
| Seasonal PEV boosters → fixed delay | [`set_pev_epi()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_pev_epi) | Seasonally-timed booster campaigns (warned) |
| Repeated campaigns: no `min_wait`, each round reaching people independently | [`set_mass_pev()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mass_pev), [`set_tbv()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_tbv) | Repeated mass PEV or TBV rounds |
| The absorbing top age group taken whole, or not at all, by a rendering band | [`set_epi_outputs()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_epi_outputs), §G | A band ending inside the top group (above 80 years on the default grid), which is given all of it |
| Refractory boosting of $`I_{CA}`$, $`I_D`$, $`I_{VA}`$ as a quasi-steady closure | §B.5 | Exposure changing sharply within a week |
| One immunity value per age group | [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md), §B.3 | Severe incidence, and clinical incidence in young children at high transmission, on a grid coarser than the default |
| Vivax: one immunity distribution per age, heterogeneity and batch cell, shared by every disease state | §V.7 | School-age LM prevalence and clinical incidence, 1–3% high at EIR 10 |
| Vivax: IAA’s 44-day refractory window as four 11-day stages | §V.3 | Infants’ IAA, 2.7% of itself; nothing else by more than 0.7% |
| Vivax: prophylaxis chains protecting fully while they last, Tr included | [`set_drugs()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_drugs), §V.4–V.5 | The shape of liver-stage protection after radical cure; under a switch from primaquine to tafenoquine, tafenoquine’s protection keeps its mean and widens its spread |
| Vivax: the day’s bites averaged before the competing draw | §V.2, §V.7 | D’s infection probability 0.1–0.5% low, and the bite share of infections 0.2–1% low |
| Vivax: refractory people infected at their cell’s average rate, Tr included | §V.3, §V.7 | IAA boosts in young children under treatment, up to an estimated 1.6% too few |
| Vivax: IAA’s within-cell spread from a renewal closure | §V.3, §V.7 | Clinical incidence and LM prevalence where daily infection probabilities pass about 0.05, at high EIR |
