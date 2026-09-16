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
and some things do not survive the averaging. This article is the
practical consequence of that: where `fleet` departs from
`malariasimulation`, what to do about it, and which results to treat
with caution.

It assumes you can already run the model. If not, start with
[`vignette("fleet")`](https://pwinskill.github.io/fleet/articles/fleet.md).
For the formal system, the state space and an argument-by-argument
account of every `set_*()` function, see
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md);
for the measured agreement behind the claims here, see
[`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md).

Section references below of the form §B.4 or §3.15 are to that
specification article,
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md),
not to this one.

## Things to do

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
the equilibrium age structure the whole run is seeded on (§3.4).

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
They have no effect; `fleet` re-derives `total_M` itself (§3.3: from
`init_EIR` under the default demography, as the IBM’s own value under a
custom one).

## Results to treat with caution

**Severe incidence.** `fleet` runs below the IBM here, by 6.2% on all
ages at EIR 20, 4.0% at EIR 50 and 4.6% at EIR 120, and by 15–25% in the
5–20 year bands; the DALYs derived from severe incidence inherit that.
It is the largest open discrepancy between the two models
([`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)),
so it is worth being exact about the cause, because the obvious
explanation is the wrong one.

It is *not* that $`\theta`$ is the steepest of the three
acquired-immunity Hill functions. It is the shallowest. From
`get_parameters()`:

|            | exponent        | scale              |
|------------|-----------------|--------------------|
| $`b`$      | $`k_b = 2.160`$ | $`I_{B0} = 43.9`$  |
| $`\phi`$   | $`k_c = 2.369`$ | $`I_{C0} = 18.02`$ |
| $`\theta`$ | $`k_v = 2.000`$ | $`I_{V0} = 1.096`$ |

What singles $`\theta`$ out is the **scale**, not the exponent. A Hill
function is flat well below its scale parameter and only approaches its
full log-log slope well above it, so its local slope is a matter of
*position*. Writing $`X =
f_{v,i}\left(\left(I_{VA}+\delta+I_{VM}\right)/I_{V0}\right)^{k_v}`$ for
the denominator term of $`\theta`$ (§B),

``` math
\frac{\partial\log\theta}{\partial\log I_{VA}} =
-\,\frac{k_v\left(1-\theta_1\right)X}{\left(1+X\right)\left(1+\theta_1X\right)},
```

which is near 0 for $`X \ll 1`$, reaches $`-1.96`$ at
$`X = \theta_1^{-1/2} \approx
92`$, and returns toward 0 once $`\theta`$ has bottomed out on its
$`\theta_0\theta_1`$ floor. The same expression with $`k_c`$ and
$`\phi_1`$ governs $`\phi`$, whose steepest attainable slope, $`-2.24`$,
is the steeper of the two.

Position is what separates them, and $`I_{V0}`$ is 16 times smaller than
$`I_{C0}`$ while the two acquired immunities take nearly the same values
at every age. At `fleet`’s own seed at the reference EIR of 20,
$`I_{VA} = 9.0`$ and $`I_{CA} = 9.5`$ at age 3, and 37.6 against 40.3 at
age 10: they differ almost entirely in what they are divided by. So
$`I_{VA}`$ is $`8\,I_{V0}`$ at age 3 and $`34\,I_{V0}`$ at age 10,
putting $`\partial\log\theta/\partial\log I_{VA}`$ at $`-1.87`$ at age 3
and holding it between $`-1.96`$ and $`-1.80`$ from age 4 to age 10,
within a few percent of the steepest slope the function can have. Over
the same ages $`I_{CA}`$ climbs only from $`0.5\,I_{C0}`$ to
$`2.2\,I_{C0}`$ and $`\phi`$’s slope from $`-0.42`$ to $`-2.05`$: at age
3, $`\phi`$ is still nearly linear, less than a quarter as steep as
$`\theta`$, and it does not overtake $`\theta`$ until about age 8.

Sitting on that tail makes $`\theta`$**convex** in $`I_{VA}`$ across
essentially the whole population, where $`\phi`$ is concave in
$`I_{CA}`$ below age 4 at this EIR and only convex above it. That fixes
the sign of the mean-field error. `fleet` evaluates $`\theta`$ once per
(age group $`\times`$ heterogeneity group) at that stratum’s mean
$`I_{VA}`$; the IBM averages $`\theta`$ over individuals whose bite
histories differ *inside* the stratum. By Jensen’s inequality
$`\mathbb{E}\left[\theta(I_{VA})\right] \geq \theta\left(\mathbb{E}[I_{VA}]\right)`$,
so the mean field returns the smaller number and severe incidence comes
out low. The magnitude scales with the residual within-stratum variance:
a coefficient of variation of 0.3 in $`I_{VA}`$ inside a stratum puts
$`\mathbb{E}[\theta]`$ 23% above $`\theta`$ at the stratum mean at age 3
and 28% above at age 5, where the same calculation on $`\phi`$ gives
$`-0.7`$% and $`+4.8`$%. `fleet` resolves the log-normal biting strata
explicitly (§2.1), so only the spread within an (age, $`\zeta`$) cell is
lost and the realised gap is a few percent rather than tens of percent,
but the direction is the same.

Two things make severe incidence uniquely exposed to this:

- **$`\theta`$ is a post-hoc multiplier.** `n_inc_severe_*` is
  $`\theta`$ times an infection count the state equations have already
  determined (§G); nothing downstream responds, so the whole averaging
  error lands in the output. $`\phi`$ is not like that. It splits the
  infection inflow between $`D`$ and $`T`$, which changes infectivity,
  treatment and the immunity trajectories themselves, so the equilibrium
  absorbs part of any error in it.
- **The bias needs the tail.** At EIR 1, $`I_{VA}`$ is only
  $`0.46\,I_{V0}`$ at age 3, $`\theta`$ sits near its $`\theta_0`$
  plateau (log-slope $`-0.09`$) and there is next to nothing to average
  badly: all-age severe incidence is within 1% of the IBM at EIR 1, 3
  and 10. The gap opens over the EIR 20–120 plateau, where the tail is
  where the population is.

None of this is a correction you can apply. It sets the sign and the
order of magnitude, and it says which comparisons to distrust first:
severe incidence in narrow age bands, and any scenario difference taken
on severe incidence.

**Anything through a nets or IRS era.** Population-averaged coverage
loses the correlation of protection within individuals across bites, so
the ODE tends to sit above the IBM while vector control is active
(§2.2).

**Repeated mass PEV or TBV campaigns.** Overlapping campaigns combine as
independent protections rather than most-recent-receipt, and `min_wait`
re-vaccination exclusion is not applied (§3.15, §3.16). A single
campaign is fine.

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

Measured with `comparison/benchmark.R` (minimum of five repeats, since
contention can only add time) on R 4.5.2, `aarch64-w64-mingw32`. Treat
them as indicative: one machine, one core, and a laptop under load will
be slower. For how this compares with the IBM’s cost, see
[`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md).

## Index of approximations

Every 🟡 in §3, and when it matters:

| Mechanism | § | Matters when |
|----|----|----|
| Weibull prophylaxis → Erlang chain (moment-matched) | 3.5, B.4 | The far tail of chemoprevention protection; for post-treatment protection the exponential $`T`$ stage in front of the chain spreads AL’s sharp 10-day protection, reproducing its integral rather than its shape |
| Multi-drug blend rather than parallel sub-populations | 3.6 | A genuinely mixed first line (a *switch* is modelled fine) |
| Single scalar slow-clearance rate blended across drugs | 3.7, B.4 | High artemisinin resistance with several treatment drugs |
| Population-averaged nets | 3.8 | Any nets era, the largest routine gap |
| Population-averaged IRS | 3.9 | Any IRS era |
| Carrying-capacity floor at $`K_0\times10^{-4}`$ | 3.10 | Modelling complete vector elimination |
| Chemoprevention as pulses; half-open age bands | 3.11, E | Very narrow target bands |
| PMC age trigger → monthly pulses | 3.12 | PMC dose timing specifically |
| Seasonal PEV boosters → fixed delay | 3.14 | Seasonally-timed booster campaigns (warned) |
| Repeated campaigns: no `min_wait`, independent combination, no ageing out of band | 3.15 | Repeated mass PEV |
| TBV campaign combination | 3.16 | Repeated TBV campaigns |
