# fleet 0.0.0.9003

`fleet` now advances one day at a time, in the order `malariasimulation`
resolves a day, instead of integrating an ODE; it runs *P. vivax* as well as
*P. falciparum*; and its default age grid is four times finer. Four falciparum
mechanisms were brought into line with the IBM along the way, and the output
table is now `malariasimulation`'s, column for column. Every model output moves.

## *P. vivax*

* **`run_simulation_ode()` runs *P. vivax*.** A list from
  `malariasimulation::get_parameters(parasite = "vivax")` runs a vivax block in
  the same daily model, dispatched on `parameters$parasite` as
  `malariasimulation` does. It carries the hypnozoite-batch dimension (0 to
  `kmax` batches, relapse at `k * f`, one batch lost a day with probability
  `1 - exp(-k * gammal)`), radical cure with liver-stage protection
  (`CQ_PQ_params_vivax`, `CQ_TQ_params_vivax`), and no severe disease. MDA, SMC
  and PMC are refused under vivax, because `malariasimulation` itself fails on
  them. `vignette("model")` §V specifies it.
* Each day is the IBM's: every bite counts, relapse adds to the bite hazard and
  PEV reduces both, and infection competes with each state's progression in
  one draw a day. People who share an age, heterogeneity group and batch count
  differ in immunity by their batch history, and the vivax curves are steep, so
  each cell carries the second moment of immunity and the curves are averaged
  over a gamma (`parameters$immunity_spread`, default `TRUE`).
* **Refractory windows are carried as stocks, moved with the people.** A bite
  boosts a person and moves them up a batch, into a cell whose relapse hazard
  is higher, so a cell's newcomers are freshly boosted; a renewal rate
  `p / (p u + 1)` cannot see that, and put young children's anti-parasite
  immunity 10-11% above the IBM's. Against four 50,000-person IBM runs at EIR
  10, IAA by age now agrees but for the age grid's share, and so does its
  within-cell spread.
* **The vivax clinical curve takes the IBM's per-person +0.5**, at quadrature
  nodes that stand for people: `acquired_immunity_offset` defaults to 0.5 for
  vivax, 0 for falciparum.
* **A bite that forms a batch keeps the day's batch decay from happening.** The
  IBM queues the decay first and the bite's new batch count later, and the later
  update overwrites the earlier. The decay is taken after the day's other
  events, from all but the day's batch-formers, so someone who relapses and
  loses a batch on the same day ends a batch down, infected, as in the IBM.
* **The radically cured carry the immunity of the nodes that fell ill**, the
  less immune, rather than their cell's mean, and ICA's within-cell spread
  gains the square of the day's boost probability, which makes it exact with no
  refractory window.
* **With heterogeneity off, vivax is seeded as the IBM seeds it**: everyone
  from the heterogeneous equilibrium's first group, among mosquitoes sized from
  the whole of it.
* Against the same IBM runs, infections, relapses, clinical incidence, hypnozoite
  carriage, the batch count, and LM and PCR prevalence agree by age to within
  about 3%. School-age LM prevalence and clinical incidence run 1-3% high: the
  IBM's immunity-dependent transitions sort people by immunity within a cell,
  which one immunity distribution per cell does not carry.
* A vivax run costs about 0.9 s per simulated year on the default grid, 2 s
  with primaquine radical cure and 3.6 s with tafenoquine: the hypnozoite
  dimension and the immunity spread make it about sixteen times a falciparum
  run.
* The output is the IBM's vivax table: `n_relapses`, `n_with_hypnozoites`,
  `n_inc_relapse_*` and `n_with_hypnozoites_*` over their rendering lists, no
  severe columns unless a severe band is set, and no `p_detect_lm_*`.

## Output

* **The immunity means are rendered**, as the IBM renders them: `ica_mean`,
  `icm_mean`, `ib_mean`, `iva_mean`, `ivm_mean`, `id_mean` for falciparum,
  `ica_mean`, `icm_mean`, `iaa_mean`, `iam_mean`, `hypnozoites_mean` for vivax,
  and each by age over its `<name>_rendering_*` list.

## The default age grid

* **`default_age_lower()` now has 118 groups (was 53), weighted towards the
  under-fives.** The groups are still shared across the anchor intervals by log
  width, now weighted 1.25 below 5 years and 0.5 from 21: the grid's error comes
  from young children at high transmission, and groups above 21 years move no
  child's incidence measurably. 118 is the smallest grid on which every
  falciparum claim in `fleetcheck` passes; the clinical age profile's 3-5 year
  band at EIR 120 decides it. Severe incidence at EIR 120 is 2% below the IBM
  median, where it was 8%, and a 30-year falciparum run costs about 1.7 s, where
  it cost 0.65 s. `default_age_lower(n_group = 53)` has the new weights, so it is
  not the old grid.

## The daily model

* **`inst/odin/malaria_daily.R` replaces `inst/odin/malaria_ode.R`.** Each day
  reads the state at its start and lands every change at its end, as the IBM
  queues its updates: immunity decays, the day's bites are drawn and
  deduplicated, infection and each state's progression are one competing draw
  per person (so an infection pre-empts that day's progression, as in
  `CompetingHazard$resolve`), the mosquito model is stepped across the day, and
  the dead are replaced by newborns. Only the day's bitten carry an infection
  hazard in the draw, as in the IBM, rather than everyone carrying the
  bite-averaged one.

* **The EIR, human-infectivity and incubation lags are the IBM's own delay
  lines**, not Erlang chains: the EIR `de` days ago and the infectivity
  `delay_gam` days ago, interpolated between two days as `LaggedValue` reads
  `delay_gam`, and the incubation flux that entered `ceiling(dem) - 1` days ago
  -- the IBM's queue holds `ceiling(dem)` days, so the default incubation is 9
  days, not 10.

* **A boost replaces the day's immunity decay.** The IBM queues each
  immunity's decay first and a boost (start-of-day value + 1) later, and the
  later update overwrites the earlier, so a boosted person does not decay that
  day. At high transmission that raises equilibrium immunity by up to a sixth,
  and the continuous-time `dI/dt = boost - I/d` misses it.

* **The larval carrying capacity is held at its start-of-day value** across the
  day, as the IBM's aquatic model holds it (it truncates the solver's time to a
  whole day).

* **The mosquito model is stepped by an exponential midpoint rule**, 32
  sub-steps a day (`ode_tuning(n_sub =)`), which is stable however stiff the
  larval equations get. There is no ODE solver and no tolerance to set.

* **Settings the IBM reads on the day take effect on their own timestep.**
  Treatment coverage, the drug mix, resistance and death rates change on the
  timestep scheduled, one day earlier than before; deployments (nets,
  spraying, vaccination, chemoprevention) act from the day after theirs, and
  their efficacy decays on the IBM's clock, one day since deployment on that
  first day.

* **Prophylaxis chains have daily-clock stage counts.** A stage left with a
  fixed probability a day lasts a geometric number of days, so matching the
  Weibull's mean `M` and variance `V` takes `M^2 / (V + M)` stages: 10 for
  SP-AQ as chemoprevention (was 14), 9 for SP-AQ after treatment (was 16), and
  the chain means are counted in whole protected days.

## Mechanisms brought into line with the IBM

* **Refractory boosting of clinical, detection and severe immunity.** The IBM
  boosts ICA, ID and IVA on an infection only once the refractory window since
  the last boost has passed. `fleet` applied the renewal rate `p / (p u + 1)` in
  every state, as though someone in S or U were reinfected as often as someone
  in A, and so under-boosted. It now carries where the recently boosted have got
  to -- D, A or U, from each day's rates -- and withholds boosts from that
  refractory stock only. Clinical incidence falls by 0.3% at EIR 1 to 3% at EIR
  120, and severe by up to 4%; prevalence moves by under 0.1%.

* **A chemoprevention round sends the detectable it treats through Tr.**
  `update_mass_drug_admin()` moves the clinical and the LM-detectable
  asymptomatic a round treats into Tr, where they stay detectable and infectious
  (at their infectivity x `drug_rel_c`) for `dt`, or `dt_slow` under
  resistance, before protection. `fleet` sent everyone straight to protection,
  and so cut prevalence in the dosed band overnight: at EIR 120 it read 0.12 the
  week after an SMC round where the IBM reads 0.42. The treated then join a
  chain carrying the protection left after Tr, sized as the post-treatment chain
  is.

* **Mass PEV and TBV protect the cohort they vaccinated, as it ages.**
  Protection was written into the targeted age groups for every later day, so
  the vaccinated lost it as they aged out of the band and children born after
  the campaign inherited it (an all-age campaign still cut infant infection by a
  fifth two years on). The band now moves up with the cohort, and repeated
  rounds combine as the most recent dose, as the IBM's overwrite of the dose
  time does. The EPI group straddling the schedule's completion age is
  vaccinated in proportion to the part of it that is old enough.

* **A treated run is seeded at raw treatment coverage**, as the IBM hands it to
  `malariaEquilibrium` and `malariaEquilibriumVivax`, for both parasites. At
  coverage times efficacy, a treated run's mosquitoes had been sized 2-5% below
  the IBM's, and its transmission ran that much lower.

* **Maternal immunity comes from mothers aged 20 to 21**, as in the IBM
  (`trunc(age / 365) == 20`), not from the 20-22.5 year group, which gave
  newborns about 7% too much. `default_age_lower()` pins an edge at 21 years,
  and a newborn inherits the population-weighted mean immunity of the groups
  making up [20, 21), as the IBM's uniform draw of a mother gives.

* **A treated run starts where the IBM starts it.** `set_equilibrium()` sizes
  the mosquitoes under the treatment coverage in force at timestep 1, but the IBM
  draws its initial humans under the coverage at timestep 0. With the usual
  `set_clinical_treatment(timesteps = 1)` the humans therefore start untreated
  among mosquitoes sized for a treated population. `fleet` seeded both treated,
  so its treated runs started apart from the IBM's and met it only after the
  burn-in; the first month of such a run is now within a few per cent of the
  IBM's on every state count.

* **Bites and infectivity are shared by zeta x psi over the live population**,
  as `human_pi()` shares them, where `fleet` assumed zeta averages 1. The
  quadrature's mean is 0.99972 at five nodes (0.975 at three).

## Output

* **The output table is `malariasimulation`'s.** A parameter list gives the
  columns the IBM would give it, under the same names and with the same
  meanings, and `tests/testthat/test-output-parity.R` checks this against the IBM
  across rendering configurations:
  - each column family is rendered over its own rendering list only, with no
    2-10 and all-age bands added to a family whose list is empty, so a bare
    `get_parameters()` gives the 2-10 prevalence columns and no incidence
    columns, as the IBM does;
  - `p_detect_lm_*` is the expected number detected, as in the IBM, not a
    prevalence: prevalence is `n_detect_lm_*` / `n_age_*`, as postie computes
    it. The incidence families gain their `p_inc_*`, `p_inc_clinical_*` and
    `p_inc_severe_*` partners (`p_inc_*` is 0, as `malariasimulation` 3.0.0
    renders it), and `n_infections` is added;
  - a band `[lower, upper]` holds the whole-day ages `lower` to `upper`, both
    ends included, as the IBM's does, where it was `[lower, upper)`: bands that
    share an edge share its day, a single day `[a, a]` is one day's cohort, and
    a band's tag is formatted as the IBM formats it (`182.5_1825.5`,
    `1e.05_2e.05`);
  - `S_count` counts the drug-protected, who are uninfected, as the IBM counts
    them; `Ph_count` says how many of them there are;
  - `ft` is rendered only when clinical treatment is deployed, as in the IBM;
  - `EIR` is the EIR biting humans today, the IBM's lagged value, where it was
    the value `de` days before humans felt it.

* **A run returns one row per day, `1..timesteps`**, as `malariasimulation`
  does, where it returned `0..timesteps`. Row `t` holds the state at the start
  of day `t` and the incidence during it, so row 1 is the seed and the tables
  line up with an IBM run row for row. Code that dropped the first row to align
  with the IBM should stop.

## API

* **`ode_tuning()` holds `age_lower`, `n_ph`, `n_phc`, `n_sub` and
  `odin_file`.** The ODE solver's `atol`, `rtol` and `step_size_max` and the
  Erlang lag counts `n_eir`, `n_foim` and `n_eip` are accepted with a warning
  and ignored, so scripts written for them keep running. An unnamed or
  non-finite field is an error, a stored tuning object is validated again, and
  a list of tuning fields passed as the third argument, `correlations`, is an
  error even when it holds a retired field.

* **The positivity guards count the day's ageing and deaths.** An age group
  loses its ageing and death fraction out of the same stock as its infection,
  progression or prophylaxis exit, so the two together must fit: an age grid
  whose narrowest group cannot hold both is an error, and the prophylaxis
  chains' default stage counts are held to what fits (8 for AL chemoprevention
  on a 209-group grid, where 9 went negative).

* **A refractory window of zero days is no window**, not a negative one.

## Agreement with the IBM

Re-measured in `fleetcheck` with the IBM rows unchanged, on the default 118-group
age grid. Falciparum tier 2 is inside the IBM replicate band throughout the
transmission grid: prevalence, under-5 and all-age clinical and all-age severe
incidence at 6 of 6 EIRs, within 2.2% of the IBM median, and the clinical and
severe age profiles in 25 of 25 and 16 of 16 bands, the clinical profile's
3-5 year band at EIR 120 by 0.04 replicate SD. Intervention impacts are outside
the band in 2 of 72 cells, both bed nets, where the best of the IBM's own
replicates is outside in 9.7%. `fleet` runs 44 times the IBM's speed per
simulated year. On 53 groups, at under half the cost, the grid's discretisation
error puts clinical incidence 3 to 4% and severe incidence up to 7% low at EIR
50 and 120, outside the band at 2 of 6 EIRs. Tier 1: an undisturbed run moves at
most 0.35% off the seed. Tier 3: across the 63-country site files `fleet` tracks
the IBM inside the claim, r 0.983 and 0.958 on clinical and severe incidence
with slopes 0.983 and 0.927, and runs 7% above it, the excess concentrated below
EIR 1 and not explained.

For *P. vivax*, tier 2 is inside the IBM replicate band throughout the
transmission grid, EIR 0.3 to 30: prevalence, under-5 and all-age clinical
incidence, relapses and hypnozoite carriage at 5 of 5 EIRs. The clinical age
profile is inside in 22 of 23 bands: 3-5 years at EIR 10 sits 0.02 replicate SD
above the band, the mean field's school-age offset plus 0.3% from the default
grid, which the falciparum claims set. Radical cure by either drug, treatment
scale-up and bed nets above EIR 1 move every outcome as the IBM does; indoor
residual spraying and bed nets at EIR 1 do not, largely the IBM's relapse
defect. A vivax run costs about what an IBM run does per simulated year, 1.07
times its speed. Tier 1: a run settles, 3 to 15% above the seed on prevalence as
the IBM's does.

## Testing

* `tests/testthat/test-output-parity.R` runs `malariasimulation` beside `fleet`
  and checks the output table column for column across seven rendering
  configurations, `fleet`'s own identities (the `p_*` partners, the whole-day
  band convention, the population counts), and the values against a 50,000-
  person IBM run.

* `tests/testthat/test-daily-mechanics.R` pins each daily-clock mechanism on its
  own: the day settings and deployments take effect, the three delay lines, the
  skipped decay, the competing draw, the seed and the two treatment clocks, the
  chemoprevention treated phase, the refractory closure against an independent
  implementation, the vaccinated cohort ageing, and the positivity guards.

* `tests/testthat/test-vivax.R` pins the vivax block's set-up: the parameter
  translation, the seed on `fleet`'s own age grid, the heterogeneity nodes the
  IBM shares, conservation, radical cure, the column set, and the refusals.
  `tests/testthat/test-vivax-mechanics.R` pins its daily mechanics one at a
  time, as `test-daily-mechanics.R` does for falciparum: the infection and
  relapse draw, PEV on relapses, the hypnozoite ladder, immunity decay, the
  refractory windows and the boost that skips a day's decay, the liver-stage
  clock and the delay lines. `test-output-parity.R` checks the vivax table
  against `malariasimulation`'s.

* `tests/testthat/reference-values.csv` and `reference-interventions.csv` are
  regenerated on the 209-group grid, and the intervention reference is sampled
  on the days each effect first reaches its quantity. `reference-vivax.csv` and
  `reference-vivax-interventions.csv` pin the vivax block the same way: its
  output levels untreated and under radical cure by primaquine and by
  tafenoquine, and the same three intervention probes on vivax's own delay
  lines, with a check that none acts before its deployment day.

# fleet 0.0.0.9002

Two changes to how age is discretised, both of which move model output. The
ageing rate is now exponentially fitted, and the default age grid is log-spaced
rather than made of fixed monthly / quarterly / yearly / 5-yearly sections.
Together they halve the discretisation error at essentially the same cost.

## The default age grid

* **`default_age_lower()` is log-spaced between pinned anchor ages.** Group
  edges are pinned at the conventional reporting boundaries — 0, 1, 2, 3, 5, 7,
  10, 15, 20, 30, 40, 60 years — so the usual output bands still fall on group
  boundaries exactly and every band weight is 0 or 1. Between two anchors the
  groups are of equal width, and the budget of groups is shared across the
  anchors in proportion to **log** width.

  Log width, because what the grid has to resolve is the rise of immunity with
  age, and that is much closer to a function of log age than of age. The old
  grid spent 12 of its 52 groups in the first year of life and 4 on the whole of
  40–60; this one spends 7 and 5.

  Measured against `fleet` run to grid convergence, at the same number of
  groups: the largest departure across age bands falls from **10.7% to 4.5%**
  and the rms from **5.4% to 3.2%**, for a few per cent of extra runtime. An
  equal-width grid at the same cost is 13.0% rms — 2.4× the old grid and 4× this
  one — so the grading matters more than the count.
  `validations/age-grid/run.R` in `fleetcheck` reproduces all of it.

  This does **not** close the gap to `malariasimulation`'s age profile of
  clinical incidence, and refining the grid further will not either. The two
  errors behave differently, which is the whole point of measuring them apart:
  `fleet`'s departure from its own converged profile falls 4.5% → 2.3% → 1.1% →
  0.6% as the group count doubles (a ratio of 2.0, clean first order), while its
  departure from the IBM falls 7.7% → 5.3% → 4.6% → 4.3% and levels off.
  Discretisation is removable; what is left is the mean-field approximation —
  one immunity value per stratum against a spread of individual infection
  histories at the same age. That spread is between people, not across an age
  group, so no grid touches it.

* **`default_age_lower()` gained `n_group`** (default 53). Refining the whole
  grid while keeping its shape is now `default_age_lower(n_group = 105)` rather
  than a bespoke script, which is what a grid-convergence check wants.

* Age boundaries that do not land on a group edge — an SMC campaign targeting 3
  to 59 months, say — are **apportioned by exact fractional overlap**, both for
  output bands and for intervention targeting. Nothing is snapped to the grid,
  so a boundary the old grid happened to have an edge on has not become less
  accurate, only differently discretised.

## Numerics

* **The ageing rate is exponentially fitted rather than `1/width`.** A linear
  chain empties an age band at rate `r`, so the stationary ratio between
  consecutive bands is `r / (r + mu)`. With the obvious `r = 1/h` that is
  `1 / (1 + mu*h)`, where the continuous solution of the McKendrick equation is
  `exp(-mu*h)` — and since `1/(1 + x) > exp(-x)` for every `x > 0`, the obvious
  rate **always** decays too slowly and always leaves too many people alive at
  old ages. It is a first-order donor-cell discretisation of an advection.

  Setting `r = mu / expm1(mu*h)` makes that ratio exact, so the stationary age
  structure is the analytic band-integrated survival curve rather than an
  approximation to it. This is exponential fitting, the standard cure for a
  first-order advection scheme.

  It tends to `1/h` as `mu*h -> 0`, so it is nearly a no-op where bands are
  narrow — 0.3% in the narrowest bands — and does its work in the wide ones,
  11.4% in the 5-year bands above 60. Against 20 `malariasimulation` replicates the
  worst age-band error fell from 7.4% to 1.3%, and the population age structure
  went from 9 of 11 bands inside the IBM's replicate band to 11 of 11.

  `r` now depends on the death rate, which reads oddly: ageing should not depend
  on dying. It is not a biological rate. It is the coefficient that makes the
  discrete scheme reproduce the continuous solution, and that solution involves
  `mu`.

  It is fitted at the baseline (`t = 0`) mortality. Under `set_demography()`
  with time-varying rates `mu` moves and this coefficient does not follow it, so
  the fit is exact at the seed and an improvement, not an identity, afterwards.
  Making `r_age` time-varying in the odin model would be the next step if a
  transient ever needs it.

  The seed is unaffected. `malariaEquilibrium` assumes `r = 1/width`, so the
  constant-hazard reference the equilibrium arrives on keeps that convention and
  is computed separately.

## Testing

* `tests/testthat/reference-values.csv` and `reference-interventions.csv` are
  regenerated. Every pinned age-resolved output moved.

* The time-varying-mortality test in `test-site-features.R` asserts a smaller
  transient, 3% rather than 5%: a more accurate age structure shrank the
  movement it is testing for from 6% to 4%. What the test asserts is unchanged.

* `test-package.R` no longer hardcodes `r_age = 1/5y` for the 5-yearly section.
  It asserts the exponential-fitting identity against the rate actually in use.

# fleet 0.0.0.9001

The package was renamed from `blink` to `fleet`. Otherwise this release is the
defects a full package review found, fixed. Three of them move model output:
the age-band rendering weights, five intervention series that took effect before
their scheduled date, and TBV / mass-PEV age targeting. The rest turn inputs that
used to be accepted silently into errors or warnings, switch the comparison CI
on, and give the test suite something to fail on. Output on bands that line up
with the model age grid, in a run with no interventions, is unchanged.

## Package rename

* **`blink` is now `fleet`, and the repository moved to
  <https://github.com/pwinskill/fleet>.** The old name collided with an existing
  CRAN package (Steorts, *Record Linkage for Empirically Motivated Priors*), and
  the collision was not cosmetic: pkgdown builds its "View on CRAN" link from the
  package name alone, so the site stood to send readers to someone else's
  release. Install with `remotes::install_github("pwinskill/fleet")` and
  replace `library(blink)` with `library(fleet)`.

  **No function signatures changed.** The three exported functions
  (`run_simulation_ode()`, `ode_tuning()`, `default_age_lower()`) keep their
  arguments, and nothing moved numerically: 372 tests passing, both reference
  CSVs byte-identical, and the drift check reproducing every committed row to
  1e-6. Three things do follow the name: the environment variables (`BLINK_LIB`,
  `BLINK_ALLOW_SKIP`, `BLINK_REGENERATE_REFERENCE`, `BLINK_VALIDATE`,
  `BLINK_DIAGRAM_OUT` and `CMP_BLINK_ONLY` are now `FLEET_*` /
  `CMP_FLEET_ONLY`), the `blink.pev_gq` option, now `fleet.pev_gq`, and the
  `model` column of `comparison/data/rep_*.csv`. `blink2_validate` is
  deliberately untouched: it names a separate repository.

## Bug fixes

* **Two coefficients in the model were frozen at the seed's age structure and
  are now taken from the live population.** Immunity ageing used
  `re[i] = r_age[i] + mu_age[i]` as the per-capita inflow into an age cell, and
  the FOIM normalisation used the seed's mean `psi`. Both are exact only while
  the age structure is stationary. Under `set_demography()` with time-varying
  death rates -- which every malariaverse site file supplies, yearly from 2000
  -- the structure drifts for decades, and `re[i]` was built from the new death
  rate while the population still had the old shape. The inflow is now
  `r_age[i-1] * N[i-1, j] / N[i, j]` (newborns at `i = 1`), which the
  intensive-variable derivation gives directly, and the FOIM denominator is the
  current `sum(psi * N) / sum(N)`, which is what malariasimulation's `human_pi`
  computes on the live population. Both equal the old forms to machine
  precision at the seed, so **nothing changes while death rates are constant**:
  every committed comparison row reproduces to 1e-6, and a synthetic run is
  identical to 6e-15 until the day its death rates change.

  Where they do change it is material, and measured against the IBM rather than
  argued from direction. Afghanistan's site-file mortality series, 2000-2029,
  no interventions, against 20 `malariasimulation` replicates: over the last
  ten years the bias against the IBM median goes from +3.6% to -0.3% on EIR,
  +0.9% to -0.3% on PfPR(2-10), +2.8% to +1.4% on all-age clinical incidence,
  +7.6% to +5.3% on under-5 clinical incidence, and +1.7% to +0.9% on under-5
  severe. Years inside the IBM's 10-90% replicate band go from 21 to 27 of 30
  on EIR, 28 to 30 on PfPR and 17 to 24 on under-5 clinical incidence; nothing
  moves out of the band. All-age severe is the one outcome that goes slightly
  the other way, -1.3% to -1.9%, inside a replicate band 30% wide and 29 of 30
  years inside it either way. A synthetic mortality step gives the same picture
  with every outcome closer. The human age structure is untouched, as it should
  be: the population total is conserved to 1e-15 and its age split moves 2e-6,
  which is solver tolerance over 30 years.

  `mean_psi` is no longer a model parameter.

* **Output age bands are rendered by exact overlap, not by age-group midpoint.**
  A model age group used to contribute its whole population to a rendering band
  if its midpoint fell inside the band, and nothing at all otherwise. A band
  narrower than the group it lands in captures no midpoint, so the whole column
  family came back as hard zeros: 52 of the 65 single-year bands from 15 to 80
  years on the default grid, where `p_detect_lm_7300_7665 = 0` reads as "no
  malaria in 20-21 year olds" rather than "no data". And a partial overlap
  counted the whole group or none of it, so the standard 6-59-month band came
  back 2.2% light on a default parameter list, flowing into every PfPR
  denominator and every postie person-day. A group now contributes the
  fraction of its own width lying inside `[min_age, max_age)`, so a group
  straddling a boundary is split between the two bands and the weights sum to
  exactly 1 per group over any set of bands that partition a span of the age
  axis. `fleet` also warns now when a rendering band overlaps no model age group
  at all. **This moves published numbers for any band that does not align with
  the model age grid**, in every scenario; bands that do align are unaffected,
  and regenerating the unit-test reference outputs produced a byte-identical
  file.

* **Bed nets, IRS, carrying capacity, PEV and TBV no longer take effect before
  their scheduled date.** All five series are interpolated linearly, which is
  right for the within-round decay each carries, but a value that changes on a
  scheduled day ramps in from whatever knot precedes it, and the background grids
  are coarse: 10 days for vector control, 30 for TBV and PEV, 5 for carrying
  capacity without seasonality. So an intervention arrived early, in some cases
  almost entirely. Nets scheduled for day 100 had reached 66% of their effect on
  EIR by day 99, before a single net existed; a TBV round on day 400 had
  delivered 88% of its transmission blocking by day 399; a mass PEV campaign
  whose efficacy onset was day 490 had delivered 34% of its FOI reduction by day
  485; a carrying-capacity halving at day 100 started moving EIR on day 99. Each
  grid now carries a knot at `onset - 1` as well as at every scheduled day,
  confining the ramp to the single step onto the scheduled day: the resolution
  the chemoprevention pulses already had. Every scenario using one of these five
  interventions moves.

* **TBV and mass-PEV age targeting weights by the fraction of each group
  covered.** Both used to select the model age groups whose midpoint fell in the
  target, so a target that straddled no midpoint selected nothing and the whole
  campaign was a silent no-op. Above age 14 the default grid is 5-yearly and only
  the years 17, 22, 27 and so on are midpoints, so `set_tbv(ages = 18:20)`
  vaccinated nobody, as did a mass PEV band of `min_ages = 20 * 365`,
  `max_ages = 21 * 365`, which falls between the midpoints at 17 and 22 years.
  Each group now carries the fraction of its own width inside the target; a
  target aligned with group edges gives weight 1 and reproduces the old result
  exactly, and one overlapping no group at all now warns.

  Fractional weights then exposed a second defect in mass PEV. The several
  `min_ages`/`max_ages` bands of a *single* campaign were folded in one at a
  time with the independent-protection rule `1 - (1-a)(1-b)`, which is the rule
  for separate campaigns, not for two bands of the same one; the midpoint rule
  could not expose it, because a group belonged to at most one band. Splitting
  one band into two at 1195 days, which should change nothing, moved a
  straddling group's covered fraction from 0.7295 to 0.6834, a 6.3% shortfall.
  The bands of one campaign now add into a single covered fraction (capped at 1)
  before coverage is applied; independent combination is kept across campaigns.

* **Inputs that could not mean what they said are rejected rather than run.**
  Each of these used to be accepted, or to fail somewhere downstream with a
  message that never named the parameter at fault:

  * `parameters$n_heterogeneity_groups` of 1 or 2. Fewer than three
    Gauss-Hermite nodes do not integrate the log-normal biting distribution,
    they just move every human off it: at `n = 1` the single node sits at
    `zeta = exp(-s2/2) = 0.434`, so every human in the model experienced 43% of
    the EIR the output column reported. For a genuine no-heterogeneity run set
    `parameters$enable_heterogeneity = FALSE`, which collapses to one stratum at
    `zeta = 1` exactly.
  * `parameters$bite_dedup` outside 0/1 (`TRUE`/`FALSE` still accepted). It is a
    switch between two hazard forms, not a dial, and 0.5 was silently a third
    thing.
  * `parameters$acquired_immunity_offset` outside `[0, 1]`.
  * `default_age_lower(max_age =)` below 20 years, which made the grid's
    5-yearly section `seq(15, max_age - 5, by = 5)` run backwards and died with
    base R's "wrong sign in 'by' argument".
  * `timesteps` that is not a single finite whole number `>= 1`. This also
    catches the swapped-argument call `run_simulation_ode(parameters,
    timesteps)`, which used to die on `parameters$init_EIR` with "$ operator is
    invalid for atomic vectors".

* **Editing the parameter list after `set_equilibrium()` now warns.**
  `set_equilibrium()` stores its back-translation as `parameters$eq_params`, and
  `build_inputs()` merges that over the live translation, so the stored copy wins
  for every shared biological constant and an edit made afterwards
  (`p$du <- 10`) is a silent no-op. That contract stands — `set_equilibrium()`
  is meant to be the last call on the parameter list — but a violation of it
  must not be silent. `fleet` now compares the two and warns, naming the
  `malariasimulation` fields whose live translation no longer agrees with the
  frozen copy, so the message points at the field you typed rather than at its
  equilibrium alias. A value that took a different arithmetic route to the same
  number cannot trigger it, and a deliberate custom `eq_params` is called out as
  expected.

## Numerics and performance

* **One `exp()` per cell where there were two.** `p_inf` now reuses `q_b`, which
  is the same `1 - exp(-EPS)`; odin2 does no common-subexpression elimination,
  so the two spellings cost two evaluations per cell per derivative call. The
  right-hand side is bound by its ~2,600 transcendentals per evaluation and this
  was one in ten of them by count -- but `exp` is cheap next to the five
  non-integer `pow` calls per cell, so the measured saving is 3-5% of the
  derivative cost (22.5 to 21.2 microseconds per evaluation on an aseasonal
  run) and 2-6% of wall time. Output is bit-identical (the same operation,
  evaluated once): the solver takes exactly the same steps, and every committed
  comparison row reproduces to 1e-6.

* **Where the solver's time goes, measured rather than asserted.** On an
  equilibrium run every accepted step is 0.62 days, and that does not move with
  `atol`, `rtol`, `step_size_max`, the output grid or the delay-chain lengths
  below 4 per day. It is an explicit-stepper stability limit: the stiffest
  eigenvalue is the late-larval density-dependent mortality
  `ml * gamma * (E + 2L) / K`, 5.3 per day at the seed, and 3.3 / 5.31 = 0.621
  days. That is intrinsic to malariasimulation's larval model, so the ~1.6
  steps per day it forces is the floor for this solver class, and the daily
  output grid adds ~2,700 rejected trial steps over 30 years on top. On seasonal
  runs the tolerance does bind: `rtol = 1e-6` is 1.44x faster with a maximum
  deviation of 4e-7 in daily clinical incidence. The `ode_tuning()`
  documentation now says this; it previously said the cap and the daily grid set
  the aseasonal step count, which was a reasonable guess that measurement did
  not support.

* **The default age grid is not converged for severe disease or for adult
  clinical incidence.** Against 20 IBM replicates at EIR 20, all-age severe
  incidence goes from -2.8% below the IBM median on the default grid to +0.2%
  on a grid with quarterly bands to 15 y and yearly bands to 40 y (110 groups,
  2.1x the run time), and adult clinical incidence in the 30-60 y bands from
  +14-16% to +6%; prevalence and under-5 clinical incidence do not move.
  Both are Jensen error from averaging steep age dependence over wide bands.
  The default is unchanged, because every published comparison is stated on
  it, but the `age_lower` documentation now says which results need a finer
  grid, and the two claims in the fleetcheck register that rest on severe
  incidence by age (`severe-allage-bias`, `age-structure`) carry a
  discretisation component of the same size as the discrepancy they report.

## Continuous integration

* **All four workflows are live**, so something now checks the model without
  being asked. The triggers that `0.0.0.9000` left staged and commented are
  switched on, as are the `paths-ignore` and `concurrency` blocks staged
  alongside them: `comparison.yaml` on every push and pull request that can touch
  the model plus a Monday 06:00 UTC schedule, and `figures.yaml` on anything
  touching the comparison data, the renderers, the theme or the shared constants.
  The Monday run is the one with no substitute: it is the only thing that catches
  `malariasimulation` moving underneath the comparison.
  `.github/workflows/CI.md` is now a description of what runs, not a go-live
  checklist, and carries the detail.

* **`R-CMD-check.yaml` gains a job that fails when the compiled model is stale.**
  `inst/odin/malaria_ode.R` is the model's source of truth, but what gets
  compiled is the committed, generated `src/malaria_ode.cpp` (with `R/dust.R`,
  `R/cpp11.R`, `src/cpp11.cpp` and `inst/dust/`), which no ordinary build
  regenerates. An edit to the odin model without a following
  `odin2::odin_package(".")` therefore left every runner compiling the *old*
  model, the tests passing, and `check_drift.R` reporting "nothing moved" —
  which reads as *numerically inert* when it means *never compiled*. The new job
  regenerates from `inst/odin` and fails if the working tree comes out dirty.

## Testing and tooling

* **The suite pins absolute output, not just relationships between outputs.**
  It could not previously detect a wrong answer: the only severe-incidence
  assertion was a ratio, invariant to any scale factor, so doubling `sev_inc_a`
  passed every test. `tests/testthat/reference-values.csv` and
  `reference-interventions.csv` now pin recorded output to 1e-6 across the EIR
  grid, and for each intervention at the day before deployment, at onset and
  after — the shape that catches an effect arriving early.

* **A missing `malariasimulation` fails the suite instead of quietly emptying
  it.** Nearly every `test_that()` block opens with
  `skip_if_not_installed("malariasimulation")`, so a machine where the package
  was missing or unloadable ran the handful of dependency-free blocks, skipped
  the rest, and reported success having tested almost nothing.
  `tests/testthat/setup.R` now errors, with `FLEET_ALLOW_SKIP` as the deliberate
  escape hatch.

* **`NAMESPACE` is roxygen-generated.** It was hand-written and held the only
  copy of `useDynLib(fleet, .registration = TRUE)`, so the first
  `devtools::document()` would have dropped it and left the package unable to
  load.

## Documentation

* **The docs are split into articles and the README is a front page again.** It
  was 250 lines, roughly 170 of which restated the articles at summary level,
  which is how its claims went stale. Each topic now has one home:
  `vignette("fleet")`, the new `vignette("using")`, `vignette("model")` (the
  specification alone), the new `vignette("parameters")` and
  `vignette("comparison")`. `.github/CI.md` moved to `.github/workflows/CI.md`,
  where pkgdown no longer renders it as the site's front page.

# fleet 0.0.0.9000

First development release: a deterministic mean-field (ODE) twin of the
`malariasimulation` individual-based model of *P. falciparum* malaria, built on
odin2/dust2.

## API changes

* **`run_simulation_ode()`'s signature now matches
  `malariasimulation::run_simulation()`.** It was
  `(timesteps, parameters, correlations, init_EIR, age_lower, n_eir, n_foim,
  n_eip, n_ph, n_phc, atol, rtol, step_size_max, odin_file)`: fourteen
  arguments, eleven of which the IBM does not have. It is now
  `(timesteps, parameters, correlations, tuning)`. Two changes, both breaking.

  **`init_EIR` is gone as an argument.** The target EIR is a model input, so it
  belongs on the parameter list, where `malariasimulation::set_equilibrium()`
  already puts it as `parameters$init_EIR` and where the IBM reads it from:

  ```r
  p <- malariasimulation::set_equilibrium(p, init_EIR = 20)
  out <- run_simulation_ode(3650, p)
  ```

  Setting `parameters$init_EIR` by hand still works, but skips the `eq_params`
  that `set_equilibrium()` stores and `fleet` honours. This is numerically inert
  for anyone already calling `set_equilibrium()`: the whole comparison harness
  reproduces bit-identically through the new API.

  **Everything numerical moved into `tuning`**, a new `ode_tuning()` object (also
  accepted as a plain named list of just the fields you want to change). That is
  `age_lower`, the five Erlang stage counts `n_eir`/`n_foim`/`n_eip`/`n_ph`/`n_phc`,
  the solver controls `atol`/`rtol`/`step_size_max`, and `odin_file`:

  ```r
  run_simulation_ode(3650, p, tuning = list(rtol = 1e-6, step_size_max = 10))
  ```

  Nothing epidemiological sits outside the parameter list any more: every field
  of `tuning` is a numerical-approximation knob. `ode_tuning()` validates its
  fields, and an unrecognised name in a `tuning` list is an error rather than a
  silently ignored field. Passing any of the removed arguments gives a targeted
  error naming the call that replaces it, rather than R's bare "unused argument".

* **`get_epi_outputs()` is removed** (#3). It was a thin wrapper around
  `postie::get_rates()` and `postie::get_prevalence()`, plus argument-splitting
  machinery (`rates_args`, `prevalence_args`, `...`) that was more to learn than
  the two calls it replaced. malariasimulation has no such wrapper and `fleet`'s
  output table is malariasimulation-shaped, so IBM post-processing pipelines
  already work on a `fleet` run unchanged. Replace `epi <- get_epi_outputs(out)`
  with `postie::get_prevalence(out, diagnostic = "lm")` and
  `postie::get_rates(out)`. `postie` remains in Suggests for the vignette and
  tests, but nothing `fleet` exports depends on it.

## Behaviour changes

* **Custom demography: mosquito sizing now replicates `set_equilibrium()`.**
  malariasimulation sizes the adult-mosquito population from the human
  equilibrium under its *default* exponential age structure (`set_equilibrium()`
  → `equilibrium_total_M()`; custom mortality never enters), so under
  `set_demography()` the IBM drifts to whatever transmission that density
  supports. `fleet` used to re-solve the equilibrium under the custom age
  structure and hold `init_EIR` exactly, so the same parameter list realised
  different transmission in the two models (EIR 20 in fleet against 13.8 in the
  IBM in the comparison scenario, PfPR(2–10) 0.55 against 0.49). It now
  reproduces the IBM's `total_M` exactly (`ibm_total_M()`), root-finds the EIR at
  which its own equilibrium under the custom demography has that density, and
  seeds there — still a fixed point, so there is no burn-in (the scenario now
  seeds at EIR 13.86, PfPR 0.489). Under the default demography nothing changes:
  fleet keeps its own sizing, which lands ~1.5% above the IBM's `total_M`, so the
  same exponential mortality expressed through `set_demography()` seeds at
  EIR 19.7 rather than 20, within the IBM's replicate noise.
  `parameters$hold_init_EIR = TRUE` restores the previous behaviour. Every
  malariaverse site file uses `set_demography()`, so this affects all country
  work: across the 63-country validation set it moves fleet's mean
  clinical-incidence excess over the IBM from +10.1% to +8.7% and the regression
  slope from 1.02 to 1.00 (r 0.980 → 0.982); severe from +9.4% to +8.8%.

* **Prophylaxis is an Erlang chain, not one exponential compartment.** The IBM
  applies the Weibull survival `W(t - t_drug)` to each treated person's infection
  probability, so a treated cohort's mean protection at lag `t` is exactly `W(t)`.
  `fleet` represented both prophylaxis compartments (`Ph` post-treatment, `Ph_c`
  chemoprevention) as a single exponential at the Weibull mean, which leaks
  protection early — for SP-AQ (shape 4.3, scale 38.1) only 42% are still
  protected at day 30 against the Weibull's 70% — and under-estimated seasonal
  SMC (under-5 clinical reduction 40% vs the IBM's 52% in the comparison article;
  now 54% against 52%, within the IBM's replicate range). Both are now Erlang
  chains moment-matched to the Weibull. The chemoprevention chain has `k = 1/CV²`
  stages (14 for SP-AQ, 15 for DHA-PQP; SP-AQ chain at day 30: 0.67). The
  post-treatment chain follows the exponential treated stage `Tr` while the IBM's
  clock runs from the dose, so its mean is the integrated protection left after
  `Tr` (`mean_W - ∫exp(-r_T t) W(t) dt`) and its length matches the variance of
  the whole `Tr + Ph` sojourn: 16 for SP-AQ, 20 for DHA-PQP, and 1 for AL, whose
  10-day protection is already less variable than `Tr` itself (measured
  output-identical to a 20-stage chain at half the run time). Counts are capped
  at 20; `ode_tuning(n_ph =, n_phc =)` overrides them and 1 recovers the old
  behaviour. The equilibrium seed distributes the prophylaxis mass across the
  stages exactly (`solve_disease_block()` generalised), so runs still hold flat,
  and runs now record only the output variables rather than the full state,
  keeping memory flat despite the larger state.

## Bug fixes

* **Severe and all-infection incidence are counted with the deduplicated
  probability, not the raw hazard.** This is the change in this release that
  moves the no-intervention numbers; the review-pass fixes in `0.0.0.9001` move
  intervention scenarios and misaligned rendering bands, and nothing else.
  `clin_inc_a` was already counted with `h_c`; `sev_inc_a` and `inc_a` used
  `FOI`. These outputs are rates that R integrates over a day, and the right rate
  depends on how fast the compartment drains. `S` and `U` drain at the full
  `FOI`, so `int FOI*X exp(-FOI t) dt = X*(1 - exp(-FOI)) = X*p_inf` — the IBM's
  count, exactly. `A` does not: a sub-clinical re-infection of an `A` leaves them
  in `A`, so `A` drains only at `h_c` and stays roughly flat over the day, making
  `int FOI*A dt ~ FOI*A` and over-counting by `FOI/p_inf = -log(1-p_inf)/p_inf`.
  `A` now takes `p_inf` directly. The IBM draws severe from the *deduplicated*
  `infections` bitset, so this is a replication defect, not a modelling choice.

  Against 3 IBM replicates of 10,000 people on the same parameter list, all-age
  `n_inc_*` ran **+2.6% (EIR 20) and +4.2% (EIR 50)** above the IBM median, both
  outside its replicate spread; after the fix, **+0.6% and +0.4%**. The bias was
  concentrated in adults, where `A` is a large share of the at-risk pool (under-5
  `n_inc_*` was only +0.3% / +1.7%), the signature of an `A`-only defect.
  Clinical incidence and prevalence are unchanged.

  It also makes one comparison worse, and that is worth stating plainly. The same
  correction lowers severe incidence slightly, and fleet already ran below the
  IBM there: all-age severe goes from −5.9% to −6.2% at EIR 20 and from −3.4% to
  −4.0% at EIR 50, dropping just outside the IBM's 10–90% band at those two EIRs
  (by 0.2% and 0.3% of the lower edge) where it had been marginally inside. The
  upward bias had been masking part of a larger, separate severe deficit. That
  compensation was accidental, so the corrected numbers are the honest ones and
  the severe gap is now the largest open discrepancy between the two models.
  `vignette("comparison")` and its tables are re-rendered accordingly.

* **`test-prophylaxis-chain.R`'s first test no longer errors when
  `malariasimulation` is absent.** It read `AL_params`/`SP_AQ_params` without a
  `skip_if_not_installed()` guard, the only one of the suite's 57 test blocks to
  do so. `malariasimulation` is a GitHub Remote in Suggests, so on any check
  machine that could not install it the other 56 blocks skipped and this one
  errored, turning a missing optional dependency into an `R CMD check` failure.
  The block is now split rather than simply guarded: its closed-form
  `erlang_stages()` checks keep running unguarded, and only the Weibull
  parameters read out of the IBM's tables sit behind the guard.

* **`build_inputs()` no longer recomputes grid-invariant PEV protection.** Two
  fixes in `pev_series()`, both verified bit-identical on the full
  `[n_age, n_time]` multiplier matrix across the EPI single-timestep path, the
  EPI multi-timestep path and mass campaigns. `.pev_eff_fun()` built its
  interpolator with `stats::approx()`, which re-runs `regularize.values()` on
  *every* scalar call; it now builds one `stats::approxfun()` closure. And the
  EPI loop called `pev_protection()` once per grid point when
  `.booster_cov_vec()`'s admin-time argument is ignored (which the single-row
  `booster_coverage` that `set_pev_epi()` builds by default always is), so that
  call is now hoisted to once per age band. Together: `pev_series()` on a 20-year
  two-booster EPI schedule drops from 2.92 s to 1.09 s (2.7x).

* **Chemoprevention pulses renew existing protection and clear `Tr_slow`.** The
  IBM's `update_mass_drug_admin()` resets `drug_time` for everyone successfully
  treated whatever their state; `fleet`'s pulse left people already in `Ph`/`Ph_c`
  decaying from their earlier dose and skipped the slow-clearance treated
  compartment. Both now move to the first stage of `Ph_c` with the rest of the
  covered fraction — material for monthly SMC rounds, where most of the previous
  round's recipients are still protected.
* **`rT_slow` uses the whole-day exit probability.** The slow-parasite-clearance
  rate under antimalarial resistance was `1/dt_slow`; it is now
  `1 - exp(-1/dt_slow)`, the same conversion `build_inputs()` applies to
  `rA`/`rD`/`rU`/`rT` (the IBM leaves states with per-day probability
  `rate_to_prob(1/d)`).
* **Custom demography: top age group's death rate.** The open-ended oldest model
  age group is represented by its lower bound, which the right-closed
  `set_demography()` bins assigned to the band *below* it whenever that bound
  coincided with a bin edge: 80 y on the default grid, a very common edge. The
  over-80s therefore died at the 60–80 rate (0.05 instead of 0.12 per year in the
  comparison scenario) and the group held ~2.4× too many people (60–85 share 18.2%
  against an analytic 13.1%). The rate is now looked up one day above the group's
  lower bound. Found by the new comparison article.

## Documentation

* New article `vignette("comparison")` — *Comparison with malariasimulation* —
  covering EIR, age, seasonality, demography, 63 country site files and five
  interventions, with the IBM as the median of 10 replicates and a 10–90% band.
  `comparison/` was rewritten around one runner, one renderer and a shared theme,
  with the per-replicate summaries committed so figures redraw without re-running
  the models; the pre-replication-pass figures (`A_`–`F_`, `inc_*`) are gone.
* The comparison surfaced three fleet-side items, all fixed above: the top age
  group's death rate under `set_demography()`, exponential prophylaxis
  under-estimating seasonal SMC, and the `set_equilibrium()` convention under
  custom demography.
* A run-time benchmark (#1), reproduced by `comparison/benchmark.R`, now the
  *Cost model* section of `vignette("using")`.
* New figure `cmp_programme_ts` and a *Programmes over fifteen years* section,
  from five new `ts_*` scenarios running 15 years past deployment at EIR 20.
* The intervention-impact and equilibrium-vs-EIR figures each gain all-age
  clinical and all-age severe panels, laid out 2×2; `summary_tables.R` reports
  both new outcomes.
* **The 63-country site comparison is now an explicit snapshot**,
  `comparison/data/site_snapshot.json`, read by `summary_tables.R` in place of a
  live run against a checkout most people do not have; `CMP_REFRESH_SITES=1`
  re-takes it.
* New workflows `comparison.yaml` (the drift check) and `figures.yaml` (do the
  committed figures and tables still match the committed data), plus
  `dependabot.yml`. All landed inert and were switched on in `0.0.0.9001`; see
  `.github/workflows/CI.md`.
* **New `comparison/check_drift.R`: check the match without re-running the IBM.**
  It separates *did anything move* from *is the match still good*, failing only
  on the second unless `CMP_STRICT=1`, and `comparison/data/ibm_reference.json`
  records what the committed IBM rows were made from. The scenarios, summariser
  and run loop moved into `comparison/scenarios.R`.
* **The model flow diagram is simplified.** Its layering now means only the human
  panel's age × heterogeneity grid, and the stage-count outlines, at-risk rings
  and immunity box are gone. The FIDELITY invariant cites `deriv(A)` and
  `deriv(U)` rather than odin line numbers that had rotted.
* Every line figure draws the IBM's dashed median *over* fleet's solid line
  rather than under it, where agreement had been hiding the IBM entirely. The
  rule is in `comparison/theme.R`; the unused `geom_ibm_envelope()` is removed.
* The country-site figure reads its 1:1 correlation properly (#2): an orange
  dashed reference line rather than a series colour, a near-white low end to the
  log10 hex ramp, and legend breaks at 1 / 10 / 100 / 1,000 / 10,000.
* Three stale comment blocks are removed, each left *above* the correction that
  replaced it: `drug_mix()`'s weights and `rP` in `R/interventions.R`,
  `beta_eff` in `inst/odin/malaria_ode.R`, and `rT_slow` in
  `tests/testthat/test-package.R`.

## Exact-replication pass (malariasimulation v3.0.0)

A systematic audit replaced every place `fleet` approximated a mechanism whose
malariasimulation implementation is known. No tuned constants were introduced.

* **Output bands — bug fix, affects every user.** `render_output()` emitted
  `n_inc_clinical_*`/`n_inc_severe_*`/`n_inc_*` over the *union* of all rendering
  bands, producing overlapping strata that `postie` double-counted, inflating
  all-age clinical by ~1.21x and severe by ~1.10x. Each family is now emitted
  over its own `*_rendering_ages`, with `n_age_*` over the union.
* **Maternal immunity source band.** An exact tie in
  `which.min(abs(age_mid - 20*365))` drew `ICM`/`IVM` from 15-20 year-olds
  instead of the 20-21y mothers the IBM uses; they ran 9-15% low. It now selects
  the band containing 20y.
* **EIP survival.** The Erlang chain carried mortality in every stage
  (`(reip/(reip+mum))^n_eip`, +4% at baseline); it is now a loss-free delay with
  the IBM's `exp(-mum*dem)` applied at the exit.
* **Slow parasite clearance.** `Tr` is split into parallel `Tr`/`Tr_slow`,
  replicating the IBM's Bernoulli assignment of each treated individual to `dt`
  or `dt_slow` rather than one exponential at the blended mean.
* **Clinical event counting.** The clinical hazard is now `-log(1 - phi*p)`,
  which integrates to exactly the IBM's `phi*p*N`; `phi*FOI` over-counted.
* **State sojourns.** Rates now use the IBM's per-day exit probability
  `1-exp(-1/d)`, so the realised mean dwell matches (5.52 d for `dd = 5`, not 5).
* **Immunity boosting.** The boost rate is now `q/(q*u_eff + 1)` with
  `q = 1-exp(-rate)` and `u_eff = ceil(u) - 1`, and `ICA`/`ID`/`IVA` are boosted
  only for individuals eligible to be infected (`S`/`A`/`U`), `IB` for everyone
  bitten.
* **PEV antibodies.** Efficacy is integrated over the antibody distribution the
  IBM samples (4-D Gauss-Hermite over `cs`/`rho`/`ds`/`dl`) rather than evaluated
  at the profile median, which overstated R21 efficacy by up to ~4.7 pp.
  Validated against the IBM's own sampler to <5e-4.
* **PEV gating** follows `pev_epi_coverages` / `pev_epi_timesteps` /
  `mass_pev_timesteps`, not `parameters$pev`, which the IBM never reads.
* **Seeding consequence.** Because `malariaEquilibrium` encodes the simplified
  forms, the analytic seed is now a close approximation rather than an exact
  fixed point: an undisturbed run relaxes by up to ~0.5% over the first years, as
  malariasimulation itself does.

### Validation after this pass

The pass was validated monthly, *P. falciparum* only on both sides, over 63
countries / 1,391 sub-sites / 450,684 sub-site-months against pre-run
malariasimulation output, and moved every headline statistic the right way.
Those numbers, the committed snapshot they were re-taken from
(`comparison/data/site_snapshot.json`, 2026-09-09) and the discrepancies still
open are in `vignette("comparison")` and in the warning at the top of the
README, characterised there rather than tuned away.

## Parameter-ingestion completeness

A systematic audit of every user-facing `malariasimulation` config function
(`get_parameters` + all `set_*`) found gaps where the model was not using the
full parameter flexibility they expose, and closed them: the PEV booster sequence
beyond the first booster, every target age band of a mass PEV campaign,
fractional-overlap age targeting for chemoprevention, antimalarial resistance
applied to chemoprevention drugs and not only to clinical treatment, IRS
accumulating across spray rounds, TBV's exact `ages` set, and
`set_equilibrium`'s custom `eq_params`. All-infection incidence (`n_inc_*`;
`incidence_rendering_*`) was added as a new output. The two remaining mean-field
approximations — the time-varying multi-drug first-line switch and seasonal PEV
boosters — are warned about rather than silently applied.

The argument-by-argument record of what is supported is `vignette("parameters")`.

## What the package is, rather than what changed

Earlier revisions of this file carried a full feature inventory here: the model
structure, the intervention modules, the interface, the validation summary and
the scope limits. That is a description of the package, not a changelog, and
keeping a second copy of it here is how it went stale. It lives where it is
maintained:

* **What is modelled, argument by argument** -- `vignette("model")`.
* **Where the mean field departs, and what to do about it** -- `vignette("using")`.
* **How well it agrees with the IBM** -- `vignette("comparison")`.
* **Scope, and the open discrepancies** -- the warning at the top of the README.

