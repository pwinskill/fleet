# Changelog

## fleet 0.0.0.9001

The defects a full package review found, fixed. Three of them move model
output: the age-band rendering weights, five intervention series that
took effect before their scheduled date, and TBV / mass-PEV age
targeting. The rest turn inputs that used to be accepted silently into
errors or warnings, switch the comparison CI on, and give the test suite
something to fail on. Output on bands that line up with the model age
grid, in a run with no interventions, is unchanged.

### Bug fixes

- **Output age bands are rendered by exact overlap, not by age-group
  midpoint.** A model age group used to contribute its whole population
  to a rendering band if its midpoint fell inside the band, and nothing
  at all otherwise. Two things went wrong with that. A band narrower
  than the group it lands in captures no midpoint, so the whole column
  family came back as hard zeros: 52 of the 65 single-year bands from 15
  to 80 years on the default grid, where `p_detect_lm_7300_7665 = 0`
  reads as “no malaria in 20-21 year olds” rather than “no data”. And a
  partial overlap counted the whole group or none of it, so the standard
  6-59-month band came back 2.2% light on a default parameter list, an
  error that flowed into every PfPR denominator and every postie
  person-day. A group now contributes the fraction of its own width
  lying inside `[min_age, max_age)`. The partition invariant is
  preserved and in fact strengthened: over a set of bands that partition
  a span of the age axis the weights sum to exactly 1 per group inside
  that span, so a group straddling a boundary is split between the two
  bands rather than handed whole to whichever one owns its midpoint.
  `fleet` now also warns when a rendering band overlaps no model age
  group at all, and reports the share of the population no band
  captured. **This moves published numbers for any band that does not
  align with the model age grid**, in every scenario, with or without
  interventions. Bands that do align are unaffected: regenerating the
  unit-test reference outputs, whose bands are grid-aligned, produced a
  byte-identical file.

- **Bed nets, IRS, carrying capacity, PEV and TBV no longer take effect
  before their scheduled date.** All five series are interpolated
  linearly, which is right for the within-round decay each of them
  carries, but it means a value that changes on a scheduled day ramps in
  from whatever knot precedes it, and the background grids are coarse:
  10 days for vector control, 30 for TBV and PEV, 5 for carrying
  capacity without seasonality. So an intervention arrived early, in
  some cases almost entirely. Nets scheduled for day 100 had reached 66%
  of their effect on EIR by day 99, before a single net existed. A TBV
  round on day 400 had delivered 88% of its transmission blocking by
  day 399. A mass PEV campaign whose efficacy onset was day 490 had
  delivered 34% of its FOI reduction by day 485. A carrying-capacity
  halving at day 100 started moving EIR on day 99. Each grid now carries
  a knot at `onset - 1` as well as at every scheduled day, pinning the
  pre-deployment value one day out and confining the ramp to the single
  step onto the scheduled day, so a round takes effect within a day of
  its schedule: the resolution the chemoprevention pulses already had.
  Every scenario using one of these five interventions moves: impact now
  starts on schedule instead of up to a grid step early, so the days
  around each deployment shift later.

- **TBV and mass-PEV age targeting weights by the fraction of each group
  covered.** Both used to select the model age groups whose midpoint
  fell in the target, so a target that straddled no midpoint selected
  nothing and the whole campaign was a silent no-op. Above age 14 the
  default grid is 5-yearly and only the years 17, 22, 27 and so on are
  midpoints, so `set_tbv(ages = 18:20)` vaccinated nobody, as did a mass
  PEV band of `min_ages = 20 * 365`, `max_ages = 21 * 365`, which falls
  between the midpoints at 17 and 22 years. Each group now carries the
  fraction of its own width inside the target, which multiplies coverage
  exactly as an age-restricted campaign should; a target aligned with
  group edges gives weight 1 and reproduces the old result exactly. A
  target that overlaps no group at all now warns instead of doing
  nothing quietly.

  Fractional weights then exposed a second defect in mass PEV. The
  several `min_ages`/`max_ages` bands of a *single* campaign were folded
  in one at a time with the independent-protection rule
  `1 - (1-a)(1-b)`, which is the rule for separate campaigns, not for
  two bands of the same one. The midpoint rule could not expose it,
  because a group belonged to at most one band. Splitting one band into
  two at 1195 days, which should change nothing, moved a straddling
  group’s covered fraction from 0.7295 to 0.6834, a 6.3% shortfall. The
  bands of one campaign now add into a single covered fraction (capped
  at 1) before coverage is applied; independent combination is kept
  across campaigns, where it is the right rule.

- **Inputs that could not mean what they said are rejected rather than
  run.** Each of these used to be accepted, or to fail somewhere
  downstream with a message that never named the parameter at fault:

  - `parameters$n_heterogeneity_groups` of 1 or 2. Fewer than three
    Gauss-Hermite nodes do not integrate the log-normal biting
    distribution, they just move every human off it: at `n = 1` the
    single node sits at `zeta = exp(-s2/2) = 0.434`, so every human in
    the model experienced 43% of the EIR the output column reported, and
    nothing said so. For a genuine no-heterogeneity run set
    `parameters$enable_heterogeneity = FALSE`, which collapses to one
    stratum at `zeta = 1` exactly.
  - `parameters$bite_dedup` outside 0/1 (`TRUE`/`FALSE` still accepted).
    It is a switch between two hazard forms, not a dial, and a value of
    0.5 was silently a third thing.
  - `parameters$acquired_immunity_offset` outside `[0, 1]`.
  - `default_age_lower(max_age =)` below 20 years, which made the grid’s
    5-yearly section `seq(15, max_age - 5, by = 5)` run backwards and
    died with base R’s “wrong sign in ‘by’ argument”.
  - `timesteps` that is not a single finite whole number `>= 1`. This
    also catches the swapped-argument call
    `run_simulation_ode(parameters, timesteps)`, which used to die on
    `parameters$init_EIR` with “\$ operator is invalid for atomic
    vectors”.

- **Editing the parameter list after `set_equilibrium()` now warns.**
  `set_equilibrium()` stores its back-translation as
  `parameters$eq_params`, and `build_inputs()` merges that over the live
  translation, so the stored copy wins for every shared biological
  constant: an edit made afterwards (`p$du <- 10`) is a silent no-op.
  That contract stands, since `set_equilibrium()` is meant to be the
  last call on the list, but a violation of it must not be silent.
  `fleet` now compares the two and warns, naming the `malariasimulation`
  fields whose live translation no longer agrees with the frozen copy,
  so the message points at the field you typed rather than at its
  equilibrium alias. Values that merely took a different arithmetic
  route to the same number cannot trigger it, and a deliberate custom
  `eq_params` is called out in the message as expected.

### Continuous integration

- **All four workflows are live**, so something now checks the model
  without being asked. The triggers that `0.0.0.9000` left staged and
  commented are switched on: `comparison.yaml` runs on every push and
  pull request that can touch the model, plus a Monday 06:00 UTC
  schedule, and `figures.yaml` on any push or pull request touching the
  comparison data, the renderers, the theme or the shared constants.
  `R-CMD-check.yaml` and `pkgdown.yaml` were already running on every
  push and pull request. All four keep their `workflow_dispatch`
  trigger, so any of them can also be run on demand. The Monday run is
  the one with no substitute: it is the only thing that catches
  `malariasimulation` moving underneath the comparison, the one
  staleness cause that never shows up in a commit here.
  `.github/workflows/CI.md` is now a description of what runs, not a
  go-live checklist.

- **`R-CMD-check.yaml` gains a job that fails when the compiled model is
  stale.** `inst/odin/malaria_ode.R` is the model’s source of truth, but
  the artifact that actually gets compiled is the committed, generated
  `src/malaria_ode.cpp` (with `R/dust.R`, `R/cpp11.R`, `src/cpp11.cpp`
  and `inst/dust/`), and no ordinary build regenerates any of it. So an
  edit to the odin model without a following `odin2::odin_package(".")`
  left every runner compiling the *old* model, the tests passing, and
  `comparison/check_drift.R` reporting “nothing moved”, which reads as
  *this change was numerically inert* when it means *this change was
  never compiled*. The new job regenerates from `inst/odin` and fails if
  the working tree comes out dirty. It runs once, on Ubuntu, with the
  hard dependencies plus the generator and nothing else, so a
  `malariasimulation` build failure on a runner cannot turn the
  staleness guard red for a reason that has nothing to do with
  staleness.

### Testing and tooling

- **The suite pins absolute output, not just relationships between
  outputs.** It could not previously detect a wrong answer: the only
  severe-incidence assertion was a ratio, invariant to any scale factor,
  so doubling `sev_inc_a` passed every test.
  `tests/testthat/reference-values.csv` and
  `reference-interventions.csv` now pin recorded output to 1e-6: across
  the EIR grid, and for each intervention at the day before deployment,
  at onset and after, which is the shape that catches an effect arriving
  early. A change that moves a number now has to say so.

- **A missing `malariasimulation` fails the suite instead of quietly
  emptying it.** Nearly every `test_that()` block opens with
  `skip_if_not_installed("malariasimulation")`, so a machine where the
  package was missing or unloadable ran the handful of dependency-free
  blocks, skipped the rest, and reported success having tested almost
  nothing. `tests/testthat/setup.R` now errors, with `FLEET_ALLOW_SKIP`
  as the deliberate escape hatch.

- **`NAMESPACE` is roxygen-generated.** It was hand-written and held the
  only copy of `useDynLib(fleet, .registration = TRUE)`, so the first
  `devtools::document()` would have dropped it and left the package
  unable to load.

## fleet 0.0.0.9000

### API changes

- **[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)’s
  signature now matches
  [`malariasimulation::run_simulation()`](https://rdrr.io/pkg/malariasimulation/man/run_simulation.html).**
  It was
  `(timesteps, parameters, correlations, init_EIR, age_lower, n_eir, n_foim, n_eip, n_ph, n_phc, atol, rtol, step_size_max, odin_file)`:
  fourteen arguments, eleven of which the IBM does not have. It is now
  `(timesteps, parameters, correlations, tuning)`.

  Two changes, both breaking.

  **`init_EIR` is gone as an argument.** The target EIR is a model
  input, so it belongs on the parameter list, and
  [`malariasimulation::set_equilibrium()`](https://rdrr.io/pkg/malariasimulation/man/set_equilibrium.html)
  already puts it there as `parameters$init_EIR`, which is where the IBM
  reads it from too. Seed the list the way you would for a matched IBM
  run:

  ``` r

  p <- malariasimulation::set_equilibrium(p, init_EIR = 20)
  out <- run_simulation_ode(3650, p)
  ```

  Setting `parameters$init_EIR` by hand still works, but skips the
  `eq_params` that `set_equilibrium()` stores and `fleet` honours. This
  is numerically inert for anyone already calling `set_equilibrium()`:
  the whole comparison harness reproduces bit-identically through the
  new API.

  **Everything numerical moved into `tuning`**, a new
  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  object (also accepted as a plain named list of just the fields you
  want to change). That is `age_lower`, the five Erlang stage counts
  `n_eir`/`n_foim`/`n_eip`/`n_ph`/`n_phc`, the solver controls
  `atol`/`rtol`/`step_size_max`, and `odin_file`:

  ``` r

  run_simulation_ode(3650, p, tuning = list(rtol = 1e-6, step_size_max = 10))
  ```

  The dividing line is that nothing epidemiological sits outside the
  parameter list any more: every field of `tuning` is a
  numerical-approximation knob, and changing one should move the answer
  only by its own error term.

  [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  validates its fields, and an unrecognised name in a `tuning` list is
  an error rather than a silently ignored field; a mistyped `rtol` that
  quietly did nothing would read as the tolerance simply not mattering.
  Passing any of the removed arguments to
  [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
  gives a targeted error naming the call that replaces it, rather than
  R’s bare “unused argument”.

- **`get_epi_outputs()` is removed**
  ([\#3](https://github.com/pwinskill/fleet/issues/3)). It was a thin
  wrapper that called
  [`postie::get_rates()`](https://rdrr.io/pkg/postie/man/get_rates.html)
  and
  [`postie::get_prevalence()`](https://rdrr.io/pkg/postie/man/get_prevalence.html)
  and returned them in a list, plus argument-splitting machinery
  (`rates_args`, `prevalence_args`, `...`) that was more to learn than
  the two calls it replaced. malariasimulation has no such wrapper (IBM
  post-processing pipelines call postie directly) and `fleet`’s output
  table is malariasimulation-shaped, so those pipelines already work on
  a `fleet` run unchanged. Replace `epi <- get_epi_outputs(out)` with
  `postie::get_prevalence(out, diagnostic = "lm")` and
  `postie::get_rates(out)`; every postie argument is now passed to
  postie directly. `postie` remains in Suggests for the vignette and
  tests, but nothing `fleet` exports depends on it.

### Behaviour changes

- **Custom demography: mosquito sizing now replicates
  `set_equilibrium()`.** malariasimulation sizes the adult-mosquito
  population from the human equilibrium under its *default* exponential
  age structure (`set_equilibrium()` → `equilibrium_total_M()`; custom
  mortality never enters), so under `set_demography()` the IBM drifts to
  whatever transmission that density supports. `fleet` used to re-solve
  the equilibrium under the custom age structure and hold `init_EIR`
  exactly, so the same parameter list realised different transmission in
  the two models (EIR 20 in fleet against 13.8 in the IBM in the
  comparison scenario, PfPR(2–10) 0.55 against 0.49). It now reproduces
  the IBM’s `total_M` exactly (`ibm_total_M()`) and root-finds the EIR
  at which its own equilibrium under the custom demography has that
  density, and seeds there. That is still a fixed point, so there is no
  burn-in (the scenario now seeds at EIR 13.86, PfPR 0.489). Under the
  default demography nothing changes: fleet keeps its own sizing there,
  which is the IBM’s formula at fleet’s own equilibrium and lands ~1.5%
  above the IBM’s `total_M` (so the same exponential mortality expressed
  through `set_demography()` seeds at EIR 19.7 rather than 20, within
  the IBM’s replicate noise). `parameters$hold_init_EIR = TRUE` restores
  the previous behaviour. Every malariaverse site file uses
  `set_demography()`, so this affects all country work: across the
  63-country validation set this release moves fleet’s mean
  clinical-incidence excess over the IBM from +10.1% to +8.7% and the
  regression slope from 1.02 to 1.00 (r 0.980 → 0.982); severe from
  +9.4% to +8.8%.

- **Prophylaxis is an Erlang chain, not one exponential compartment.**
  The IBM applies the Weibull survival `W(t - t_drug)` to each treated
  person’s infection probability, so a treated cohort’s mean protection
  at lag `t` is exactly `W(t)`. `fleet` represented both prophylaxis
  compartments (`Ph` post-treatment, `Ph_c` chemoprevention) as a single
  exponential at the Weibull mean, which leaks protection early — for
  SP-AQ (shape 4.3, scale 38.1) only 42% are still protected at day 30
  against the Weibull’s 70% — and under-estimated seasonal SMC (under-5
  clinical reduction 40% vs the IBM’s 52% in the comparison article; now
  54% against 52%, within the IBM’s replicate range). Both are now
  Erlang chains moment-matched to the Weibull. The chemoprevention chain
  has `k = 1/CV²` stages (14 for SP-AQ, 15 for DHA-PQP; SP-AQ chain at
  day 30: 0.67). The post-treatment chain follows the exponential
  treated stage `Tr` while the IBM’s clock runs from the dose, so its
  mean is the integrated protection left after `Tr`
  (`mean_W - ∫exp(-r_T t) W(t) dt`) and its length matches the variance
  of the whole `Tr + Ph` sojourn: 16 for SP-AQ, 20 for DHA-PQP, and 1
  for AL, whose 10-day protection is already less variable than `Tr`
  itself (measured output-identical to a 20-stage chain at half the run
  time). Counts are capped at

  20. `run_simulation_ode(tuning = ode_tuning(n_ph =, n_phc =))`
      overrides them; 1 recovers the old behaviour. The equilibrium seed
      distributes the prophylaxis mass across the stages exactly
      (`solve_disease_block()` generalised), so runs still hold flat.
      Runs now record only the output variables rather than the full
      state, which keeps memory flat despite the larger state.

### Bug fixes

- **Severe and all-infection incidence are counted with the deduplicated
  probability, not the raw hazard.** This is the change in this release
  that moves the no-intervention numbers; the review-pass fixes in
  `0.0.0.9001` move intervention scenarios and misaligned rendering
  bands, and nothing else. `clin_inc_a` was already counted with `h_c`;
  `sev_inc_a` and `inc_a` used `FOI`. These outputs are rates that R
  integrates over a day, and the right rate depends on how fast the
  compartment drains. `S` and `U` drain at the full `FOI`, so
  `int FOI*X exp(-FOI t) dt = X*(1 - exp(-FOI)) = X*p_inf` — the IBM’s
  count, exactly. `A` does not: a sub-clinical re-infection of an `A`
  leaves them in `A`, so `A` drains only at `h_c` and stays roughly flat
  over the day, making `int FOI*A dt ~ FOI*A` and over-counting by
  `FOI/p_inf = -log(1-p_inf)/p_inf`. `A` now takes `p_inf` directly. The
  IBM draws severe from the *deduplicated* `infections` bitset
  (`update_severe_disease`), so this is a replication defect, not a
  modelling choice.

  Measured against 3 IBM replicates of 10,000 people on the same
  parameter list: all-age `n_inc_*` ran **+2.6% (EIR 20) and +4.2% (EIR
  50)** above the IBM median, both outside its replicate spread; after
  the fix, **+0.6% and +0.4%**. The bias was concentrated in adults,
  where `A` is a large share of the at-risk pool (under-5 `n_inc_*` was
  only +0.3% / +1.7%), the signature of an `A`-only defect. Clinical
  incidence and prevalence are unchanged.

  It also makes one comparison worse, and that is worth stating plainly.
  The same correction lowers severe incidence slightly, and fleet
  already ran below the IBM there: all-age severe goes from −5.9% to
  −6.2% at EIR 20 and from −3.4% to −4.0% at EIR 50, dropping just
  outside the IBM’s 10–90% band at those two EIRs (by 0.2% and 0.3% of
  the lower edge) where it had been marginally inside. The upward bias
  had been masking part of a larger, separate severe deficit. That
  compensation was accidental, so the corrected numbers are the honest
  ones and the severe gap is now the largest open discrepancy between
  the two models.
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)
  and its tables are re-rendered accordingly; the comparison harness
  does not render all-infection incidence, so those figures show this
  change’s cost without showing its benefit.

- **`test-prophylaxis-chain.R`’s first test no longer errors when
  `malariasimulation` is absent.** It read `AL_params`/`SP_AQ_params`
  without a `skip_if_not_installed()` guard, the only one of the suite’s
  57 test blocks to do so. `malariasimulation` is a GitHub Remote in
  Suggests, so on any check machine that could not install it the other
  56 blocks skipped and this one errored, turning a missing optional
  dependency into an `R CMD check` failure. The block is now split
  rather than simply guarded: its closed-form `erlang_stages()` checks
  need no drug tables and keep running unguarded, so the suite retains
  real assertions without the dependency, and only the Weibull
  parameters read out of the IBM’s tables sit behind the guard.

- **`build_inputs()` no longer recomputes grid-invariant PEV
  protection.** Two fixes in `pev_series()`, both verified bit-identical
  ([`identical()`](https://rdrr.io/r/base/identical.html) on the full
  `[n_age, n_time]` multiplier matrix, across the EPI single-timestep
  path, the EPI multi-timestep path and mass campaigns).
  `.pev_eff_fun()` built its interpolator with
  [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html), which
  re-runs `regularize.values()` (a sortedness check and a copy of the
  whole daily grid) on *every* scalar call; it now builds one
  [`stats::approxfun()`](https://rdrr.io/r/stats/approxfun.html)
  closure. And the EPI loop called `pev_protection()` once per grid
  point when `.booster_cov_vec()`’s admin-time argument is ignored
  (which the single-row `booster_coverage` that `set_pev_epi()` builds
  by default always is), so that call is now hoisted to once per age
  band behind a guard mirroring `.booster_cov_vec()`’s own. Together:
  `pev_series()` on a 20-year two-booster EPI schedule drops from 2.92 s
  to 1.09 s (2.7x).

- **Chemoprevention pulses renew existing protection and clear
  `Tr_slow`.** The IBM’s `update_mass_drug_admin()` resets `drug_time`
  for everyone successfully treated whatever their state; `fleet`’s
  pulse left people already in `Ph`/`Ph_c` decaying from their earlier
  dose and skipped the slow-clearance treated compartment. Both now move
  to the first stage of `Ph_c` with the rest of the covered fraction.
  That is material for monthly SMC rounds, where most of the previous
  round’s recipients are still protected.

- **`rT_slow` uses the whole-day exit probability.** The
  slow-parasite-clearance rate under antimalarial resistance was
  `1/dt_slow`; it is now `1 - exp(-1/dt_slow)`, the same conversion
  `build_inputs()` applies to `rA`/`rD`/`rU`/`rT` (the IBM leaves states
  with per-day probability `rate_to_prob(1/d)`).

- **Custom demography: top age group’s death rate.** The open-ended
  oldest model age group is represented by its lower bound, which the
  right-closed `set_demography()` bins assigned to the band *below* it
  whenever that bound coincided with a bin edge: 80 y on the default
  grid, a very common edge. The over-80s therefore died at the 60–80
  rate (0.05 instead of 0.12 per year in the comparison scenario) and
  the group held ~2.4× too many people (60–85 share 18.2% against an
  analytic 13.1%). The death rate is now looked up one day above the
  group’s lower bound. Found by the new comparison article.

### Documentation

- README gains a **Run times** section
  ([\#1](https://github.com/pwinskill/fleet/issues/1)): seconds per
  complete
  [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
  call across six scenarios and three horizons, plus the
  population-independence check (flat from 1,000 to 10,000,000 people)
  and what the solver presets cost. Reproduced by
  `comparison/benchmark.R`, which reports the minimum of five repeats;
  contention can only add time, so the fastest repeat is the least
  contaminated estimate.

- The intervention-impact figure gains two panels (clinical incidence at
  all ages, and severe incidence at all ages) alongside the existing
  PfPR(2–10) and under-5 clinical panels, now laid out 2×2. This needed
  all-age rendering bands in `comparison/run_replicates.R` and a re-run
  of the five intervention scenarios. Together the four panels show
  *who* each intervention protects: SMC and RTS,S roughly halve their
  effect when measured over all ages rather than under-5s, while nets
  and IRS read the same either way. fleet lands inside the IBM’s 10–90%
  replicate range in 16 of the 20 scenario × outcome cells; the noisiest
  is treatment scale-up on severe incidence, where the IBM’s own
  replicates span −12% to +6% and the two models differ in sign.

- **The 63-country site comparison is now an explicit snapshot.** It was
  the one figure whose inputs live outside this repository: a roughly
  seven-hour run against the malariaverse site files, in a separate,
  unversioned checkout. Both the figure and its statistics are therefore
  frozen rather than kept in step, and both now say so — in the README
  caption, in a callout above the figure in
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md),
  and in `comparison/data/site_snapshot.json`, which records when the
  snapshot was taken and against which `fleet` version.

  This also closes a hole. `summary_tables.R` used to compute those
  statistics live and silently omit the whole section when the
  validation checkout was absent, which is the normal case for anyone
  else and for CI. That would have dropped the article’s site numbers
  out of the reproducible record, and failed the new figures job for a
  reason unrelated to the model. It now reads the committed snapshot.
  `render_figures.R` likewise leaves the committed PNG alone and says
  which flag redraws it. `CMP_REFRESH_SITES=1` re-takes both, and is the
  only thing that touches either.

  The alternative — wiring the validation run into the harness — was
  considered and rejected: it would put a seven-hour dependency in the
  path of a two-minute check, for a comparison that is meant to be
  repeated at a release rather than on every change.

- **CI for the comparison.** Two new workflows, `comparison.yaml` (the
  drift check) and `figures.yaml` (do the committed figures and tables
  still match the committed data), plus `dependabot.yml` for the Actions
  themselves. Both workflows landed here carrying only a
  `workflow_dispatch` trigger, with their real triggers staged directly
  above as comments, and `0.0.0.9001` switched them on: all four
  workflows are live now, as `.github/workflows/CI.md` describes. Read
  that file for what each one runs on and why each piece is shaped the
  way it is.

  The design choices that matter for keeping it low-maintenance.
  `paths-ignore` rather than `paths` wherever the question is “could
  this break the model”, so a new source directory is covered by default
  instead of being silently dropped from an allow-list.
  `concurrency: cancel-in-progress`, so pushing again supersedes an
  in-flight run. Dependabot monthly and grouped into one PR, because
  runner images and actions are deprecated on GitHub’s schedule and the
  failure mode otherwise is a red build one morning for a reason
  unrelated to the package. And `figures.yaml` fails on `tables.md` but
  only warns on the PNGs, since a PNG can differ byte-for-byte from font
  hinting alone and a check that cries wolf is one you learn to ignore.

  `R-CMD-check.yaml` gains `paths-ignore` and `concurrency` blocks in
  the same staged form, left commented here and switched on in
  `0.0.0.9001` with the rest.

- **New `comparison/check_drift.R`: check the match without re-running
  the IBM.** The IBM does not depend on `fleet`, so its committed rows
  stay valid for any fleet-side change, and there was never a reason to
  re-run a 25-minute IBM sweep to find out whether the match still
  holds. The check re-runs fleet alone (~2 min for all 18 scenarios, ~20
  s for a `CMP_ONLY` subset) and compares against the frozen reference.

  It answers two questions separately, because they mean different
  things. **Did anything move?** fleet now against fleet’s committed
  rows, to 1e-6. A moved number is not automatically wrong — a
  deliberate model fix moves numbers — but it must be seen; silent
  movement is how a regression ships. **Is the match still good?** fleet
  now against the committed IBM medians and 10–90% bands, at per-outcome
  thresholds set from what the models achieve today with headroom (max
  \|fleet − IBM\| of 2% for PfPR, 5% for both clinical measures, 10% for
  severe; and a tolerated count of grid points falling outside the
  replicate band, expressed as failures allowed rather than successes
  required so a subset run does not make the limit unreachable). Only
  the second fails the run; `CMP_STRICT=1` fails on any movement too.

  Supporting this, `comparison/data/ibm_reference.json` now records what
  the committed IBM rows were actually made from: date,
  malariasimulation version, `N_REP`/`POP`/burn-in, and a digest of the
  scenario definitions — every parameter that differs from a bare
  `get_parameters()`, plus each EIR and horizon, deliberately not the
  whole list, which carries defaults that shift between
  malariasimulation versions for unrelated reasons. `run_replicates.R`
  writes it whenever it actually runs the IBM. The check warns when the
  installed malariasimulation or the scenario digest no longer matches,
  which is the only time the expensive re-run is genuinely required — so
  the script says when, rather than leaving it to memory.

  The scenario definitions, the shared summariser and the fleet-run loop
  moved into `comparison/scenarios.R`, sourced by both scripts.
  Otherwise a drift check could quietly test a different set of
  scenarios, or a different solver tuning, from the one the reference
  was built on.

- **The model flow diagram is simplified.** It had accumulated more
  visual grammar than it could carry, and the worst of it was that the
  *same* layered idiom meant two unrelated things: offset layers behind
  a panel meant “array dimension”, while a stacked outline behind a box
  meant “Erlang chain”. The layering now means one thing. The deck
  behind the human panel carries one dimension on each visible edge: age
  *i* = 1 … 52 and biting heterogeneity *j* = 1 … 5, so the layers
  illustrate the grid and nothing else.

  Three elements are gone. The stacked outlines on `P`, `P_c` and the
  EIP, with their three stage-count labels: the chains are disclosed in
  the caption and
  [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
  instead, where the counts belong anyway since they vary by drug. The
  dotted at-risk rings on `S`, `A` and `U`, and the key row explaining
  them: a whole grammar element, plus a padding offset on every flow
  that touched those boxes, for one fact now stated in words beside the
  hazard. And the immunity box with its coupling arrow, which sat in the
  busiest part of the panel to represent something algebraic rather than
  a state. The key is down from six rows to five, and the captions in
  the README and
  [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
  are rewritten to match; they described idioms the figure no longer
  uses.

  Also fixed while in there: the FIDELITY invariant cited two odin line
  numbers that had rotted by ~18 lines, so it now cites the equations
  (`deriv(A)`, `deriv(U)`) instead. Bare line numbers into an actively
  edited file rot, and these already had.

- New figure `cmp_programme_ts` and a *Programmes over fifteen years*
  section in
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md).
  Every other intervention figure isolates one builder over six years,
  which is the shape for attributing a difference but not for seeing
  what a programme does. Five new harness scenarios (`ts_*`) run 15
  years past deployment at EIR 20 in the seasonal setting, all carrying
  20% baseline case management so each is the one before it with one
  more thing added: no interventions, bed nets alone (5 campaigns,
  3-yearly), seasonal SMC alone (4 rounds a year for 15 years), case
  management alone (20% → 60%), and all three together. Four outcomes
  per scenario: PfPR(2–10), clinical incidence in under-5s and over all
  ages, and all-age severe incidence.

  Through five net distributions, sixty SMC rounds and a treatment
  scale-up the two models stay together: fleet sits inside the IBM’s
  10–90% replicate band in 78–84% of the 915 scenario-months per outcome
  (an 80% band contains a perfectly-tracking deterministic mean about
  80% of the time, so that is at the target, not short of it), and the
  largest disagreement in 15-year mean burden reduction across the 16
  scenario × outcome cells is 1.8 percentage points. The figure needs
  four separate faceted columns rather than one `facet_grid`, because
  `facet_grid` frees the y scale by row and here the scales differ by
  column. `summary_tables.R` reports both statistics; it deliberately
  does *not* report a per-month relative error, because the seasonal
  dry-season trough goes to within rounding of zero and \|fleet − IBM\|
  / IBM then reaches billions of percent on months carrying no burden.

  The IBM envelope on this figure is drawn heavier than on the rest of
  the set (0.40 against `ENV_ALPHA`’s 0.16, and the lines are thinned to
  match), because fifteen years of monthly points is about three pixels
  per month and at 0.16 the band was invisible even where it is wide. It
  has real width to disclose: over the months carrying the top quartile
  of burden the 10–90% replicate range is 26% of the panel peak for
  all-age severe incidence, against 3–6% for prevalence and the two
  clinical measures. Severe is the rarest of the four, so a 30-day bin
  holds only ~20 severe episodes at the seasonal peak in a population of
  10,000. The caption and the article now give those widths, since the
  figure can only hint at them. Raised locally rather than in `theme.R`,
  so the four already-reviewed figures are not changed unseen.

- The equilibrium-vs-EIR figure gains the same two panels, also laid out
  2×2: all-age clinical incidence and all-age severe incidence against
  EIR, beside the existing PfPR(2–10) and under-5 clinical panels. The
  all-age views carry the *shape* of the relationship rather than its
  level: clinical incidence spans 3.8-fold across the EIR grid against
  18-fold in under-5s, and severe incidence is not monotonic at all,
  plateauing around EIR 20–50 in both models before falling by EIR 120.
  Agreement is looser here than on the under-5 metrics and leans the
  other way: fleet runs −1.6% to +2.1% of the IBM median on all-age
  clinical (inside the replicate band at five of six EIRs) and −5.9% to
  +0.1% on all-age severe (inside the band at all six), so the small
  excess fleet carries in under-5s is more than repaid in the 5–20 year
  bands. `summary_tables.R` reports both new outcomes over the grid.

- Three stale comment blocks are removed, each of which had been left
  sitting *above* the corrected version that replaced it, so the wrong
  one is the one a reader hits first. `R/interventions.R` documented
  `drug_mix()`’s weights as each drug’s peak coverage (now only the
  `t = NULL` fallback) and `rP` as `mean_W - 1/rT`, the exact
  approximation the comment 33 lines above it exists to say is wrong
  (19% error for AL). `inst/odin/malaria_ode.R` described `beta_eff` as
  a per-species interpolated series that varies under nets/IRS; it is a
  plain constant `parameter()`, as the comment 17 lines below it already
  said. `tests/testthat/test-package.R` called `rT_slow` a reciprocal
  time three lines above the comment (and the assertion) that correctly
  make it a whole-day exit probability.

- Every line figure now draws the IBM’s dashed median *over* fleet’s
  solid line rather than under it. Where the models agree — which, on
  these figures, is nearly everywhere — the line drawn last is the only
  one visible, so the old order quietly hid the IBM and made a
  two-series panel read as one. Dashes let the solid line show through
  the gaps, so both read. Points keep the opposite order for the same
  reason: fleet’s hollow marker hides less than the IBM’s filled one, so
  it goes on top. The rule is written down in `comparison/theme.R` and
  applies to `core_eir`, `core_age`, `core_demography`, `core_seasonal`
  and `int_timeseries`; the effect is largest on the low-noise
  intervention panels, which previously looked like a single curve. An
  unused `geom_ibm_envelope()` helper that bundled the old order is
  removed.

- The country-site comparison figure reads its 1:1 correlation properly
  ([\#2](https://github.com/pwinskill/fleet/issues/2)). Two changes. The
  reference line is a plain orange dash, hue-opposed to the indigo ramp,
  and not a series colour, so it cannot be misread as a model. And the
  low end of the hex fill ramp is now near-white: the counts are wildly
  skewed, with 52% of the cells carrying 0.2% of the sub-site-months
  while the top 5% of cells carry 90% of them, so the old saturated low
  end spent most of the plot’s ink on almost none of the data and buried
  the ridge it was meant to show. The scale is still log10, so sparse
  cells stay visible (the scatter is real and worth seeing); they simply
  no longer out-shout the ridge. The legend gained intermediate breaks
  (1 / 10 / 100 / 1,000 / 10,000) so the ramp can be decoded.

- New article
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)
  — *Comparison with malariasimulation* — with a redesigned figure set:
  equilibrium PfPR and clinical incidence across EIR 1–120, age profiles
  of prevalence / clinical / severe incidence, the settled seasonal
  cycle, custom demography, the 63-country monthly site-file comparison,
  and five interventions (treatment scale-up, RTS,S via EPI, seasonal
  SMC, IRS, bed nets) as time series and a %-reduction summary. The IBM
  is drawn as the median of 10 replicates with a 10–90% band rather than
  as a single realisation. `comparison/` was rewritten around one runner
  (`run_replicates.R`), one renderer and a shared theme; the
  per-replicate summaries are committed so the figures can be redrawn
  without re-running the models. The pre-replication-pass figures
  (`A_`–`F_`, `inc_*`) and the scripts that made them are removed.

- The comparison surfaced three fleet-side items, all addressed in this
  release (see *Behaviour changes* and *Bug fixes* above): the top age
  group’s death rate under `set_demography()`, exponential prophylaxis
  under-estimating seasonal SMC, and the `set_equilibrium()` convention
  under custom demography.

### Exact-replication pass (malariasimulation v3.0.0)

A systematic audit replaced every place `fleet` approximated a mechanism
whose malariasimulation implementation is known. No tuned constants were
introduced.

- **Output bands — bug fix, affects every user.** `render_output()`
  emitted `n_inc_clinical_*`/`n_inc_severe_*`/`n_inc_*` over the *union*
  of all rendering bands plus hardcoded extras, producing OVERLAPPING
  strata (e.g. a 2-10y prevalence band alongside a 0-5/5-15/15-100
  clinical partition, and the 2-10y band twice). `postie` treats each
  column as an independent stratum, so any person-day-weighted aggregate
  double-counted the overlapped ages, inflating all-age clinical by
  ~1.21x and severe by ~1.10x. Each family is now emitted only over its
  own `*_rendering_ages`, exactly as malariasimulation does, with
  `n_age_*` over the union.
- **Maternal immunity source band.** `which.min(abs(age_mid - 20*365))`
  hit an exact tie on the default grid (17.5y and 22.5y are equidistant
  from the 20y band edge) and silently took the first, drawing
  `ICM`/`IVM` from 15-20 year-olds instead of the 20-21y mothers the IBM
  uses; they ran 9-15% low. Now selects the band containing 20y.
- **EIP survival.** The Erlang chain carried mortality in every stage,
  giving through-survival `(reip/(reip+mum))^n_eip` instead of the IBM’s
  `exp(-mum*dem)` (+4% at baseline, worsening as vector control raises
  `mum`). The chain is now a loss-free delay with the survival applied
  at the exit, matching `adult_mosquito_eqs.cpp`; `total_M` needs no
  compensation.
- **Slow parasite clearance.** The IBM assigns each treated individual
  to `dt` or `dt_slow` by a Bernoulli draw — a two-component mixture.
  `fleet` collapsed this into one exponential at the blended mean; `Tr`
  is now split into parallel `Tr`/`Tr_slow`.
- **Clinical event counting.** The IBM resolves at most one infection
  outcome per person per day, so clinical episodes are `phi*p*N`. Using
  `phi*FOI` over-counted (the ODE re-exposes within the day); the
  clinical hazard is now `-log(1 - phi*p)`, which integrates to exactly
  `phi*p*N`.
- **State sojourns.** The IBM exits states with per-day probability
  `1-exp(-1/d)`, so its realised mean dwell is `1/(1-exp(-1/d))` (5.52 d
  for `dd = 5`, not 5). Rates now use that form.
- **Immunity boosting.** Boost rate is now `q/(q*u_eff + 1)` with the
  deduplicated per-day event probability `q = 1-exp(-rate)` and the
  integer refractory window `u_eff = ceil(u) - 1`: the IBM boosts on a
  per-day event and tests `(timestep - last_boosted) >= u` with an
  integer timestep. `ICA`/`ID`/`IVA` are boosted only for individuals
  eligible to be infected (`S`/`A`/`U`), as in the IBM; `IB` is boosted
  for everyone bitten.
- **PEV antibodies.** Efficacy is integrated over the per-individual
  antibody distribution the IBM samples (4-D Gauss-Hermite over
  `cs`/`rho`/`ds`/`dl`), instead of being evaluated at the profile
  median — the latter overstated R21 efficacy by up to ~4.7
  pp. Validated against a Monte-Carlo of the IBM’s own sampler to
  \<5e-4.
- **PEV gating.** malariasimulation never reads `parameters$pev`; it
  gates on `pev_epi_coverages`/`pev_epi_timesteps`/`mass_pev_timesteps`.
  `fleet` gated on `pev$pev`, so the same list could mean “PEV on” to
  one model and “off” to the other.
- **Seeding consequence.** Because `malariaEquilibrium` encodes the
  simplified forms, the analytic seed is now a close approximation
  rather than an exact fixed point: an undisturbed run relaxes by up to
  ~0.5% over the first years, as malariasimulation itself does. See the
  README “Seeding” note.

#### Validation after this pass

Monthly, *P. falciparum* only on both sides, 63 countries / 1,391
sub-sites / 450,684 sub-site-months, against pre-run malariasimulation
output:

| metric         | before the pass | after the pass | now (snapshot) |
|----------------|-----------------|----------------|----------------|
| clinical slope | 1.16            | 1.02           | **0.997**      |
| clinical r     | 0.982           | 0.980          | **0.982**      |
| clinical bias  | +0.070          | **+0.030**     | not re-taken   |
| clinical RMSE  | 0.125           | **0.084**      | not re-taken   |
| severe slope   | 1.02            | 0.967          | **0.946**      |
| severe r       | 0.935           | 0.954          | **0.958**      |

The “now” column is the committed snapshot,
`comparison/data/site_snapshot.json`, taken on 2026-09-09 against the
whole of this release rather than against this pass alone: the
prophylaxis chains and the mosquito-sizing convention landed after it
(see *Behaviour changes*) and moved the clinical slope from 1.02 to
0.997 and clinical r from 0.980 to 0.982. Bias and RMSE in absolute
incidence units were measured at the pass and are not carried in the
snapshot, which records relative bias instead.

Per-site median monthly clinical r = 0.99. Note the “before the pass”
column was measured through the overlapping-band artefact, which
inflated fleet’s side by ~1.21x (clinical) and ~1.10x (severe); the
apparently unbiased severe slope of 1.02 was two opposite-signed errors
cancelling.

Known remaining discrepancies, all characterised rather than tuned away:

- **A positive offset against the site files.** fleet runs 8.7% above
  the IBM on average for monthly clinical incidence and 8.8% for severe,
  in the snapshot above. With the slopes at 1.00 and 0.95 that excess is
  an intercept, a few hundredths of an episode per person-year, so it
  matters in low-incidence sub-sites and dry-season months and hardly at
  all where incidence is high. Neither the prophylaxis chains nor the
  mosquito-sizing convention explains more than about a point of it (see
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)),
  and the rest is not yet understood.
- **Severe incidence on the controlled EIR grid runs the other way**,
  which is worth holding alongside the site-file offset rather than
  instead of it. On the EIR 1-120 equilibrium sweep in
  `comparison/data/rep_eq.csv`, all-age severe incidence sits within
  0.8% of the IBM median at EIR 1, 3 and 10, and then falls away over
  the plateau: -6.2% at EIR 20, -4.0% at EIR 50 and -4.6% at EIR 120. By
  age band it is 15-25% low over 5-20 years. The cause is mean-field
  averaging of a convex Hill function, described under *Diagnostic
  conventions* in the README.
- **Stochastic elimination.** Where the IBM’s transmission dies out
  (e.g. Guinea-Bissau, ratios 2.3-2.9 with r ~ 0.53), a deterministic
  mean field cannot follow it to zero. This is structural, not a defect.
- **Very low transmission (EIR ~ 1)** over-predicts against the pre-run
  files, but matches a *live* malariasimulation 3.0.0 run at the same
  site — check the baseline before reading this as a model error.
- **Very high transmission (EIR \> 200)** under-predicts (~0.78-0.91).

First development release: a deterministic mean-field (ODE) twin of the
`malariasimulation` individual-based model of *P. falciparum* malaria,
built on odin2/dust2.

### Parameter-ingestion completeness

A systematic audit of every user-facing `malariasimulation` config
function (`get_parameters` + all `set_*`) closed these gaps so the model
uses the full parameter flexibility they expose:

- **PEV** now models a full **booster sequence** (not just the first
  booster), with per-booster profiles, spacing, and (conditional)
  coverage, for both `set_pev_epi` and `set_mass_pev`; `min_wait` is
  honoured.
- **Mass PEV** applies **every** target age band at **every** campaign
  (multiple `min_ages`/`max_ages` are no longer collapsed to the
  campaign count).
- **Chemoprevention** (MDA/SMC/PMC) targets age bands by **fractional
  overlap** with the model age groups, so narrow bands (e.g. PMC dose
  ages) are no longer silently dropped and coarse groups are no longer
  over-treated. Co-deployed chemoprevention types share a
  coverage-weighted prophylaxis duration.
- **Antimalarial resistance** (early-treatment-failure) now applies to a
  drug used for **chemoprevention** (SMC/MDA/PMC), not only clinical
  treatment.
- **IRS** now **accumulates across spray rounds** (sprayed protection
  never expires in the IBM), matching the bed-net accumulation.
- **TBV** targets the exact `ages` set (non-contiguous ages honoured).
- `set_equilibrium`’s custom **`eq_params`** are read explicitly.
- New **all-infection incidence** output (`n_inc_*`;
  `incidence_rendering_*`).
- Remaining mean-field approximations (time-varying multi-drug
  first-line switch; seasonal PEV boosters) are **warned**, not silently
  applied.

### What the package is, rather than what changed

Earlier revisions of this file carried a full feature inventory here:
the model structure, the intervention modules, the interface, the
validation summary and the scope limits. That is a description of the
package, not a changelog, and keeping a second copy of it here is how it
went stale. It lives where it is maintained:

- **What is modelled, argument by argument** –
  [`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md).
- **Where the mean field departs, and what to do about it** –
  [`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md).
- **How well it agrees with the IBM** –
  [`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md).
- **Scope, and the open discrepancies** – the warning at the top of the
  README.
