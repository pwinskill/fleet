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

