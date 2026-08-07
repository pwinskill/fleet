# blink 0.0.0.9000

## Exact-replication pass (malariasimulation v3.0.0)

A systematic audit replaced every place `blink` approximated a mechanism whose
malariasimulation implementation is known. No tuned constants were introduced.

* **Output bands — bug fix, affects every user.** `render_output()` emitted
  `n_inc_clinical_*`/`n_inc_severe_*`/`n_inc_*` over the *union* of all rendering
  bands plus hardcoded extras, producing OVERLAPPING strata (e.g. a 2-10y prevalence
  band alongside a 0-5/5-15/15-100 clinical partition, and the 2-10y band twice).
  `postie` treats each column as an independent stratum, so any person-day-weighted
  aggregate double-counted the overlapped ages — inflating all-age clinical by ~1.21x
  and severe by ~1.10x. Each family is now emitted only over its own
  `*_rendering_ages`, exactly as malariasimulation does, with `n_age_*` over the union.
* **Maternal immunity source band.** `which.min(abs(age_mid - 20*365))` hit an exact
  tie on the default grid (17.5y and 22.5y are equidistant from the 20y band edge) and
  silently took the first, drawing `ICM`/`IVM` from 15-20 year-olds instead of the
  20-21y mothers the IBM uses; they ran 9-15% low. Now selects the band containing 20y.
* **EIP survival.** The Erlang chain carried mortality in every stage, giving
  through-survival `(reip/(reip+mum))^n_eip` instead of the IBM's `exp(-mum*dem)`
  (+4% at baseline, worsening as vector control raises `mum`). The chain is now a
  loss-free delay with the survival applied at the exit, matching
  `adult_mosquito_eqs.cpp`; `total_M` needs no compensation.
* **Slow parasite clearance.** The IBM assigns each treated individual to `dt` or
  `dt_slow` by a Bernoulli draw — a two-component mixture. `blink` collapsed this into
  one exponential at the blended mean; `Tr` is now split into parallel `Tr`/`Tr_slow`.
* **Clinical event counting.** The IBM resolves at most one infection outcome per
  person per day, so clinical episodes are `phi*p*N`. Using `phi*FOI` over-counted
  (the ODE re-exposes within the day); the clinical hazard is now
  `-log(1 - phi*p)`, which integrates to exactly `phi*p*N`.
* **State sojourns.** The IBM exits states with per-day probability `1-exp(-1/d)`, so
  its realised mean dwell is `1/(1-exp(-1/d))` (5.52 d for `dd = 5`, not 5). Rates
  now use that form.
* **Immunity boosting.** Boost rate is now `q/(q*u_eff + 1)` with the deduplicated
  per-day event probability `q = 1-exp(-rate)` and the integer refractory window
  `u_eff = ceil(u) - 1` — the IBM boosts on a per-day event and tests
  `(timestep - last_boosted) >= u` with an integer timestep. `ICA`/`ID`/`IVA` are
  boosted only for individuals eligible to be infected (`S`/`A`/`U`), as in the IBM;
  `IB` is boosted for everyone bitten.
* **PEV antibodies.** Efficacy is integrated over the per-individual antibody
  distribution the IBM samples (4-D Gauss-Hermite over `cs`/`rho`/`ds`/`dl`), instead
  of being evaluated at the profile median — the latter overstated R21 efficacy by up
  to ~4.7 pp. Validated against a Monte-Carlo of the IBM's own sampler to <5e-4.
* **PEV gating.** malariasimulation never reads `parameters$pev`; it gates on
  `pev_epi_coverages`/`pev_epi_timesteps`/`mass_pev_timesteps`. `blink` gated on
  `pev$pev`, so the same list could mean "PEV on" to one model and "off" to the other.
* **Seeding consequence.** Because `malariaEquilibrium` encodes the simplified forms,
  the analytic seed is now a close approximation rather than an exact fixed point: an
  undisturbed run relaxes by up to ~0.5% over the first years, as malariasimulation
  itself does. See the README "Seeding" note.

First development release: a deterministic mean-field (ODE) twin of the
`malariasimulation` individual-based model of *P. falciparum* malaria, built on
odin2/dust2.

## Model

* Age- and biting-heterogeneity-structured human model (`S`/`D`/`A`/`U`/`Tr`
  plus treatment (`Ph`) and chemoprevention (`Ph_c`) prophylaxis compartments)
  with the full immunity system (`IB`/`ICA`/`ICM`/`ID`/`IVA`/`IVM`), coupled to
  the compartmental mosquito model (`E`/`L`/`P`/`Sm`/EIP-chain/`Im`) per species.
* Seeded at the `malariaEquilibrium` fixed point; holds flat at equilibrium to
  machine precision, including under clinical treatment (`ft > 0`), via a
  corrected prophylaxis-aging recursion.
* Lags (human EIR/FOIM, mosquito EIP) implemented as tunable Erlang chains,
  exact at equilibrium for any stage count.
* Infection hazard reproduces the IBM's **per-timestep bite deduplication**: the
  IBM collects the day's bitten individuals in a bitset, so a person bitten
  repeatedly in one timestep is infected at most once. The daily infection
  probability is therefore `(1 - exp(-EPS)) * b`, converted to a hazard with
  `-log(1 - p)`. This saturates, where the unbounded `b * EPS` over-predicts
  infection as exposure approaches one bite/person/day — i.e. at seasonal peaks
  and in high-`zeta` strata. Across five countries (36,936 sub-site-months) it
  moves the monthly clinical regression slope vs the IBM from 1.14 to 1.10, with
  the effect correctly gated by exposure: largest in the highest-EIR setting
  (BFA slope 1.10 -> 1.04) and nil at low transmission (MMR unchanged).
  `parameters$bite_dedup = 0` restores the linear form, for which the
  `malariaEquilibrium` seed is an exact fixed point.
* `parameters$acquired_immunity_offset` exposes the IBM's `+0.5` acquired-immunity
  offset in the `b`/`phi`/`theta` Hill functions. Default `0`: an A/B against a
  live malariasimulation 3.0.0 ensemble mean showed `0` matches clinical (~2.4x
  closer) and severe (~9x closer) incidence better than `0.5`, because the offset
  is a per-individual detail that does not carry over to a stratum mean.

## Interventions

* Clinical treatment (time-varying coverage) and antimalarial resistance
  (early-treatment-failure and slow-parasite-clearance arms, coverage-weighted
  across drugs).
* Bed nets and indoor residual spraying (mean-field vector control).
* Chemoprevention: MDA, SMC and PMC, applied as pulsed mass drug administration.
* Vaccines: pre-erythrocytic (RTS,S / R21; EPI and mass, with boosters and
  time-varying EPI coverage) and transmission-blocking (TBV).
* Seasonality (Fourier rainfall driving larval carrying capacity) and flexible
  carrying capacity.
* Custom demography (`set_demography`): **time-varying** age-specific mortality
  (`mu_age(t)` interpolated over `deathrate_timesteps`) and the resulting
  equilibrium age structure.
* Bed nets support both exponential (log-uniform) and **logistic** net retention,
  and accumulate repeated distributions as the full mean-field usage mixture
  (most-recent-receipt x retention survival, per-distribution net-age decay) —
  needed to run real `site::site_parameters()` inputs.
* `run_simulation_ode()` exposes the dust2 solver controls (`atol`, `rtol`,
  `step_size_max`) for faster long dynamic projections.

## Parameter-ingestion completeness

A systematic audit of every user-facing `malariasimulation` config function
(`get_parameters` + all `set_*`) closed these gaps so the model uses the full
parameter flexibility they expose:

* **PEV** now models a full **booster sequence** (not just the first booster),
  with per-booster profiles, spacing, and (conditional) coverage — for both
  `set_pev_epi` and `set_mass_pev`; `min_wait` is honoured.
* **Mass PEV** applies **every** target age band at **every** campaign (multiple
  `min_ages`/`max_ages` are no longer collapsed to the campaign count).
* **Chemoprevention** (MDA/SMC/PMC) targets age bands by **fractional overlap**
  with the model age groups, so narrow bands (e.g. PMC dose ages) are no longer
  silently dropped and coarse groups are no longer over-treated. Co-deployed
  chemoprevention types share a coverage-weighted prophylaxis duration.
* **Antimalarial resistance** (early-treatment-failure) now applies to a drug used
  for **chemoprevention** (SMC/MDA/PMC), not only clinical treatment.
* **IRS** now **accumulates across spray rounds** (sprayed protection never
  expires in the IBM), matching the bed-net accumulation.
* **TBV** targets the exact `ages` set (non-contiguous ages honoured).
* `set_equilibrium`'s custom **`eq_params`** are read explicitly.
* New **all-infection incidence** output (`n_inc_*`; `incidence_rendering_*`).
* Remaining mean-field approximations (time-varying multi-drug first-line switch;
  seasonal PEV boosters) are **warned**, not silently applied.

## Interface

* `run_simulation_ode()` accepts an unmodified `malariasimulation` parameter
  list and returns a wide, `malariasimulation`-style daily count table.
* `get_epi_outputs()` returns `postie`-format rates and prevalence.
* `default_age_lower()` provides the default graded age grid.

## Validation

* Compared against `malariasimulation` across PfPR–EIR, age-prevalence profile,
  bed nets, seasonality, SMC and custom demography, plus clinical & severe
  incidence time series (see `comparison/`). LM prevalence and clinical
  incidence agree with the IBM to within ~1%.

## Scope

* *P. falciparum* only; *P. vivax* and the individual-mosquito code path are
  intentionally not supported. Severe incidence (and derived DALYs) are
  indicative (can differ from the IBM by ~5–20%, largest at low EIR).
