# malariaode 0.0.0.9000

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
