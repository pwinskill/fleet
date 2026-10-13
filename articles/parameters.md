# Parameter reference

`fleet` is not yet validated for real use; see the [status
note](https://pwinskill.github.io/fleet/#not-ready-for-real-use).

This article records, argument by argument, how `fleet` handles the
user-facing parameterisation functions of the `malariasimulation`
individual-based model (IBM), v3.0.0: whether each argument is
reproduced exactly, what is lost if not, and what to do about it. The
`parameters` argument, the list being modified, is omitted throughout.
Section references (§) are to
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md);
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
covers what to do differently because `fleet` is a mean field.

## Status labels

| Label | Meaning | Action |
|----|----|----|
| Exact | Modelled exactly, or algebraically equivalently for the mean field | None. Differences from the IBM fall within its Monte Carlo noise |
| Approximate | Modelled, but the mean field loses something; the row says what | Suitable for comparing scenarios that share the approximation. Check against the IBM if the reported result is dominated by this mechanism |
| Inert | Accepted without error but has no effect on the run, by design rather than by omission | None. Do not tune it |
| Ignored | Not modelled; the run continues, in some cases with a warning | Decide whether the run is valid without it |
| Unsupported | Errors at input, or never reaches the model | Remove it from the parameter list |

Where a row carries two labels, the second is a caveat on the first.

## Index

| Function | Status |
|----|----|
| [`get_parameters()`](#get_parameters) | Exact |
| [`set_species()`](#set_species) | Exact |
| [`set_equilibrium()`](#set_equilibrium) | Exact |
| [`set_demography()`](#set_demography) | Approximate |
| [`set_drugs()`](#set_drugs) | Approximate |
| [`set_clinical_treatment()`](#set_clinical_treatment) | Approximate |
| [`set_antimalarial_resistance()`](#set_antimalarial_resistance) | Approximate |
| [`set_bednets()`](#set_bednets) | Approximate |
| [`set_spraying()`](#set_spraying) | Approximate |
| [`set_carrying_capacity()`](#set_carrying_capacity) | Approximate |
| [`set_mda()`, `set_smc()`](#set_mda) | Approximate |
| [`set_pmc()`](#set_pmc) | Approximate |
| [`create_pev_profile()`](#create_pev_profile) | Exact |
| [`set_pev_epi()`](#set_pev_epi) | Approximate |
| [`set_mass_pev()`](#set_mass_pev) | Approximate |
| [`set_tbv()`](#set_tbv) | Approximate |
| [`set_epi_outputs()`](#set_epi_outputs) | Exact |
| [`set_parameter_draw()`](#set_parameter_draw) | Exact |
| [Helpers and run functions](#helpers) | Mixed |

## `get_parameters()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `overrides` | Named list overriding any default parameter | *Exact.* Passed through whole; see the breakdown by group below |
| `parasite` | `"falciparum"` or `"vivax"` | *Exact.* Both. `"vivax"` runs the *P. vivax* block (§V). MDA, SMC and PMC are refused under vivax, as `malariasimulation` cannot run them |

### `overrides` by group

| Group | Parameters | `fleet` |
|----|----|----|
| Initial state proportions | `s_proportion`, `d_proportion`, `a_proportion`, `u_proportion`, `t_proportion` | *Unsupported.* Not read. `fleet` seeds at the `malariaEquilibrium` fixed point for `init_EIR` (or, under a custom demography, for the EIR the IBM’s mosquito sizing supports; see [`set_equilibrium()`](#set_equilibrium)), not at user-supplied proportions |
| Initial immunity | `init_ib`, `init_ica`, `init_iva`, `init_icm`, `init_ivm`, `init_id`, and for vivax `init_iaa`, `init_iam`, `init_hyp` | *Unsupported.* Not read, for the same reason (vivax seeds at the `malariaEquilibriumVivax` solution) |
| Population size | `human_population` | *Exact.* Scales output counts only. The model is per capita, so run time does not depend on it |
|  | `human_population_timesteps` | *Unsupported.* Not read. Population is conserved; there is no dynamic population size |
| Baseline demography | `average_age` | *Exact.* Sets $`\eta = 1/\text{average_age}`$, the constant death hazard and the equilibrium age structure when `custom_demography = FALSE` |
|  | `custom_demography` | *Exact.* Switches to the [`set_demography()`](#set_demography) path |
| Biting heterogeneity | `a0`, `rho` | *Exact.* $`\psi_i = 1 - \rho e^{-a_i/a_0}`$ |
|  | `sigma_squared` | *Exact.* Log-normal variance of $`\zeta`$ |
|  | `n_heterogeneity_groups` | *Exact.* Number of Gauss–Hermite nodes $`n_z`$. Must be a whole number $`\geq 3`$: fewer nodes are not a quadrature on the log-normal at all and only move every human off it |
|  | `enable_heterogeneity` | *Exact.* `FALSE` collapses to a single node with $`\zeta = 1`$. Under vivax the seed is the IBM’s for this case: everyone from the heterogeneous equilibrium’s first group, among mosquitoes sized from the whole of it (§V.6) |
| Aquatic mosquito | `del`, `dl`, `dpl`, `me`, `ml`, `mup`, `gamma` | *Exact.* Used as given in the $`E`$/$`L`$/$`P_L`$ equations |
| Adult mosquito | `mum`, `beta`, `total_M`, `blood_meal_rates`, `Q0`, `foraging_time`, `species`, `species_proportions` | *Exact.* Per species. `total_M` is re-derived: from `init_EIR` under the default demography, and as the IBM’s own value under a custom one ([`set_equilibrium()`](#set_equilibrium)) |
|  | `init_foim` | *Unsupported.* Not read. `fleet` recomputes FOIM from its own seeded human infectivity, which must be consistent with its seed |
| Seasonality | `model_seasonality`, `g0`, `g`, `h`, `rainfall_floor` | *Exact.* Truncated Fourier rainfall drives $`K_s(t) = K_{0,s}\,\text{scaler}(t)\,\text{rainfall}(t)/\bar{R}`$ on a daily grid |
| Flexible carrying capacity | `carrying_capacity`, `carrying_capacity_timesteps` | *Exact.* See [`set_carrying_capacity()`](#set_carrying_capacity). `carrying_capacity_scalers` is not a `get_parameters()` default and cannot be passed in `overrides`; only `set_carrying_capacity()` creates it |
|  | `carrying_capacity_values` | *Unsupported.* Not read, by `fleet` or by `malariasimulation`; only the `_scalers` form is used |
| Human disease and immunity constants | `dd`, `dt`, `da`, `du` | *Exact.* Back-translated to $`r_D, r_T, r_A, r_U`$, then converted to the IBM’s per-day exit probability $`1-e^{-1/d}`$ so the realised dwell matches |
|  | `rb`, `rc`, `rva`, `rid`, `rm`, `rvm` | *Exact.* Immunity decay time constants $`d_{I_B}, d_{I_{CA}}, d_{I_{VA}}, d_{I_D}`$ and the maternal decay $`d_m, d_{vm}`$ |
|  | `ub`, `uc`, `uv`, `ud` | *Exact.* Refractory periods, applied as $`u_{\text{eff}} = \lceil u\rceil - 1`$ (§B.5) |
|  | `pcm`, `pvm` | *Exact.* $`P_M`$, $`P_{VM}`$ maternal transfer fractions |
|  | `b0`, `b1`, `ib0`, `kb` | *Exact.* $`b`$ Hill function |
|  | `phi0`, `phi1`, `ic0`, `kc` | *Exact.* $`\phi`$ Hill function |
|  | `theta0`, `theta1`, `iv0`, `kv`, `fv0`, `av`, `gammav` | *Exact.* $`\theta`$ Hill function and its age modifier $`f_v`$ |
|  | `d1`, `id0`, `kd`, `fd0`, `ad`, `gammad` | *Exact.* $`q`$ Hill function and its age modifier $`f_d`$ |
|  | `cd`, `cu`, `ct`, `gamma1` | *Exact.* Infectivity by state. `ct` is superseded by the drug-linked $`c_T(t)`$ whenever clinical drugs are set |
|  | `de`, `delay_gam`, `dem` | *Exact.* $`\tau_E`$, $`\tau_l`$, $`\tau_M`$: the three delay lines, read as the IBM reads them |
| Vector-control shape | `phi_bednets`, `phi_indoors`, `k0` | *Exact.* Per species, in the $`a(t)`$/$`\mu(t)`$ algebra |
| Output bands | `age_group_rendering_*`, `incidence_rendering_*`, `clinical_incidence_rendering_*`, `severe_incidence_rendering_*`, `prevalence_rendering_*` | *Exact.* Each drives its own output family, as in `malariasimulation` (§G) |
|  | `ib_`/`id_`/`ica_`/`iva_`/`icm_`/`ivm_rendering_*`, and for vivax `iaa_`/`iam_`/`hypnozoites_rendering_*` | *Exact.* Mean immunity by band, `<name>_mean_*`, with the band’s population as `n_*`, as `malariasimulation` renders it. `idm_rendering_*` names no immunity the IBM has, and renders nothing in either model |
| Metapopulation test-and-treat | `rdt_intercept`, `rdt_coeff` | *Unsupported.* Not read. They parameterise the PCR-to-RDT conversion used by `run_metapop_simulation()`’s mixing, which `fleet` does not support ([Helpers and run functions](#helpers)) |
| Engine and solver | `mosquito_limit`, `individual_mosquitoes`, `r_tol`, `a_tol`, `ode_max_steps`, `progress_bar` | *Unsupported.* Not read. `fleet` is always compartmental and advances one day at a time, with no ODE solver to control |
| Vivax only | `b`, `f`, `gammal`, `kmax`, `philm_*`, `alm50`, `klm`, `dpcr_*`, `apcr50`, `kpcr`, `ca`, `ra`, `ua` | *Exact* under `parasite = "vivax"` (§V). *Approximate:* for `kmax`, the batch counts from 4 to `kmax` - 1 share levels in pairs, with a Poisson distribution inside each (§V.1); for `ua`, IAA’s 44-day refractory window runs as four 11-day stages (§V.3). The shared names (`da`, `dd`, `dt`, `cd`, `cu`, `ct`, `rc`, `uc`, `rm`, `pcm`, `phi0`, `phi1`, `ic0`, `kc`) carry the vivax values on a vivax list |
|  | `n_with_hypnozoites_rendering_*`, and `incidence_relapse_rendering_*` set on the list | *Exact.* `n_with_hypnozoites_*`, with the band’s population as `n_*`, and `n_inc_relapse_*` |

### Fields specific to `fleet`

These four fields are not `malariasimulation` parameters. Each may be
added to the list, and each defaults to the validated choice.

| Field | Default | Effect |
|----|----|----|
| `bite_dedup` | `1` | *P. falciparum* only; vivax counts every bite, as the IBM does. `1` reproduces the IBM’s per-timestep bite deduplication (saturating hazard, §B.2). `0` uses the linear $`b\varepsilon`$ that `malariaEquilibrium` assumes, which reduces the drift of prevalence off the seed from about 0.4% to about 0.1% at EIR 20. It does not remove it, since it changes only the hazard and not the boosting or sojourn forms. Values other than 0 or 1 (`TRUE`/`FALSE` accepted) are an error |
| `acquired_immunity_offset` | `0`, or `0.5` for vivax | $`\delta`$ in the $`b`$/$`\phi`$/$`\theta`$ Hill functions. The IBM adds 0.5 to each person’s positive acquired immunity. Falciparum reads its curves at a stratum mean, where `0` is the better match against the IBM ensemble mean; the vivax curves are read at quadrature nodes that stand for people, and take the IBM’s `0.5` (and `0` with `immunity_spread = FALSE`, back at the mean). Must lie in $`[0, 1]`$ |
| `immunity_spread` | `TRUE` | *P. vivax* only. Carries the second moment of immunity in each age, heterogeneity and batch cell, and averages the immunity curves over a gamma distribution with that mean and variance (§V.3); `FALSE` reads them at the cell mean, with no offset. Must be `TRUE` or `FALSE` |
| `hold_init_EIR` | `FALSE` | Has an effect only under [`set_demography()`](#set_demography). By default `fleet` matches the IBM’s mosquito density rather than its nominal EIR, so the same list realises the same transmission in both models; `init_EIR` is then the EIR the IBM was sized at, not the one `fleet` realises. `TRUE` seeds at `init_EIR` directly. Must be `TRUE` or `FALSE` |

## `set_species()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `species` | List of species parameter lists (`gamb_params`, `arab_params`, `fun_params`, `kol_params`, `steph_params`, or custom) | *Exact.* Each entry’s `blood_meal_rates`, `foraging_time`, `Q0`, `phi_bednets`, `phi_indoors` and `mum` become per-species vectors driving an independent mosquito sub-model, coupled to the others only through the shared human population |
| `proportions` | Relative abundance, summing to 1 | *Exact.* Splits `total_M` and hence each species’ baseline carrying capacity. A species at proportion 0 stays inert; its $`K_0`$ is floored at a negligible positive value so the larval term does not evaluate $`0/0`$ |

## `set_equilibrium()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `init_EIR` | Target EIR to seed from | *Exact.* Stored as `parameters$init_EIR`, which is where [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md) reads it; the function has no `init_EIR` argument of its own. It means the same in both models: under the default demography it is the EIR `fleet` realises; under `set_demography()` it sizes the mosquito population exactly as the IBM’s does, and `fleet` seeds at the EIR that density supports ([`set_demography()`](#set_demography)) |
| `eq_params` | Custom `malariaEquilibrium` parameter set | *Exact.* Merged over `fleet`’s own back-translation, taking precedence. `set_equilibrium()` writes this field even when passed `NULL` (storing its own back-translation), so after any call to it that stored set, not `R/translate_params.R`, supplies every shared constant (§F) |
| `EIR_population_input` | `"adult"` (default) or `"total"` | *Exact.* Handled upstream: `malariasimulation` converts a total-population EIR to adult EIR before storing it, and `fleet` always interprets `init_EIR` as adult EIR (bites per adult per year) |

`set_equilibrium()` also calls `parameterise_mosquito_equilibrium()`,
which sizes the adult mosquito population `total_M` from the human
equilibrium under the IBM’s default exponential age structure
(`equilibrium_total_M()`; custom mortality never enters). `fleet` never
reads the stored `parameters$total_M`. It re-derives the quantity
itself, in one of two ways:

- **Default demography.** `total_M` is chosen so that
  $`\sum_s a_s I_{M,s} = \text{init_EIR}/365`$ holds exactly under the
  human infectivity of `fleet`’s own equilibrium at the treatment
  coverage in force at timestep 1. This is the IBM’s formula evaluated
  on `fleet`’s grid, so the two agree up to discretisation and
  `init_EIR` is the EIR the mosquitoes deliver. The humans themselves
  are seeded under the coverage at timestep 0, as the IBM seeds them
  (§F).
- **Custom demography** (`set_demography()`). Holding `init_EIR` here
  would give the two models different mosquito populations, and so
  different transmission, from the same parameter list: under the custom
  age structure the IBM’s mosquito density supports less (or more)
  transmission, and the IBM drifts there. `fleet` therefore reproduces
  the IBM’s `total_M` exactly (`ibm_total_M()`: the same
  `human_equilibrium()` call on the IBM’s 0–99.9 y grid and the same
  `equilibrium_total_M()` arithmetic), root-finds the EIR at which its
  own equilibrium under the custom age structure has that density, and
  seeds there. `parameters$hold_init_EIR = TRUE` seeds at `init_EIR`
  instead. Either way the seed is the `malariaEquilibrium` solution at
  that EIR, so it relaxes slightly, as every seed does (§F).

Calling `parameterise_total_M()` or
`parameterise_mosquito_equilibrium()` directly does no harm and has no
effect.

## `set_demography()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `agegroups` | Upper edges of the death-rate age groups, in days | *Exact.* Model age groups are binned into them by midpoint using right-closed intervals, matching `.bincode(age, c(0, agegroups))` |
| `timesteps` | When each death-rate row takes effect | *Exact.* Knots of a constant-interpolated $`\mu_i(t)`$, so custom demography is time-varying: a demographic transition is modelled, not frozen. Each row applies from its own timestep, as the IBM reads it (§H) |
| `deathrates` | `[length(timesteps) × length(agegroups)]` daily death rates | *Exact.* Used directly as $`\mu_i(t)`$. The $`t = 0`$ row also sets the equilibrium age structure onto which the human seed is rescaled. The open-ended top model group is looked up one day above its lower edge, so it takes the rate of the band it belongs to rather than the one below. *Approximate:* model ages above the top `agegroups` edge take the top rate, whereas the IBM removes them. A warning is given; raise `default_age_lower(max_age =)` or extend `agegroups` to match |

With a custom demography the mosquito population is sized as the IBM
sizes it (see [`set_equilibrium()`](#set_equilibrium)): `fleet` takes
the IBM’s `total_M` and seeds at the EIR that its own equilibrium under
this age structure supports, which is generally not `init_EIR`. An older
population sustains less transmission from the same mosquitoes: the
comparison article’s scenario seeds at EIR 13.9 for `init_EIR = 20`,
matching the IBM’s realised 13.8.

## `set_drugs()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `drugs` | List of 4-element vectors `c(efficacy, rel_c, prophylaxis_shape, prophylaxis_scale)` (for example `AL_params`, `DHA_PQP_params`, `SP_AQ_params`) | *Approximate.* Each element is handled as below |
|   `drug_efficacy` | Probability treatment clears the infection | *Exact.* Multiplies coverage into $`f_t^{\text{eff}}`$, and multiplies chemoprevention pulse fractions |
|   `drug_rel_c` | Infectivity of a treated case relative to `cd` | *Exact.* $`c_T = c_d \times \text{rel_c}`$, coverage-share-weighted across active drugs and time-varying |
|   `drug_prophylaxis_shape`, `drug_prophylaxis_scale` | Weibull prophylaxis survival | *Exact* in both moments, on the protected days a dose gives, $`W(1), W(2), \dots`$ The post-treatment chain carries the protected days left once out of $`T`$, $`\sum_n (1-(1-r_T)^{n-1})W(n)`$ (§B.4), and its length $`k_P`$ matches the variance of the whole $`T + P`$ sojourn: 9 stages for SP-AQ, 10 for DHA-PQP, 1 for AL. The chemoprevention chain, with no $`T`$ stage in front, matches the dose’s own mean and variance: 10 stages for SP-AQ, 9 for DHA-PQP. No chain has more stages than its mean has days, and none more than 20. `ode_tuning(n_ph =)` overrides the count. *Approximate:* a stage chain is not a Weibull, and it protects fully while it lasts, where the IBM’s protection is partial and fades. The $`T`$ stage adds variability the IBM’s parallel clock does not have, so for AL the protection curve is reproduced in its integral, not its shape |
|   `drug_hypnozoite_efficacy`, `drug_hypnozoite_prophylaxis_shape`, `drug_hypnozoite_prophylaxis_scale` | *P. vivax* radical cure | *Exact.* The treated lose every batch with the hypnozoite efficacy, drawn independently of blood-stage clearance, and start liver-stage protection: a chain of levels matched to the protected days, $`M^2/(V+M)`$ stages, eleven for tafenoquine (§V.5). *Approximate:* a chain protects fully while it lasts; there is none where the protection ends inside the blood-stage protection given with it, as primaquine’s does inside chloroquine’s |

## `set_clinical_treatment()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `drug` | Index into the drug table | *Exact.* Identifies which drug’s efficacy, `rel_c` and prophylaxis enter the blend |
| `timesteps` | When coverage changes | *Exact.* Knots of the step-interpolated $`f_t(t)`$ and of the drug-mix series |
| `coverages` | Fraction of clinical cases treated with this drug | *Exact.* $`f_t(t) = \min\left(1, \sum_{\text{drugs}} \text{cov}_d(t)\right)`$. `set_clinical_treatment()` errors if the summed coverage exceeds 1 at any timestep, so the cap at 1 is defensive and never binds on a valid parameter set. *Approximate:* the drug-linked quantities $`\text{eff}(t)`$, $`c_T(t)`$, $`r_P(t)`$ are blended by each drug’s instantaneous coverage share, so a first-line switch is modelled, but a mixed first line is represented by its blend rather than as parallel sub-populations |

## `set_antimalarial_resistance()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `drug` | Which clinical drug carries resistance | *Exact.* Matched to the clinical-treatment drug index, and to the MDA, SMC and PMC drug, so a resistant chemoprevention drug clears fewer infections |
| `timesteps` | When resistance levels change | *Exact.* Knots of $`\text{ETF}(t)`$ and $`\text{SPC}(t)`$, merged with the treatment-coverage change times |
| `artemisinin_resistance_proportion` | Fraction of infections that are artemisinin-resistant | *Exact.* Multiplies both the ETF and SPC probabilities |
| `partner_drug_resistance_proportion` | Partner-drug resistance | *Unsupported.* `malariasimulation` requires 0, so it cannot reach the model |
| `slow_parasite_clearance_probability` | Probability a treated resistant case clears slowly | *Exact.* $`\text{SPC}(t)`$, the split of the treated inflow between $`T`$ and $`T_s`$, coverage-share-weighted across drugs |
| `early_treatment_failure_probability` | Probability treatment fails early | *Exact.* $`\text{ETF}(t)`$ reduces $`f_t^{\text{eff}}`$, diverting those cases to $`D`$ |
| `late_clinical_failure_probability` | Late clinical failure | *Unsupported.* `malariasimulation` requires 0 |
| `late_parasitological_failure_probability` | Late parasitological failure | *Unsupported.* `malariasimulation` requires 0 |
| `reinfection_during_prophylaxis_probability` | Reinfection while prophylactic | *Unsupported.* `malariasimulation` requires 0 |
| `slow_parasite_clearance_time` | Mean duration of slow clearance | *Approximate.* $`r_T^{\text{slow}} = 1 - e^{-1/\text{dt_slow}}`$, a single scalar with the same whole-day conversion as the other sojourns (§B.4). With one resistant drug the blend is exact; with several it is the coverage-weighted blend evaluated at peak resistance (peak rather than final, so a rise-then-fall schedule is not collapsed to zero) |

## `set_bednets()`

Status: approximate. Net efficacy is averaged over the net-using
fraction into per-species $`a_s(t)`$ and $`\mu_s(t)`$, so protection is
not correlated within individuals across bites. Within that, every
argument is handled exactly.

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `timesteps` | Distribution dates | *Exact.* Every past distribution contributes: at each grid time the population is the mixture over all rounds, each weighted by its most-recent-receipt probability $`\text{cov}_d\prod_{j>d}(1-\text{cov}_j)`$ |
| `coverages` | Fraction receiving a net at each distribution | *Exact.* As above; repeated distributions accumulate rather than only the latest applying |
| `dn0` | `[n_timesteps × n_species]` probability a net kills a mosquito | *Exact.* Decayed by net age: $`d_{n}(s_n) = d_{n0}e^{-s_n/\gamma_n}`$ |
| `rn` | `[n_timesteps × n_species]` initial repelling probability | *Exact.* Decayed towards `rnm`: $`r_n(s_n) = (r_n - r_{nm})e^{-s_n/\gamma_n} + r_{nm}`$ |
| `rnm` | `[n_timesteps × n_species]` minimum (asymptotic) repelling probability | *Exact.* As above |
| `gamman` | Insecticide decay time constant per distribution | *Exact.* The $`\gamma_n`$ above; each distribution decays on its own clock |
| `retention` | Mean net-retention time (log-uniform model) | *Exact.* Exponential usage survival $`e^{-s_n/\text{retention}}`$ |
| `logistic_half_life`, `logistic_k` | Alternative logistic retention | *Exact.* The survival of the IBM’s logistic retention time: $`S(s_n) = \exp\!\left(-k\frac{r^2}{1-r^2}\right)`$ for $`r = s_n/l < 1`$ (else 0), with $`l = \text{half_life}/\sqrt{1 - k/(k-\log 0.5)}`$, which gives $`S(\text{half_life}) = 0.5`$ exactly |

## `set_spraying()`

Status: approximate. As for nets, protection is population-averaged, and
exact within that. Sprayed protection never expires in the IBM, so the
mixture over past rounds has no retention factor.

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `timesteps` | Spray dates | *Exact.* Mixture over all past rounds, weighted by most-recent-spray probability |
| `coverages` | Fraction of houses sprayed | *Exact.* As above |
| `ls_theta`, `ls_gamma` | `[n_timesteps × n_species]` mortality logistic parameters | *Exact.* $`l_s(s_s) = \text{logit}^{-1}(\theta + \gamma s_s)`$ |
| `ks_theta`, `ks_gamma` | Feeding-success logistic parameters | *Exact.* $`k_s = k_0\,\text{logit}^{-1}(\theta + \gamma s_s)`$ |
| `ms_theta`, `ms_gamma` | Deterrence logistic parameters | *Exact.* $`m_s = \text{logit}^{-1}(\theta + \gamma s_s)`$ |

The repellency $`r_s`$ and survival $`s_s`$ derived from these are the
same spray outcome per individual, so `fleet` forms the joint mean
$`\overline{(1-r_s)s_s}`$ rather than multiplying two independent
averages. The IBM recomputes these every timestep; `fleet` reads
$`a_s(t)`$ and $`\mu_s(t)`$ once a day from series that follow the
within-round logistic decay.

## `set_carrying_capacity()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `timesteps` | When each scaler row takes effect | *Exact.* Step change in $`K_s(t)`$, composed multiplicatively with seasonality. $`K_s`$ is read at the start of each day, as the IBM reads it, so a change scheduled for timestep $`t`$ applies from day $`t + 1`$, as in the IBM (§D, §H) |
| `carrying_capacity_scalers` | `[n_timesteps × n_species]` multipliers on baseline $`K_0`$ | *Exact.* Multiplies baseline $`K_0`$ per species. *Approximate:* a scaler of exactly 0 is floored at $`K_0 \times 10^{-4}`$ to keep the larval term finite, so complete vector elimination is approached but not reproduced exactly |

## `set_mda()` and `set_smc()`

`set_smc()` takes the same arguments as `set_mda()` and is handled
identically.

Status: approximate. A mass campaign is applied as a pulse rather than
as per-individual events, and its protection is all or nothing while it
lasts. As in the IBM, the clinical and LM-detectable asymptomatic cases
a round treats pass through the treated state first, detectable and
infectious for the Tr stay.

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `drug` | Drug administered | *Exact.* Its `drug_efficacy` scales the treated fraction. Its Weibull protection sets the chains’ means and stage counts: $`k = M^2/(V+M)`$ for the protection from the dose, and the post-treatment rule for the protection left after the Tr stay (§B.4). Its `drug_rel_c` scales the infectivity of the treated. Antimalarial-resistance ETF on this drug reduces the treated fraction, and SPC sends its share of the detectable to the slow-clearing treated state |
| `timesteps` | Round dates | *Exact.* The run pauses after each round’s day; the pulse takes effect the day after the scheduled timestep |
| `coverages` | Fraction of the target band reached | *Exact.* Treated fraction $`= \text{coverage}\times\text{efficacy}\times(1-\text{ETF})\times o_i`$. A zero-coverage round has no effect |
| `min_ages`, `max_ages` | Per-round vectors (one entry per `timesteps` entry) giving the target band in days, inclusive at both ends | *Exact.* Per-round bands are honoured, each mapped onto the model age grid by fractional overlap $`o_i`$, so narrow bands are not dropped and coarse groups are not over-treated. The absorbing top group is treated as fully covered if and only if the band reaches its lower edge. *Approximate:* `fleet` treats the band as half-open, $`[\ell, u)`$, one day narrower than the IBM’s; this is negligible except for very narrow (PMC-style) bands |

Co-deployed chemoprevention types (SMC, MDA, PMC) share the two
chemoprevention chains, moment-matched to the mixture of their drugs’
protection. Each drug is weighted by its type’s total scheduled coverage
(the sum over all its rounds), so a type with more rounds carries
proportionally more weight, whatever its per-round coverage.

## `set_pmc()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `drug` | Drug administered | *Exact.* As for MDA and SMC |
| `timesteps` | Coverage-change schedule (not dose dates; PMC delivery is age-triggered) | *Exact.* Read as a coverage schedule |
| `coverages` | Fraction of infants receiving each dose | *Exact.* The coverage in force at each pulse time is used |
| `ages` | Ages (days) at which doses are given | *Approximate.* The IBM gives a dose when an individual reaches each dose age. `fleet` approximates this as pulses every 30 days over a 30-day band starting at each dose age, with per-group overlap weighting, keeping the effective dose at about coverage × efficacy |

## `create_pev_profile()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `vmax` | Maximum efficacy | *Exact.* |
| `alpha`, `beta` | Hill shape and scale mapping antibody titre to efficacy | *Exact.* |
| `cs` | `c(mu, sigma)` of $`\log`$ peak antibody titre | *Exact*, including the sigma |
| `rho` | `c(mu, sigma)` of the logit short-lived antibody fraction | *Exact*, including the sigma |
| `ds`, `dl` | `c(mu, sigma)` of $`\log`$ short and long antibody half-lives | *Exact*, including the sigmas |

`malariasimulation` re-draws all four antibody parameters from their
profile distributions at every timestep, for every vaccinated
individual; they are not stored per person. Population efficacy is
therefore exactly $`\mathbb{E}\left[\text{Hill}(\text{Ab})\right]`$, not
$`\text{Hill}(\text{Ab at the median})`$. The sigmas are large and the
Hill function is nonlinear, so the two differ materially: evaluating at
the median overstated R21 efficacy by up to about 4.7 percentage points.
`fleet` integrates over the 4-D distribution with a tensor Gauss–Hermite
rule (`options(fleet.pev_gq =)`, default 7 nodes per axis), validated
against a Monte Carlo of the IBM’s own sampler to under
$`5\times10^{-4}`$. Because the IBM’s draws are independent between days
and between people, there is no within-person correlation for the mean
field to lose, and the quadrature is exact rather than approximate. The
mean-efficacy curve is cached per profile on a daily grid.

## `set_pev_epi()`

Status: approximate. The routine schedule and the full booster sequence
are modelled; seasonal boosters are approximated, with a warning.

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `profile` | Primary-series PEV profile | *Exact.* |
| `coverages` | Coverage of the routine programme over time | *Exact.* Time-varying: each cohort is protected at the coverage in force on its own first-dose date, not at a single coverage |
| `timesteps` | When each coverage takes effect | *Exact.* Also gates eligibility: a cohort is protected only once its first dose falls inside the programme. Gating follows `malariasimulation`, which keys on `pev_epi_coverages`/`pev_epi_timesteps` and never reads `parameters$pev` |
| `age` | Age (days) at the first dose | *Exact.* Efficacy starts at `age + max(pev_doses)`, after the final primary dose |
| `min_wait` | Minimum time since the last vaccination | *Inert.* For a routine EPI programme this guards re-vaccination and seasonal-booster timing, not the primary-plus-booster schedule `fleet` models, so it has no effect on the result. For mass campaigns it does matter; see [`set_mass_pev()`](#set_mass_pev) |
| `booster_spacing` | Days from the final primary dose to each booster | *Exact.* The full booster sequence is modelled: the vaccinated are partitioned by the most recent booster each has reached, and each stratum carries its own decayed efficacy |
| `booster_coverage` | Matrix of conditional booster coverages | *Exact.* Each booster’s coverage is read at its own administration date (matching `coverage[match_timestep(timesteps, admin_time), booster]`), which reduces to the single-row lookup in the usual case |
| `booster_profile` | Per-booster profiles | *Exact.* Each booster uses its own profile, integrated over the antibody distribution as in [`create_pev_profile()`](#create_pev_profile) |
| `seasonal_boosters` | Times the first booster relative to the start of the calendar year: `booster_spacing[1]` becomes a day of the year. It makes no reference to the fitted seasonality; the user chooses the day, for example with `peak_season_offset()` | *Approximate.* Treated as a fixed days-since-primary schedule. `fleet` warns when this is set |

## `set_mass_pev()`

Status: approximate. Repeated campaigns are the weak point: `min_wait`
re-vaccination exclusion is not applied, and each campaign reaches
people independently of the last. A campaign’s protection ages with the
cohort it vaccinated, and campaigns combine as the most recent dose, as
in the IBM. A single campaign is unaffected.

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `profile` | Primary-series profile | *Exact.* |
| `timesteps` | Campaign dates | *Exact.* The IBM keeps only the most recent vaccination per person, so the population protection is a most-recent-dose mixture, $`\sum_k c_k w_k \prod_{j>k}(1 - c_j w_j)\,e_k`$ over completed campaigns, with $`w`$ each group’s share of a campaign’s cohort. This is the structure `fleet` uses for [bed nets](#set_bednets). *Approximate:* reaching someone in one campaign does not make them more likely to be reached in the next |
| `coverages` | Coverage per campaign (recycled if length 1) | *Exact.* |
| `min_ages`, `max_ages` | Target age bands | *Exact.* Every band is applied at every campaign, each mapped onto the model age grid by fractional overlap as for [MDA and SMC](#set_mda), so a band narrower than the groups it falls between still vaccinates its share rather than nobody. The vaccinated then age: at each later day the band has moved up by the time since the campaign, so protection stays with the cohort that received it and nobody born since shares it. The bands of one campaign are one campaign: their weights add into a single covered fraction per group, capped at 1. A band overlapping no group gives a warning |
| `min_wait` | Minimum time since the last vaccination | *Ignored*, and it matters here. In `malariasimulation` each campaign excludes anyone vaccinated within `min_wait` of it, so with repeated campaigns this is the main control on who is re-vaccinated. `fleet` treats every campaign as reaching its whole target band independently, so a `min_wait` longer than the campaign spacing leads to over-vaccination |
| `booster_spacing` | Days from the final primary dose to each booster | *Exact.* Full sequence, as in [`set_pev_epi()`](#set_pev_epi) |
| `booster_coverage` | Conditional booster coverage matrix | *Exact.* Read at each booster’s administration date |
| `booster_profile` | Per-booster profiles | *Exact.* |

## `set_tbv()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `timesteps` | Vaccination dates | *Exact.* Antibody titre decays from each round’s date. The IBM overwrites each person’s vaccination date, and rounds combine as the most recent dose, as for mass PEV. *Approximate:* each round reaches people independently of the last |
| `coverages` | Fraction vaccinated | *Exact.* |
| `ages` | Integer ages (years) vaccinated | *Exact.* Each year $`y`$ is the day interval $`[365y, 365(y+1))`$, matched to the age groups by fractional overlap, and the vaccinated then age, the year set moving up with them. Overlap has to be fractional because, between its pinned ages, the default grid puts edges wherever its group budget falls (above 15, every 0.38 years), so a whole-year target need not line up with any group edge. A set overlapping no group gives a warning |

Transmission-reducing activity is mapped to state-specific
transmission-blocking activity by the IBM’s own transform
$`\text{TBA}(m_x, k, \text{TRA})`$, with per-state $`m_x \in
\{\texttt{tbv_mu}, \texttt{tbv_ma}, \texttt{tbv_md}, \texttt{tbv_mt}\}`$,
giving four independent per-age multipliers on $`U`$, $`A`$, $`D`$ and
$`T`$ infectivity. The `tbv_*` shape constants (`tbv_tau`, `tbv_rho`,
`tbv_ds`, `tbv_dl`, `tbv_tra_mu`, `tbv_gamma1`, `tbv_gamma2`, `tbv_k`)
are used as given.

## `set_epi_outputs()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `age_group` | Bands for population counts | *Exact.* `n_age_*`, also emitted over the union of all other families’ bands, sorted, as `malariasimulation` does |
| `incidence` | Bands for all-infection incidence | *Exact.* `n_inc_*`, and `p_inc_*` (0, as `malariasimulation` 3.0.0 renders it) |
| `clinical_incidence` | Bands for clinical incidence | *Exact.* `n_inc_clinical_*`, `p_inc_clinical_*` |
| `severe_incidence` | Bands for severe incidence | *Exact.* `n_inc_severe_*`, `p_inc_severe_*` |
| `prevalence` | Bands for prevalence | *Exact.* `n_detect_lm_*`, `p_detect_lm_*` (the expected number detected, not a prevalence), `n_detect_pcr_*` |
| `ica`, `icm`, `iva`, `ivm`, `id`, `ib` | Bands for mean immunity | *Exact.* `<name>_mean_*`, with the band’s population as `n_*` |
| `n_with_hypnozoites`, `hypnozoites`, `iaa`, `iam` | *P. vivax* only | *Exact.* `n_with_hypnozoites_*` and `<name>_mean_*`, each with the band’s population as `n_*` |

A family whose band list is empty renders nothing, as in the IBM. A bare
`get_parameters()` leaves four of the five families empty and sets the
prevalence bands to 2–10 y, so it renders the 2–10 prevalence columns
and no incidence. A band `[lower, upper]` holds the whole-day ages
`lower` to `upper`, both ends included, and its tag is formatted as the
IBM formats it (`n_age_182.5_1825.5`, `n_inc_clinical_1e.05_2e.05`).
Every `n_*`/`p_*` pair is the IBM’s sampled count and its expectation; a
mean field’s expected count is both.

## `set_parameter_draw()`

| Argument | malariasimulation meaning | `fleet` |
|----|----|----|
| `draw` | Index 1–1000 into the fitted posterior draws | *Exact.* `malariasimulation` overwrites the immunity and disease constants in the parameter list in place; they then pass through `fleet`’s back-translation like any other parameter. `iv0` is read from the live list: `set_equilibrium()`’s stored copy holds `malariaEquilibrium`’s default whatever the draw, because `malariasimulation`’s back-translation table has that entry the wrong way round, and the IBM’s severe disease reads the live value. Must be called before `set_equilibrium()`, as in the IBM |

## Helpers and run functions

| Function | malariasimulation meaning | `fleet` |
|----|----|----|
| `peak_season_offset(parameters)` | Day of peak seasonal transmission | *Exact.* A pure helper on `g0`, `g` and `h`, which `fleet` reads; useful for timing SMC and booster schedules |
| `get_correlation_parameters(parameters)` | Correlation between intervention recipients | *Unsupported.* Individual-level, with no mean-field analogue |
| `run_simulation(timesteps, parameters, correlations)` | Run the IBM | *Exact.* `run_simulation_ode(timesteps, parameters, correlations, tuning)` mirrors the first three arguments; `tuning` carries the discretisation settings, so nothing epidemiological sits outside the parameter list. A non-`NULL` `correlations` is ignored with a warning: the model assumes independence between interventions |
| `run_simulation_with_repetitions(...)` | Repeat stochastic runs | *Unsupported.* Not applicable: `fleet` is deterministic, so one run is the ensemble mean |
| `run_metapop_simulation(timesteps, parameters, ..., export_mixing, import_mixing)` | Coupled patches | *Unsupported.* `fleet` is a single patch. There is no guard, so a metapopulation setup passed to it is not honoured |
| `run_resumable_simulation(...)` | Resume from saved state | *Unsupported.* Not implemented. `fleet` runs are cheap enough to re-run from the seed |
