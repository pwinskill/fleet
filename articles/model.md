# The blink model: formal specification

This article is the reference description of the system of ordinary
differential equations that `blink` integrates, and of how every
user-facing `malariasimulation` parameterisation function maps onto it.
It is written to be read alongside the model source in
`inst/odin/malaria_ode.R` and the input builders in `R/build_inputs.R` /
`R/interventions.R`.

For a task-oriented introduction with runnable code, see [Getting
started with
`blink`](https://pwinskill.github.io/blink/articles/blink.md). **This
article is the canonical description of the model and of what it
approximates:** §3 says what the state space does and does not carry, §4
covers every `malariasimulation` argument one by one, and §5 lists what
to do differently because this is an ODE and not the IBM. If you read
only one section, read §5.

Everything below describes `malariasimulation` **v3.0.0**.

## 1. Scope and design

`blink` is a deterministic mean-field (ODE) counterpart to the
`malariasimulation` individual-based model (IBM) of *Plasmodium
falciparum*.

*Mean field* means the model tracks the **distribution of the population
over states** rather than the people themselves. Any quantity that
depends on two attributes belonging to the *same* person — owning a net
*and* being vaccinated, say — is therefore unavailable by construction.
Almost everything in §3.2 and §3.4 follows from that one fact.

`blink` targets the IBM’s *code*, not the published equations.
`malariasimulation`’s implementation departs from the analytic
Griffin-style model in four places: it deduplicates bites within a day,
rounds immunity-boosting refractory windows to whole days, holds people
in a state for whole days, and applies extrinsic-incubation survival at
a fixed delay. Wherever the two differ, `blink` reproduces what the IBM
does.

The design constraint that shapes everything below is that **the state
must not carry per-individual attributes**. The IBM tracks, for each
person, their net, their spray status, their vaccine doses and antibody
titres, their drug and the date they took it, and the date each of their
four immunity functions last boosted. A mean-field model cannot afford a
state dimension for each. `blink` therefore keeps only two human
structuring dimensions — age and biting heterogeneity — and represents
everything else either as a **time-varying coefficient** in the ODE, or
as a **discrete state jump** between integration segments. §3.1–3.4 set
out exactly which is which.

## 2. The model at a glance

![Flow diagram of the blink model. The upper half shows the human
compartments — susceptible, untreated clinical disease, two treated
compartments for standard and slow parasite clearance, asymptomatic and
sub-patent infection, and two prophylaxis compartments — with the
infection hazard splitting four ways out of susceptible and recovery
returning everyone to susceptible. A dashed arrow marks the
chemoprevention pulse. The lower half shows the mosquito compartments:
eggs, larvae, pupae, then susceptible, exposed and infectious adults.
Two arrows couple the halves: the entomological inoculation rate driving
human infection, and human infectivity driving the mosquito force of
infection.](model_flow.png)

Flow diagram of the blink model. The upper half shows the human
compartments — susceptible, untreated clinical disease, two treated
compartments for standard and slow parasite clearance, asymptomatic and
sub-patent infection, and two prophylaxis compartments — with the
infection hazard splitting four ways out of susceptible and recovery
returning everyone to susceptible. A dashed arrow marks the
chemoprevention pulse. The lower half shows the mosquito compartments:
eggs, larvae, pupae, then susceptible, exposed and infectious adults.
Two arrows couple the halves: the entomological inoculation rate driving
human infection, and human infectivity driving the mosquito force of
infection.

Every box is a compartment the model integrates; every solid arrow is a
rate in the differential equations. Reading it:

- **Top half, humans.** Infection out of $`S`$ arrives at rate
  $`\Lambda`$ and immediately splits four ways: into untreated clinical
  disease $`D`$, into one of the two treated compartments $`T`$ /
  $`T_s`$, or — if the infection is not clinical — straight into
  asymptomatic infection $`A`$. Recovery runs left to right,
  $`D \to A \to U`$, and everyone eventually returns to $`S`$: from
  $`U`$ directly, from $`T`$ / $`T_s`$ via post-treatment prophylaxis
  $`P`$, and from chemoprevention prophylaxis $`P_c`$.
- **The dashed red arrow** is the one flow that is *not* a rate. MDA,
  SMC and PMC are applied as instantaneous jumps between integration
  segments, clearing a fraction of the target age band into $`P_c`$
  (§E).
- **Bottom half, mosquitoes.** Eggs $`E`$ → larvae $`L`$ → pupae $`P_L`$
  → susceptible adults $`S_M`$, which become exposed $`E_M`$ and then
  infectious $`I_M`$ after surviving the extrinsic incubation period.
- **The two purple arrows are the transmission coupling**, and they are
  what makes this one model rather than two. Infectious mosquitoes drive
  the human infection hazard through the EIR; infectious humans drive
  the mosquito force of infection $`\Lambda^M`$.

What the diagram deliberately does not show is the *structure within
each box*. Every human compartment is really an array over age and
biting heterogeneity, and each of the four immunity functions is another
such array — which is the subject of §3. Nor does it show where
interventions act: bed nets, IRS, vaccines, treatment and seasonality do
not add boxes, they modify the coefficients on the arrows already drawn
(§3.2).

## 3. State space

### 3.1 What has a state dimension

Every human compartment is a two-dimensional array over age group $`i`$
and biting-heterogeneity node $`j`$; every mosquito compartment is
indexed by species $`s`$. Human compartments are **per-capita
densities** (fractions of the total human population), which is why run
time is independent of `human_population`.

| Dimension | Symbol | Default size | Set by | What it indexes |
|----|----|----|----|----|
| Age group | $`i = 1,\dots,n_a`$ | 52 | `age_lower =` / [`default_age_lower()`](https://pwinskill.github.io/blink/reference/default_age_lower.md) | Graded age grid: monthly to 1 y, quarterly to 5 y, yearly to 15 y, 5-yearly to 80 y, plus an absorbing top group. Individuals age by a linear-chain (McKendrick) flow $`r_i`$ between adjacent groups. |
| Biting heterogeneity | $`j = 1,\dots,n_z`$ | 5 | `parameters$n_heterogeneity_groups`, `enable_heterogeneity` | Gauss–Hermite nodes $`\zeta_j`$ and weights $`w_j`$ of the log-normal relative-biting-rate distribution with variance $`\sigma^2`$ (`sigma_squared`); $`\zeta_j = \exp(z_j\sigma - \sigma^2/2)`$. Fixed at birth, never mixed — a person stays in their stratum for life, as in the IBM. |
| Vector species | $`s = 1,\dots,n_v`$ | 1 | `set_species()` | Independent aquatic + adult mosquito populations, each with its own biting rate, mortality, carrying capacity and vector-control response, coupled only through the shared human population. |
| EIR lag stage | $`1,\dots,n_E`$ | 10 (`n_eir`) | `run_simulation_ode(n_eir =)` | Erlang chain approximating the sporozoite-to-infection delay $`\tau_E`$ (`de`). |
| FOIM lag stage | $`1,\dots,n_F`$ | 10 (`n_foim`) | `run_simulation_ode(n_foim =)` | Erlang chain approximating the gametocyte delay $`\tau_l`$ (`delay_gam`). |
| EIP stage | $`1,\dots,n_P`$ | 20 (`n_eip`) | `run_simulation_ode(n_eip =)` | Loss-free Erlang chain approximating the fixed extrinsic incubation period $`\tau_M`$ (`dem`). |

Human compartments, all $`[n_a \times n_z]`$. The last column gives the
column name in the output table, which is where you will meet these in
practice — note that eight compartments are reported in six columns.

| Compartment | odin | Output column | Meaning |
|----|----|----|----|
| $`S`$ | `S` | `S_count` | Susceptible |
| $`D`$ | `D` | `D_count` | Clinical disease, untreated |
| $`T`$ | `Tr` | `Tr_count` (with $`T_s`$) | Successfully treated, standard clearance |
| $`T_s`$ | `Tr_slow` | folded into `Tr_count` | Successfully treated, **slow parasite clearance**. A separate compartment, so the two clearance speeds form a genuine mixture of exponentials rather than one exponential at the blended mean (§B.4) |
| $`A`$ | `A` | `A_count` | Asymptomatic patent infection |
| $`U`$ | `U` | `U_count` | Sub-patent infection |
| $`P`$ | `Ph` | `Ph_count` (with $`P_c`$) | Post-treatment prophylaxis |
| $`P_c`$ | `Ph_c` | folded into `Ph_count` | Chemoprevention prophylaxis (MDA/SMC/PMC), separate because its drug — and hence its decay rate — differs from the clinical first line |
| $`I_B`$ | `IB` | — | Pre-erythrocytic (anti-infection) immunity |
| $`I_{CA}`$ | `ICA` | — | Acquired clinical immunity |
| $`I_D`$ | `ID` | — | Anti-parasite (detection) immunity |
| $`I_{VA}`$ | `IVA` | — | Acquired severe-disease immunity |

The four immunity arrays hold the **mean immunity of the individuals
currently in that (age, heterogeneity) cell**, not an extensive density;
their aging term is correspondingly a difference
$`\rho_i (I_{i-1} - I_i)`$ rather than a flux.

Mosquito compartments, all $`[n_v]`$ except the EIP chain:

| Compartment | Meaning |
|----|----|
| $`E`$, $`L`$, $`P_L`$ | Aquatic early larval, late larval, pupal (odin names `ME`, `ML`, `MP`) |
| $`S_M`$ | Susceptible adult females |
| $`X^{(P)}_{s,1..n_P}`$ | EIP delay chain (odin `Xi`) — a pure delay line, carries no mortality |
| $`E_M`$ | Incubating (exposed) adults (odin `Em_inc`) |
| $`I_M`$ | Infectious adults |

**Total state size.** With the defaults ($`n_a = 52`$, $`n_z = 5`$,
$`n_v = 1`$, $`n_E = n_F = 10`$, $`n_P = 20`$):

``` math
\underbrace{12\,n_a n_z}_{3120\ \text{human}} \;+\; \underbrace{n_E + n_F}_{20\ \text{lags}} \;+\; \underbrace{n_v(6 + n_P)}_{26\ \text{mosquito}} \;=\; 3166 .
```

### 3.2 What is captured *without* a dimension

These are the mechanisms the IBM carries as per-individual attributes
and `blink` carries as **time-varying (and sometimes age-varying)
coefficients**. Each row is a dimension that was *not* added.

| Mechanism | IBM representation | `blink` representation | Cost of the collapse |
|----|----|----|----|
| **Clinical treatment coverage** | Per-person Bernoulli draw at the moment of a clinical episode | Scalar $`f_t(t)`$, step-interpolated over the union of all drugs’ `timesteps`; splits the clinical inflow between $`T/T_s`$ and $`D`$ | None at the flow level |
| **Which drug a person received** | Drug index stored per treated individual, driving their efficacy, infectivity and prophylaxis duration | Three scalar series $`\text{eff}(t)`$, $`c_T(t)`$, $`r_P(t)`$, each a coverage-share-weighted blend across the active clinical drugs, step-interpolated over the coverage change times | A first-line **switch** is modelled (the blend moves in time); a genuinely *mixed* first line is represented by its mean, not its mixture |
| **Post-treatment prophylaxis clock** | Weibull hazard on days-since-treatment, per person | One compartment $`P`$ with exponential exit at $`r_P = 1/(\bar{d}_W - 1/r_T)`$, where $`\bar{d}_W`$ is the Weibull mean | Exact in the equilibrium mean; the transient decay is exponential rather than the sharper Weibull |
| **Early treatment failure** | Per-person draw diverting a treated case back to clinical | Folded into $`f_t^{\text{eff}}(t) = f_t\cdot\text{eff}\cdot(1-\text{ETF}(t))`$ | None |
| **Bed net ownership, net age, retention, repeat rounds** | Per-person net-receipt timestep; efficacy decays from that date; nets lost stochastically | Two scalars per species, $`a_s(t)`$ and $`\mu_s(t)`$, computed from the full mean-field mixture over *all past distributions* (each weighted by most-recent-receipt probability $`\times`$ retention survival, each decaying its own $`d_{n0}`$ / $`r_n`$) | Protection is population-averaged rather than correlated within individuals across bites, so nets tend to **under**-suppress transmission relative to the IBM |
| **IRS spray status and spray age** | Per-person house-spray timestep | Folded into the same $`a_s(t)`$, $`\mu_s(t)`$; mixture over all past rounds weighted by most-recent-spray probability, with no retention factor (sprayed protection never expires in the IBM) | As above |
| **Vaccination status, dose count, booster stratum, antibody titre** | Per-person dose history and last-vaccination date. Antibody parameters are *not* stored — the IBM re-draws all four from their profile distributions every timestep, for every vaccinated person | One $`[n_a \times n_t]`$ multiplier $`v_i(t)`$ on the infection hazard. The vaccinated are partitioned analytically by the most recent booster each has reached; each stratum’s efficacy is the **average of the Hill efficacy over the IBM’s own 4-D antibody distribution** (§4.13) | The antibody average is exact, since the IBM’s draws are independent. What is lost is the dose *history*: vaccinated status is uncorrelated with anything else (nets, treatment), `min_wait` re-vaccination exclusion is not applied (§4.15), and seasonal boosters become a fixed days-since-primary schedule (warned) |
| **TBV antibody titre** | Per-person, mapped to transmission-reducing activity | Four $`[n_a \times n_t]`$ infectivity multipliers, one per infectious state ($`U`$, $`A`$, $`D`$, $`T`$), because the IBM maps transmission-reducing activity (an antibody-level quantity) to transmission-blocking activity (the per-state reduction in infectivity) with a state-specific transform | Mean-field only |
| **Immunity boosting refractory clock** | Per-person `last_boosted` timestep; a boost fires only if $`\ge u`$ days have passed | An analytic renewal rate $`q/(q\,u_{\text{eff}} + 1)`$ with $`u_{\text{eff}} = \lceil u\rceil - 1`$ and $`q`$ the *deduplicated per-day* event probability | Exact for the mean inter-boost interval |
| **Age-specific mortality / demographic transition** | Per-person death hazard by age band and year | $`\mu_i(t)`$, a constant-interpolated $`[n_a \times n_t]`$ series; the $`t = 0`$ row also sets the equilibrium age structure of the seed | Ages above the top `deathrate_agegroups` take the top rate (the IBM removes them); a warning fires |
| **Maternal immunity** | Sampled from a mother drawn at birth | Purely algebraic: $`I_{CM,ij} = P_M \cdot \overline{I_{CA}}_{20,j}\cdot\kappa_i`$ with $`\kappa_i`$ the exact within-band average of $`e^{-a/d_m}`$ — **no state variable at all** | None beyond the mean-field average over mothers |
| **Seasonality** | Daily rainfall drives larval carrying capacity | $`K_s(t)`$, linearly interpolated on a daily grid, from the same truncated Fourier series | None (identical functional form) |
| **Repeated bites in one day** | Bitten individuals are collected in a set, so a person bitten several times counts once and can be infected at most once per day | Saturating hazard $`-\log(1-(1-e^{-\varepsilon})b)`$ instead of the linear $`b\varepsilon`$ | None at the daily-flow level; it is the *linear* form that would be wrong |

### 3.3 What is captured as a discrete event

Mass drug administration cannot be written as a smooth rate: it is a
short, near-instantaneous clearance of a target age band. `blink` splits
the integration at every scheduled MDA/SMC/PMC timestep and applies a
**state jump** (§E), then resumes. This keeps chemoprevention out of the
ODE right-hand side entirely.

### 3.4 What is not captured at all

| Feature | Why | Behaviour |
|----|----|----|
| *P. vivax* | A different model (hypnozoites, relapse) | Error at input |
| Individual mosquitoes (`individual_mosquitoes = TRUE`) | `blink` is always compartmental | Ignored |
| Intervention correlation (`get_correlation_parameters()`, `run_simulation(correlations =)`) | Purely individual-level; the mean field assumes independence between interventions | Warned and ignored |
| Per-individual stochastic variation | Deterministic model | By construction |
| Partner-drug resistance, late clinical / parasitological failure, reinfection during prophylaxis | `malariasimulation` itself rejects non-zero values for these upstream in `set_antimalarial_resistance()` | Cannot reach the model |

## 4. Parameter support, function by function

This section covers every user-facing `malariasimulation`
parameterisation function, argument by argument, as of
**`malariasimulation` v3.0.0**.

Arguments named `parameters` are omitted throughout — they are the list
being modified.

### Legend

| Mark | Meaning | What to do about it |
|----|----|----|
| ✅ **Exact** | Modelled exactly, or algebraically equivalently for the mean field | Nothing. Differences from the IBM will sit inside its Monte-Carlo noise. |
| 🟡 **Approx.** | Modelled, but the mean field loses something. The row says what. | Safe for comparing scenarios that share the approximation. Check against the IBM if the number you report is dominated by this mechanism. |
| ➖ **Inert** | Accepted without error, but has no effect on the run — by design, not by omission | Nothing. Don’t tune it. |
| ⚠ **Ignored** | Not modelled; the run continues (sometimes with a warning) | Decide whether you can live without it before trusting the run. |
| ⛔ **Rejected** | Errors at input, or cannot reach the model at all | Remove it from the parameter list. |

Where a row carries two marks, the second is a caveat on the first.

### Index

| Function                        | §    | Status |
|---------------------------------|------|--------|
| `get_parameters()`              | 4.1  | ✅     |
| `set_species()`                 | 4.2  | ✅     |
| `set_equilibrium()`             | 4.3  | ✅     |
| `set_demography()`              | 4.4  | ✅     |
| `set_drugs()`                   | 4.5  | 🟡     |
| `set_clinical_treatment()`      | 4.6  | 🟡     |
| `set_antimalarial_resistance()` | 4.7  | 🟡     |
| `set_bednets()`                 | 4.8  | 🟡     |
| `set_spraying()`                | 4.9  | 🟡     |
| `set_carrying_capacity()`       | 4.10 | 🟡     |
| `set_mda()`, `set_smc()`        | 4.11 | 🟡     |
| `set_pmc()`                     | 4.12 | 🟡     |
| `create_pev_profile()`          | 4.13 | ✅     |
| `set_pev_epi()`                 | 4.14 | 🟡     |
| `set_mass_pev()`                | 4.15 | 🟡     |
| `set_tbv()`                     | 4.16 | 🟡     |
| `set_epi_outputs()`             | 4.17 | ✅     |
| `set_parameter_draw()`          | 4.18 | ✅     |
| Helpers and run functions       | 4.19 | mixed  |

### 4.1 `get_parameters(overrides, parasite)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `overrides` | Named list overriding any default parameter | ✅ Passed through wholesale; see the group-by-group breakdown below |
| `parasite` | `"falciparum"` or `"vivax"` | ✅ falciparum. ⛔ vivax — `build_inputs()` stops with an error |

#### Breakdown of `overrides`

| Group | Parameters | `blink` |
|----|----|----|
| Initial state proportions | `s_proportion`, `d_proportion`, `a_proportion`, `u_proportion`, `t_proportion` | ⛔ Ignored — `blink` seeds at the `malariaEquilibrium` fixed point for `init_EIR`, not at user-supplied proportions |
| Initial immunity | `init_ib`, `init_ica`, `init_iva`, `init_icm`, `init_ivm`, `init_id` | ⛔ Ignored — same reason |
| Population size | `human_population` | ✅ Scales output counts only; the ODE is per-capita, so run time is independent of it |
|  | `human_population_timesteps` | ⛔ Ignored — population is conserved; no dynamic population size |
| Baseline demography | `average_age` | ✅ → $`\eta = 1/\text{average_age}`$: the constant death hazard and the equilibrium age structure when `custom_demography = FALSE` |
|  | `custom_demography` | ✅ Switches to the `set_demography()` path (§4.4) |
| Biting heterogeneity | `a0`, `rho` | ✅ $`\psi_i = 1 - \rho e^{-a_i/a_0}`$ |
|  | `sigma_squared` | ✅ Log-normal variance of $`\zeta`$ |
|  | `n_heterogeneity_groups` | ✅ Number of Gauss–Hermite nodes $`n_z`$ |
|  | `enable_heterogeneity` | ✅ `FALSE` collapses to a single node with $`\zeta = 1`$ |
| Aquatic mosquito | `del`, `dl`, `dpl`, `me`, `ml`, `mup`, `gamma` | ✅ Used verbatim in the $`E`$/$`L`$/$`P_L`$ equations |
| Adult mosquito | `mum`, `beta`, `total_M`, `blood_meal_rates`, `Q0`, `foraging_time`, `species`, `species_proportions` | ✅ Per-species; `total_M` is re-derived (§4.3) |
|  | `init_foim` | ⛔ Ignored — `blink` recomputes FOIM from its own seeded human infectivity, which must be self-consistent with its seed |
| Seasonality | `model_seasonality`, `g0`, `g`, `h`, `rainfall_floor` | ✅ Truncated Fourier rainfall drives $`K_s(t) = K_{0,s}\,\text{scaler}(t)\,\text{rainfall}(t)/\bar{R}`$, on a daily grid |
| Flexible carrying capacity | `carrying_capacity`, `carrying_capacity_timesteps` | ✅ See §4.10. (`carrying_capacity_scalers` is *not* a `get_parameters()` default and cannot be passed in `overrides` — it is created only by `set_carrying_capacity()`.) |
|  | `carrying_capacity_values` | ⛔ Not read (nor by `malariasimulation`; only the `_scalers` form is used) |
| Human disease & immunity constants | `dd`, `dt`, `da`, `du` | ✅ Back-translated to $`r_D, r_T, r_A, r_U`$, then converted to the IBM’s per-day exit probability $`1-e^{-1/d}`$ so the realised dwell matches |
|  | `rb`, `rc`, `rva`, `rid`, `rm`, `rvm` | ✅ Immunity decay time constants $`d_{I_B}, d_{I_{CA}}, d_{I_{VA}}, d_{I_D}`$ and the maternal decay $`d_m, d_{vm}`$ |
|  | `ub`, `uc`, `uv`, `ud` | ✅ Refractory periods, applied as $`u_{\text{eff}} = \lceil u\rceil - 1`$ (§B.5) |
|  | `pcm`, `pvm` | ✅ $`P_M`$, $`P_{VM}`$ maternal transfer fractions |
|  | `b0`, `b1`, `ib0`, `kb` | ✅ $`b`$ Hill function |
|  | `phi0`, `phi1`, `ic0`, `kc` | ✅ $`\phi`$ Hill function |
|  | `theta0`, `theta1`, `iv0`, `kv`, `fv0`, `av`, `gammav` | ✅ $`\theta`$ Hill function and its age modifier $`f_v`$ |
|  | `d1`, `id0`, `kd`, `fd0`, `ad`, `gammad` | ✅ $`q`$ Hill function and its age modifier $`f_d`$ |
|  | `cd`, `cu`, `ct`, `gamma1` | ✅ Infectivity by state; `ct` is superseded by the drug-linked $`c_T(t)`$ whenever clinical drugs are set |
|  | `de`, `delay_gam`, `dem` | ✅ $`\tau_E`$, $`\tau_l`$, $`\tau_M`$ — the three Erlang chains |
| Vector-control shape | `phi_bednets`, `phi_indoors`, `k0` | ✅ Per species, in the $`a(t)`$/$`\mu(t)`$ algebra |
| Output bands | `age_group_rendering_*`, `incidence_rendering_*`, `clinical_incidence_rendering_*`, `severe_incidence_rendering_*`, `prevalence_rendering_*` | ✅ Each drives its own output family, exactly as in `malariasimulation` (§G) |
|  | `ib_`/`id_`/`ica_`/`iva_`/`idm_`/`icm_`/`ivm_rendering_*` | ⛔ `blink` does not render mean immunity by band |
| Metapopulation test-and-treat | `rdt_intercept`, `rdt_coeff` | ⛔ Ignored — they parameterise the PCR→RDT conversion used by `run_metapop_simulation()`’s mixing, which `blink` does not support (§4.19) |
| Engine / solver | `mosquito_limit`, `individual_mosquitoes`, `r_tol`, `a_tol`, `ode_max_steps`, `progress_bar` | ⛔ Ignored — `blink` is always compartmental and exposes its own `atol`/`rtol`/`step_size_max` on [`run_simulation_ode()`](https://pwinskill.github.io/blink/reference/run_simulation_ode.md) |
| Vivax-only | All hypnozoite parameters, `drug_hypnozoite_*`, `n_with_hypnozoites_rendering_*` | ⛔ Not applicable |

#### `blink`-only extension fields

Neither is a `malariasimulation` parameter; both may be added to the
list and both default to the validated choice.

| Field | Default | Effect |
|----|----|----|
| `bite_dedup` | `1` | `1` reproduces the IBM’s per-timestep bite deduplication (saturating hazard, §B.2). `0` uses the linear $`b\varepsilon`$ that `malariaEquilibrium` assumes, cutting the residual drift off the seed to well under 1% — but not to zero, since it gates only the hazard and not the boosting or sojourn forms |
| `acquired_immunity_offset` | `0` | $`\delta`$ in the $`b`$/$`\phi`$/$`\theta`$ Hill functions. `0.5` reproduces the IBM’s literal per-individual offset; `0` is the better mean-field match against the IBM ensemble mean |

### 4.2 `set_species(species, proportions)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `species` | List of species parameter lists (`gamb_params`, `arab_params`, `fun_params`, `kol_params`, `steph_params`, or custom) | ✅ Each entry’s `blood_meal_rates`, `foraging_time`, `Q0`, `phi_bednets`, `phi_indoors` and `mum` become per-species vectors driving an independent mosquito sub-model, coupled only through the shared human population |
| `proportions` | Relative abundance, summing to 1 | ✅ Splits `total_M` and hence each species’ baseline carrying capacity. A species at proportion 0 stays inert; its $`K_0`$ is floored to a negligible positive value so the larval term does not evaluate $`0/0`$ |

### 4.3 `set_equilibrium(init_EIR, eq_params, EIR_population_input)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `init_EIR` | Target EIR to seed from | ✅ Stored as `parameters$init_EIR` and read by [`run_simulation_ode()`](https://pwinskill.github.io/blink/reference/run_simulation_ode.md) when its own `init_EIR` argument is `NULL` |
| `eq_params` | Custom `malariaEquilibrium` parameter set | ✅ Merged over `blink`’s own back-translation and takes precedence. Note `set_equilibrium()` writes this field **even when passed `NULL`** — storing its own back-translation — so after any call to it, that stored set, not `R/translate_params.R`, supplies every shared constant (§F) |
| `EIR_population_input` | `"adult"` (default) or `"total"` | ✅ Handled upstream — `malariasimulation` converts a total-population EIR to adult EIR before storing it, and `blink` always interprets `init_EIR` as adult EIR (bites per adult per year) |

`set_equilibrium()` also calls `parameterise_mosquito_equilibrium()`
internally. That part has **no effect** on `blink`: `parameters$total_M`
is never read, and every species’ baseline carrying capacity is computed
from a `total_M` that `blink` re-derives so that
$`\sum_s a_s I_{M,s} = \text{init_EIR}/365`$ holds exactly under its own
seeded human infectivity. Calling `parameterise_total_M()` or
`parameterise_mosquito_equilibrium()` directly is harmless but changes
nothing.

### 4.4 `set_demography(agegroups, timesteps, deathrates)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `agegroups` | Upper edges of the death-rate age groups, in days | ✅ Model age groups are binned into them by midpoint using right-closed intervals, matching `.bincode(age, c(0, agegroups))` |
| `timesteps` | When each death-rate row takes effect | ✅ Knots of a **constant**-interpolated $`\mu_i(t)`$, so custom demography is time-varying (a demographic transition is modelled, not frozen) |
| `deathrates` | `[length(timesteps) × length(agegroups)]` daily death rates | ✅ Used directly as $`\mu_i(t)`$. The $`t = 0`$ row additionally sets the equilibrium age structure the human seed is rescaled onto. 🟡 Model ages above the top `agegroups` edge take the top rate here, whereas the IBM removes them — a warning fires; raise `default_age_lower(max_age =)` or extend `agegroups` to match |

### 4.5 `set_drugs(drugs)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `drugs` | List of 4-element vectors `c(efficacy, rel_c, prophylaxis_shape, prophylaxis_scale)` (e.g. `AL_params`, `DHA_PQP_params`, `SP_AQ_params`) | 🟡 Each element is handled as below |
|   `drug_efficacy` | Probability treatment clears the infection | ✅ Multiplies coverage into $`f_t^{\text{eff}}`$; also multiplies chemoprevention pulse fractions |
|   `drug_rel_c` | Infectivity of a treated case relative to `cd` | ✅ $`c_T = c_d \times \text{rel_c}`$, coverage-share-weighted across active drugs and time-varying |
|   `drug_prophylaxis_shape`, `drug_prophylaxis_scale` | Weibull prophylaxis hazard | 🟡 Represented by the Weibull **mean** $`\bar{d}_W`$: $`r_P = 1/(\bar{d}_W - 1/r_T)`$, subtracting the days already spent refractory in $`T`$. Equilibrium-exact, transient-approximate (exponential rather than Weibull decay) |
|   `drug_hypnozoite_*` | *P. vivax* only | ⛔ Not applicable |

### 4.6 `set_clinical_treatment(drug, timesteps, coverages)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `drug` | Index into the drug table | ✅ Identifies which drug’s efficacy / `rel_c` / prophylaxis enter the blend |
| `timesteps` | When coverage changes | ✅ Knots of the step-interpolated $`f_t(t)`$, and of the drug-mix series |
| `coverages` | Fraction of clinical cases treated with this drug | ✅ $`f_t(t) = \min\left(1, \sum_{\text{drugs}} \text{cov}_d(t)\right)`$. `set_clinical_treatment()` *errors* if the summed coverage exceeds 1 at any timestep, so the `min(1, ·)` `blink` applies is a defensive floor that never binds on a valid parameter set. 🟡 The drug-linked quantities $`\text{eff}(t)`$, $`c_T(t)`$, $`r_P(t)`$ are blended by each drug’s *instantaneous* coverage share, so a first-line **switch** is modelled; a genuinely mixed first line is represented by its blend rather than as parallel sub-populations |

### 4.7 `set_antimalarial_resistance(...)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `drug` | Which clinical drug carries resistance | ✅ Matched to the clinical-treatment drug index; also matched against the MDA/SMC/PMC drug so a resistant chemoprevention drug clears fewer infections |
| `timesteps` | When resistance levels change | ✅ Knots of $`\text{ETF}(t)`$ and $`\text{SPC}(t)`$, merged with the treatment-coverage change times |
| `artemisinin_resistance_proportion` | Fraction of infections that are artemisinin-resistant | ✅ Multiplies both the ETF and SPC probabilities |
| `partner_drug_resistance_proportion` | Partner-drug resistance | ⛔ `malariasimulation` requires 0; cannot reach the model |
| `slow_parasite_clearance_probability` | Probability a treated resistant case clears slowly | ✅ $`\text{SPC}(t)`$ — the split of the treated inflow between $`T`$ and $`T_s`$. Coverage-share-weighted across drugs |
| `early_treatment_failure_probability` | Probability treatment fails early | ✅ $`\text{ETF}(t)`$ — reduces $`f_t^{\text{eff}}`$, diverting those cases to $`D`$ |
| `late_clinical_failure_probability` | Late clinical failure | ⛔ `malariasimulation` requires 0 |
| `late_parasitological_failure_probability` | Late parasitological failure | ⛔ `malariasimulation` requires 0 |
| `reinfection_during_prophylaxis_probability` | Reinfection while prophylactic | ⛔ `malariasimulation` requires 0 |
| `slow_parasite_clearance_time` | Mean duration of slow clearance | 🟡 $`r_T^{\text{slow}} = 1/\text{dt_slow}`$, a **single scalar**. With one resistant drug the blend is exact; with several it is the coverage-weighted blend evaluated at *peak* resistance (peak, not final, so a rise-then-fall schedule is not collapsed to zero). Note this rate does **not** carry the whole-day $`1-e^{-1/d}`$ conversion the other sojourns get, so the realised slow-clearance dwell runs ~1 day short of the IBM’s (§B.4) |

### 4.8 `set_bednets()`

`set_bednets(parameters, timesteps, coverages, dn0, rn, rnm, gamman, retention, logistic_half_life, logistic_k)`

**Module status: 🟡** — net efficacy is averaged over the net-using
fraction into per-species $`a_s(t)`$ and $`\mu_s(t)`$, so protection is
not correlated within individuals across bites. Within that, every
argument below is handled exactly.

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `timesteps` | Distribution dates | ✅ Every past distribution contributes: at each grid time the population is the mixture over all rounds, each weighted by most-recent-receipt probability $`\text{cov}_d\prod_{j>d}(1-\text{cov}_j)`$ |
| `coverages` | Fraction receiving a net at each distribution | ✅ As above; repeated distributions accumulate correctly rather than only the latest applying |
| `dn0` | `[n_timesteps × n_species]` probability a net kills a mosquito | ✅ Decayed by net age: $`d_{n}(s_n) = d_{n0}e^{-s_n/\gamma_n}`$ |
| `rn` | `[n_timesteps × n_species]` initial repelling probability | ✅ Decayed toward `rnm`: $`r_n(s_n) = (r_n - r_{nm})e^{-s_n/\gamma_n} + r_{nm}`$ |
| `rnm` | `[n_timesteps × n_species]` minimum (asymptotic) repelling probability | ✅ As above |
| `gamman` | Insecticide decay time constant per distribution | ✅ The $`\gamma_n`$ above; each distribution decays on its own clock |
| `retention` | Mean net-retention time (log-uniform model) | ✅ Exponential usage survival $`e^{-s_n/\text{retention}}`$ |
| `logistic_half_life`, `logistic_k` | Alternative logistic retention | ✅ Exact survival of the IBM’s logistic retention time: $`S(s_n) = \exp\!\left(-k\frac{r^2}{1-r^2}\right)`$ for $`r = s_n/l < 1`$ (else 0), with $`l = \text{half_life}/\sqrt{1 - k/(k-\log 0.5)}`$. Verified $`S(\text{half_life}) = 0.5`$ exactly |

### 4.9 `set_spraying()`

`set_spraying(parameters, timesteps, coverages, ls_theta, ls_gamma, ks_theta, ks_gamma, ms_theta, ms_gamma)`

**Module status: 🟡** — as for nets: population-averaged, exact within.
Sprayed protection never expires in the IBM, so the mixture over past
rounds carries no retention factor.

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `timesteps` | Spray dates | ✅ Mixture over all past rounds, weighted by most-recent-spray probability |
| `coverages` | Fraction of houses sprayed | ✅ As above |
| `ls_theta`, `ls_gamma` | `[n_timesteps × n_species]` mortality logistic parameters | ✅ $`l_s(s_s) = \text{logit}^{-1}(\theta + \gamma s_s)`$ |
| `ks_theta`, `ks_gamma` | Feeding-success logistic parameters | ✅ $`k_s = k_0\,\text{logit}^{-1}(\theta + \gamma s_s)`$ |
| `ms_theta`, `ms_gamma` | Deterrence logistic parameters | ✅ $`m_s = \text{logit}^{-1}(\theta + \gamma s_s)`$ |

The repellency $`r_s`$ and survival $`s_s`$ derived from these are the
**same** spray outcome per individual, so `blink` forms the joint mean
$`\overline{(1-r_s)s_s}`$ rather than multiplying two independent
averages. Because the IBM recomputes these every timestep,
$`a_s(t)`$/$`\mu_s(t)`$ are linearly interpolated so within-round
logistic decay is a ramp, not a step.

### 4.10 `set_carrying_capacity(timesteps, carrying_capacity_scalers)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `timesteps` | When each scaler row takes effect | ✅ Step change in $`K_s(t)`$, composed multiplicatively with seasonality |
| `carrying_capacity_scalers` | `[n_timesteps × n_species]` multipliers on baseline $`K_0`$ | ✅ Multiplies baseline $`K_0`$ per species. 🟡 A scaler of exactly 0 is floored at $`K_0 \times 10^{-4}`$ to keep the larval term finite, so complete vector elimination is approached but not reproduced exactly |

### 4.11 `set_mda()` and `set_smc()`

`set_mda(parameters, drug, timesteps, coverages, min_ages, max_ages)` —
and `set_smc()` with the identical signature and identical handling.

**Module status: 🟡** — a mass campaign is applied as a pulse rather
than as per-individual events, and the brief treated-infectious phase of
cleared cases is omitted (negligible for the fast-clearing drugs used).

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `drug` | Drug administered | ✅ Its `drug_efficacy` scales the cleared fraction; its Weibull mean sets $`r_{P_c}`$. Antimalarial-resistance ETF on this drug reduces the cleared fraction |
| `timesteps` | Round dates | ✅ Integration is split at each; the pulse takes effect the day **after** the scheduled timestep |
| `coverages` | Fraction of the target band reached | ✅ Cleared fraction $`= \text{coverage}\times\text{efficacy}\times(1-\text{ETF})\times o_i`$. A zero-coverage round is a no-op |
| `min_ages`, `max_ages` | Per-round vectors (one entry per `timesteps` entry) giving the target band in days, **inclusive at both ends** | ✅ Per-round bands are honoured. Mapped onto the model age grid by **fractional overlap** $`o_i`$, so narrow bands are not dropped and coarse groups are not over-treated. The absorbing top group is treated as fully covered iff the band reaches its lower edge. 🟡 `blink` treats the band as half-open $`[\ell, u)`$, one day narrower than the IBM — negligible except for very narrow (PMC-style) bands |

Co-deployed chemoprevention types (SMC + MDA + PMC) share the single
$`P_c`$ compartment. Its decay rate is the mean of their drugs’
protection durations, weighted by each type’s **total scheduled
coverage** — the sum over all of its rounds — so a type with more rounds
carries proportionally more weight regardless of its per-round coverage.

### 4.12 `set_pmc(drug, timesteps, coverages, ages)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `drug` | Drug administered | ✅ As for MDA/SMC |
| `timesteps` | Coverage-change schedule (not dose dates — PMC delivery is age-triggered) | ✅ Read as a coverage schedule |
| `coverages` | Fraction of infants receiving each dose | ✅ The coverage in force at each pulse time is used |
| `ages` | Ages (days) at which doses are given | 🟡 The IBM triggers on an individual reaching each dose age. `blink` approximates this as pulses at a ~30-day cadence over a 30-day band starting at each dose age, with per-group overlap weighting, keeping the effective dose at approximately coverage × efficacy |

### 4.13 `create_pev_profile(vmax, alpha, beta, cs, rho, ds, dl)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `vmax` | Maximum efficacy | ✅ |
| `alpha`, `beta` | Hill shape and scale mapping antibody titre to efficacy | ✅ |
| `cs` | `c(mu, sigma)` of $`\log`$ peak antibody titre | ✅ **Including the sigma** |
| `rho` | `c(mu, sigma)` of the logit short-lived antibody fraction | ✅ Including the sigma |
| `ds`, `dl` | `c(mu, sigma)` of $`\log`$ short/long antibody half-lives | ✅ Including the sigmas |

`malariasimulation` re-draws all four antibody parameters from their
profile distributions **at every timestep, for every vaccinated
individual** — nothing is stored per person. Population efficacy is
therefore exactly $`\mathbb{E}\left[\text{Hill}(\text{Ab})\right]`$, not
$`\text{Hill}(\text{Ab at the median})`$: the sigmas are large and the
Hill function is nonlinear, so the two differ materially (evaluating at
the median overstated R21 efficacy by up to ~4.7 percentage points).
`blink` integrates over the 4-D distribution with a tensor Gauss–Hermite
rule (`options(blink.pev_gq =)`, default 7 nodes per axis), validated
against a Monte Carlo of the IBM’s own sampler to under
$`5\times10^{-4}`$. Because the IBM’s draws are independent between days
and between people, there is no within-person correlation for the mean
field to lose here — the quadrature is exact rather than approximate.
The resulting mean-efficacy curve is cached per profile on a daily grid.

### 4.14 `set_pev_epi()`

`set_pev_epi(parameters, profile, coverages, timesteps, age, min_wait, booster_spacing, booster_coverage, booster_profile, seasonal_boosters)`

**Module status: 🟡** — the routine schedule and full booster sequence
are modelled, but seasonal boosters are approximated (warned).

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `profile` | Primary-series PEV profile | ✅ |
| `coverages` | Coverage of the routine programme over time | ✅ **Time-varying**: each cohort is protected at the coverage in force on *its own first-dose date*, not at a single coverage |
| `timesteps` | When each coverage takes effect | ✅ Also gates eligibility — a cohort is protected only once its first dose falls inside the programme. Gating follows `malariasimulation`, which keys on `pev_epi_coverages`/`pev_epi_timesteps` and never reads `parameters$pev` |
| `age` | Age (days) at the first dose | ✅ Efficacy starts at `age + max(pev_doses)`, i.e. after the final primary dose |
| `min_wait` | Minimum time since the last vaccination | ➖ For a routine EPI programme this guards re-vaccination and seasonal-booster timing, not the primary-plus-booster schedule `blink` models, so it has no effect on the result. (For **mass** campaigns it does bite — see §4.15.) |
| `booster_spacing` | Days from the final primary dose to each booster | ✅ The **full** booster sequence is modelled: the vaccinated are partitioned by the most recent booster each has reached, and each stratum carries its own decayed efficacy |
| `booster_coverage` | Matrix of conditional booster coverages | ✅ Each booster’s coverage is read at *its own* administration date (matching `coverage[match_timestep(timesteps, admin_time), booster]`), reducing to the single-row lookup in the usual case |
| `booster_profile` | Per-booster profiles | ✅ Each booster uses its own profile, integrated over the antibody distribution as in §4.13 |
| `seasonal_boosters` | Time the first booster relative to the **start of the calendar year** — `booster_spacing[1]` becomes a day-of-year. It makes no reference to the fitted seasonality; you choose the day yourself, e.g. via `peak_season_offset()` | 🟡 Approximated as a fixed days-since-primary schedule. `blink` **warns** when this is set |

### 4.15 `set_mass_pev()`

`set_mass_pev(parameters, profile, timesteps, coverages, min_ages, max_ages, min_wait, booster_spacing, booster_coverage, booster_profile)`

**Module status: 🟡** — repeated campaigns are the weak point:
`min_wait` re-vaccination exclusion is not applied, campaigns combine as
independent protections rather than most-recent-receipt, and the
vaccinated cohort does not age out of its band. A single campaign is
modelled faithfully.

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `profile` | Primary-series profile | ✅ |
| `timesteps` | Campaign dates | 🟡 Every campaign contributes, combined as $`1-\prod_k(1-e_k)`$. The IBM instead keeps only the **most recent** vaccination per person, so its population protection is a most-recent-receipt mixture $`\sum_k \text{cov}_k \prod_{j>k}(1-\text{cov}_j)\,e_k`$ — the same structure `blink` already uses for bed nets (§4.8). The two agree for a single campaign and diverge when campaigns overlap |
| `coverages` | Coverage per campaign (recycled if length 1) | ✅ |
| `min_ages`, `max_ages` | Target age bands | ✅ **Every** band is applied at **every** campaign |
| `min_wait` | Minimum time since the last vaccination | ⚠ **Ignored, and it matters here.** In `malariasimulation` each campaign excludes anyone vaccinated within `min_wait` of it, so with repeated campaigns this is the primary control on who gets re-vaccinated. `blink` treats every campaign as reaching its whole target band independently, so a `min_wait` longer than the campaign spacing will over-vaccinate |
| `booster_spacing` | Days from the final primary dose to each booster | ✅ Full sequence, as in §4.14 |
| `booster_coverage` | Conditional booster coverage matrix | ✅ Read at each booster’s administration date |
| `booster_profile` | Per-booster profiles | ✅ |

🟡 The vaccinated cohort does not age out of its target band: efficacy
decays in place with time since vaccination, but the protected fraction
stays attached to the age band rather than moving up the age grid with
the cohort.

### 4.16 `set_tbv(timesteps, coverages, ages)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `timesteps` | Vaccination dates | 🟡 Antibody titre decays from each campaign date; campaigns combine as $`1-\prod_k(1-e_k)`$. As for mass PEV, the IBM overwrites each person’s vaccination date, giving a most-recent-receipt mixture instead. Identical for a single campaign; divergent when campaigns overlap |
| `coverages` | Fraction vaccinated | ✅ |
| `ages` | Whole years of age targeted | ✅ Exact year set: $`\lfloor a_i/365\rfloor \in \text{ages}`$, matching the IBM |

Transmission-reducing activity is mapped to state-specific
transmission-blocking activity via the IBM’s own transform
$`\text{TBA}(m_x, k, \text{TRA})`$ with per-state $`m_x \in
\{\texttt{tbv_mu}, \texttt{tbv_ma}, \texttt{tbv_md}, \texttt{tbv_mt}\}`$,
giving four independent per-age multipliers on $`U`$, $`A`$, $`D`$ and
$`T`$ infectivity. All the `tbv_*` shape constants (`tbv_tau`,
`tbv_rho`, `tbv_ds`, `tbv_dl`, `tbv_tra_mu`, `tbv_gamma1`, `tbv_gamma2`,
`tbv_k`) are used verbatim.

### 4.17 `set_epi_outputs(...)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `age_group` | Bands for population counts | ✅ → `n_age_*` (also emitted over the union of all other families’ bands, as `malariasimulation` does) |
| `incidence` | Bands for all-infection incidence | ✅ → `n_inc_*` |
| `clinical_incidence` | Bands for clinical incidence | ✅ → `n_inc_clinical_*` |
| `severe_incidence` | Bands for severe incidence | ✅ → `n_inc_severe_*` |
| `prevalence` | Bands for prevalence | ✅ → `n_detect_lm_*`, `p_detect_lm_*`, `n_detect_pcr_*` |
| `ica`, `icm`, `iva`, `ivm`, `id`, `ib` | Bands for mean immunity | ⛔ Not rendered by `blink` |
| `n_with_hypnozoites`, `hypnozoites`, `iaa`, `iam` | *P. vivax* only | ⛔ Not applicable |

When a family’s band list is empty, `blink` falls back to the
convenience bands 2–10 y and all-ages for that family only. A bare
`get_parameters()` leaves four of the five families empty; it defaults
the **prevalence** bands to 2–10 y.

### 4.18 `set_parameter_draw(draw)`

| Argument | malariasimulation meaning | `blink` |
|----|----|----|
| `draw` | Index 1–1000 into the fitted posterior draws | ✅ `malariasimulation` overwrites the immunity/disease constants in the parameter list in place; they then flow through `blink`’s back-translation like any other parameter. Must be called before `set_equilibrium()`, as in the IBM |

### 4.19 Helpers and run functions

| Function / argument | malariasimulation meaning | `blink` |
|----|----|----|
| `peak_season_offset(parameters)` | Day of peak seasonal transmission | ✅ Pure helper on `g0`/`g`/`h`, which `blink` reads; useful for timing SMC/booster schedules |
| `get_correlation_parameters(parameters)` | Correlation between intervention recipients | ⛔ Individual-level; no mean-field analogue |
| `run_simulation(timesteps, parameters, correlations)` | Run the IBM | ✅ `run_simulation_ode(timesteps, parameters, correlations, ...)` mirrors the first three arguments exactly. A non-`NULL` `correlations` is **warned and ignored** — the ODE assumes independence between interventions |
| `run_simulation_with_repetitions(...)` | Repeat stochastic runs | ⛔ Not applicable — `blink` is deterministic, so one run is the ensemble mean |
| `run_metapop_simulation(timesteps, parameters, ..., export_mixing, import_mixing)` | Coupled patches | ⛔ Not supported. `blink` is a single patch; there is no guard, so do not pass a metapopulation setup expecting it to be honoured |
| `run_resumable_simulation(...)` | Resume from saved state | ⛔ Not implemented. `blink` runs are cheap enough to re-run from the seed |

## 5. Using `blink` well: what to do differently

### Things to do

**Burn in before calibrating.** The seed is a close approximation, not
the exact fixed point (§F), and a seasonal run needs ~10 years to settle
onto its limit cycle. The seasonal annual-mean EIR sits a few percent
below the aseasonal `init_EIR` target through nonlinear averaging.
Discard the transient.

**Match the age grid to your demography.** If the oldest model age group
runs past the top `deathrate_agegroups`, `blink` applies the top death
rate to those ages while the IBM removes them. Fix it either way round —
widen the grid, or extend the death rates:

``` r

run_simulation_ode(t, p, init_EIR = 20,
                   age_lower = default_age_lower(max_age = 100))
```

A warning fires when they diverge; don’t ignore it, because it changes
the equilibrium age structure the whole run is seeded on (§4.4).

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

run_simulation_ode(t, p, init_EIR = 20,
                   atol = 1e-6, rtol = 1e-6, step_size_max = 10)
```

**Sanity-check conservation.**
`S_count + D_count + A_count + U_count + Tr_count + Ph_count` should
equal `human_population` to numerical tolerance.

### Things to leave alone

**`acquired_immunity_offset`.** Leave it at 0. Setting 0.5 reproduces
the IBM’s literal per-individual offset but is a *worse* match to the
IBM ensemble mean (§B.3).

**`bite_dedup`.** Leave it at 1 for any scientific run. Set 0 only for
an equilibrium test, where it cuts the drift off the seed to well under
1% (§B.2).

**Erlang stage counts (`n_eir`, `n_foim`, `n_eip`).** The defaults are
fine. Raising them sharpens the gamma-shaped lags toward the IBM’s fixed
delays at linear cost in state size; equilibrium is exact for any value.

**`parameterise_total_M()` / `parameterise_mosquito_equilibrium()`.**
They have no effect — `blink` re-derives `total_M` from `init_EIR`
(§4.3).

### Results to treat with caution

**Severe incidence.** The steep $`\theta`$ Hill function amplifies small
differences in acquired immunity and disease-pool composition between
the mean field and the IBM, so severe incidence — and the DALYs derived
from it — can differ by more than prevalence or clinical incidence does,
largest at low EIR.

**Anything through a nets or IRS era.** Population-averaged coverage
loses the correlation of protection within individuals across bites, so
the ODE tends to sit above the IBM while vector control is active
(§3.2).

**Repeated mass PEV or TBV campaigns.** Overlapping campaigns combine as
independent protections rather than most-recent-receipt, and `min_wait`
re-vaccination exclusion is not applied (§4.15, §4.16). A single
campaign is fine.

### Cost model

Run time scales with $`n_a \times n_z`$ (the human state), and —
importantly — with the **number of chemoprevention pulses**, because §E
splits the integration at every one. A 20-year monthly PMC schedule
means ~240 solver restarts; prefer coarser cadences when the extra
resolution isn’t buying you anything.

### Index of approximations

Every 🟡 in §4, and when it matters:

| Mechanism | § | Matters when |
|----|----|----|
| Weibull prophylaxis → exponential compartment | 4.5 | Looking at the shape of post-treatment protection, not its mean |
| Multi-drug blend rather than parallel sub-populations | 4.6 | A genuinely mixed first line (a *switch* is modelled fine) |
| Single scalar slow-clearance rate, no whole-day conversion | 4.7, B.4 | High artemisinin resistance with several treatment drugs |
| Population-averaged nets | 4.8 | Any nets era — the largest routine gap |
| Population-averaged IRS | 4.9 | Any IRS era |
| Carrying-capacity floor at $`K_0\times10^{-4}`$ | 4.10 | Modelling complete vector elimination |
| Chemoprevention as pulses; half-open age bands | 4.11, E | Very narrow target bands |
| PMC age trigger → monthly pulses | 4.12 | PMC dose timing specifically |
| Seasonal PEV boosters → fixed delay | 4.14 | Seasonally-timed booster campaigns (warned) |
| Repeated campaigns: no `min_wait`, independent combination, no ageing out of band | 4.15 | Repeated mass PEV |
| TBV campaign combination | 4.16 | Repeated TBV campaigns |

## Appendix: the equations

The rest of this article is reference material: the full system as
implemented, for readers who need the exact form of a term rather than
its shape.

### A. Notation

#### Abbreviations

|  |  |
|----|----|
| **EIR** | Entomological inoculation rate — infectious bites per person per year. `blink` reports it per *adult* per year |
| **FOI** / **FOIM** | Force of infection (on humans) / force of infection on mosquitoes |
| **EIP** | Extrinsic incubation period — the delay in the mosquito between an infectious blood meal and becoming infectious |
| **LM** / **PCR** | Light microscopy / polymerase chain reaction — the two prevalence diagnostics reported |
| **ETF** / **SPC** | Early treatment failure / slow parasite clearance — the two antimalarial-resistance arms `malariasimulation` implements |
| **PEV** / **TBV** | Pre-erythrocytic vaccine (RTS,S, R21) / transmission-blocking vaccine |
| **TRA** / **TBA** | Transmission-reducing activity (an antibody-level quantity) / transmission-blocking activity (the resulting per-state reduction in infectivity) |
| **MDA / SMC / PMC** | Mass drug administration / seasonal malaria chemoprevention / perennial malaria chemoprevention |
| **EPI** | The Expanded Programme on Immunization routine schedule — unrelated to `set_epi_outputs()`, where “epi” means *epidemiological* |
| **IRS** | Indoor residual spraying |
| **Erlang chain** | A linear chain of $`n`$ sequential exponential stages. Its total has mean $`\tau`$ and variance $`\tau^2/n`$, so larger $`n`$ sharpens it toward the IBM’s fixed delay |
| **Gauss–Hermite quadrature** | A deterministic stand-in for Monte Carlo: it replaces random per-person draws with a small fixed set of representative values and weights that reproduce the distribution’s moments. “Tensor” means the product grid over several such axes |
| **knots** | The times at which an interpolated series is allowed to change value |
| **seed** | Throughout this article, the *initial state vector* — never a random-number seed. `blink` has no random component |

#### Symbols

| Symbol | odin name | Meaning |
|----|----|----|
| $`r_i`$ | `r_age` | Aging rate out of age group $`i`$ = 1/band width in days ($`r_{n_a} = 0`$, absorbing) |
| $`\mu_i(t)`$ | `mu_age` | Age- and time-specific death rate |
| $`\rho_i(t) = r_i + \mu_i(t)`$ | `re` | Combined exit rate |
| $`\psi_i`$ | `psi` | Relative biting rate by age, $`1 - \rho\,e^{-a_i/a_0}`$ |
| $`\zeta_j`$, $`w_j`$ | `zeta`, `het_wt` | Heterogeneity node and quadrature weight |
| $`N_{ij}`$ | `Npop` | Total population density in cell $`(i,j)`$ |
| $`\varepsilon_{ij}`$ | `EPS` | Expected infectious bites per person per day |
| $`\Lambda_{ij}`$ | `FOI` | Force of infection (hazard per day) |
| $`h^c_{ij}`$, $`h^a_{ij}`$ | `h_c`, `h_a` | Clinical and non-clinical components of $`\Lambda`$ |
| $`b_{ij}`$ | `b` | Probability an infectious bite establishes infection |
| $`\phi_{ij}`$ | `phi` | Probability an infection is clinical |
| $`\theta_{ij}`$ | `theta` | Probability an infection is severe |
| $`q_{ij}`$ | `q` | Probability an asymptomatic infection is LM-detectable |
| $`c_A, c_D, c_U, c_T`$ | `cA`, `cD`, `cU`, `cT` | Onward infectivity to mosquitoes, by state |
| $`f_t(t)`$ | `ft` | Clinical treatment coverage |
| $`v_i(t)`$ | `pev_factor` | PEV hazard multiplier |
| $`a_s(t)`$, $`\mu_s(t)`$ | `a_spp`, `mum` | Human blood-meal rate and adult mortality, per species |
| $`K_s(t)`$ | `Kcap` | Larval carrying capacity |
| $`\tau_E, \tau_l, \tau_M`$ | `de`, `tl`, `dem` | EIR lag, gametocyte lag, extrinsic incubation period |

### B. The human model

#### B.1 Demography

Aging is a linear chain over the age grid; deaths are recycled as births
into age group 1, distributed across heterogeneity nodes by the
quadrature weights, so the total population is conserved:

``` math
N_{ij} = S_{ij}+D_{ij}+A_{ij}+U_{ij}+T_{ij}+T_{s,ij}+P_{ij}+P_{c,ij},
\qquad
B = \sum_{i,j}\mu_i(t)\,N_{ij}.
```

Every human compartment $`Y`$ carries the same structural terms

``` math
\dot{Y}_{ij} \;\supset\; r_{i-1}Y_{i-1,j} \;-\; \rho_i(t)\,Y_{ij},
```

with the aging inflow replaced by $`B\,w_j`$ for $`Y = S`$, $`i = 1`$,
and by zero otherwise.

#### B.2 Exposure and the infection hazard

The species-summed infectious biting rate, lagged by the Erlang chain
$`X^{(E)}`$, gives the population EIR; individual exposure scales it by
the age and heterogeneity biting weights:

``` math
\bar\varepsilon(t) = \sum_s a_s(t)\,I_{M,s}(t),
\qquad
\dot{X}^{(E)}_1 = \tfrac{n_E}{\tau_E}\!\left(\bar\varepsilon - X^{(E)}_1\right),
\quad
\dot{X}^{(E)}_k = \tfrac{n_E}{\tau_E}\!\left(X^{(E)}_{k-1} - X^{(E)}_k\right),
```

``` math
\varepsilon_{ij} = X^{(E)}_{n_E}\,\zeta_j\,\psi_i .
```

The IBM draws the day’s bites from a Poisson and collects the bitten in
a set, so a person bitten several times in one timestep can be infected
at most once. `blink` reproduces that cap; the daily infection
*probability* and the equivalent hazard are

``` math
p_{ij} = \left(1 - e^{-\varepsilon_{ij}}\right) b_{ij}\, v_i(t),
\qquad
\Lambda_{ij} = -\log\!\left(1 - p_{ij}\right),
```

matching the IBM’s `prob_to_rate`. Setting `parameters$bite_dedup = 0`
replaces this with the unbounded linear form
$`\Lambda_{ij} = b_{ij}\varepsilon_{ij}v_i`$, which is the hazard
`malariaEquilibrium` assumes. The two agree to under 2% at low exposure
and diverge at seasonal peaks and in high-$`\zeta`$ strata. Note that
`bite_dedup = 0` does **not** make the seed an exact fixed point: it
gates only the hazard, while the boosting terms of §B.5 keep the
deduplicated event probability, the integer refractory window
$`\lceil u\rceil - 1`$ and the at-risk restriction, and the state exit
rates keep the whole-day form $`1 - e^{-1/d}`$. What it does is reduce
the residual drift to well under 1%, which is why the package’s
equilibrium tests use it.

**Hazard splitting.** In the IBM each person resolves at most one
infection outcome per day, so the daily clinical count is $`\phi p N`$
for susceptibles. (For people already in $`A`$ or $`U`$ the infection
outcome competes with disease progression inside the *same* hazard
resolution, which scales their realised infection probability by
$`\frac{\Lambda}{\Lambda + r}\cdot\frac{1 - e^{-(\Lambda+r)}}{1 - e^{-\Lambda}}`$.
`blink`’s ODE reproduces that competition automatically, since
$`\Lambda`$ and $`r_A`$/$`r_U`$ act simultaneously on the same
compartment.) Clinical infections leave the at-risk pool (to $`D`$ or
$`T`$) while non-clinical ones do not ($`S, U \to A`$; $`A \to A`$), so
the pool depletes only through the clinical route. Taking

``` math
h^c_{ij} = -\log\!\left(1 - \phi_{ij}p_{ij}\right),
\qquad
h^a_{ij} = \Lambda_{ij} - h^c_{ij}
```

integrates to exactly $`\phi p N`$ clinical episodes over a day. **Why
not just $`\phi\Lambda`$?** An ODE re-exposes people continuously
through the day, so $`\phi\Lambda`$ would start a second episode in
someone the IBM counts once — over-counting by a factor
$`\Lambda/(1-e^{-\Lambda})`$. By construction $`h^c + h^a = \Lambda`$,
so occupancies and immunity boosting are unaffected.

#### B.3 Immunity-dependent probabilities

All four are Hill functions of the immunity states. $`\delta`$ is
`acquired_immunity_offset` (default 0; set 0.5 to reproduce the IBM’s
literal per-individual `+0.5`). The IBM applies that offset only where
acquired immunity is **strictly positive**, so never-boosted individuals
get none — one of the reasons a flat `+0.5` on a stratum mean is the
worse match, and why 0 is the default.

``` math
b_{ij} = b_0\left[b_1 + \frac{1-b_1}{1+\left(\frac{I_{B,ij}+\delta}{I_{B0}}\right)^{k_b}}\right]
```

``` math
\phi_{ij} = \phi_0\left[\phi_1 + \frac{1-\phi_1}{1+\left(\frac{I_{CA,ij}+\delta+I_{CM,ij}}{I_{C0}}\right)^{k_c}}\right]
```

``` math
\theta_{ij} = \theta_0\left[\theta_1 + \frac{1-\theta_1}{1+f_{v,i}\left(\frac{I_{VA,ij}+\delta+I_{VM,ij}}{I_{V0}}\right)^{k_v}}\right]
```

``` math
q_{ij} = d_1 + \frac{1-d_1}{1+\left(\frac{I_{D,ij}}{I_{D0}}\right)^{k_d} f_{d,i}},
\qquad
c_{A,ij} = c_U + (c_D-c_U)\,q_{ij}^{\gamma_{\text{inf}}}
```

with the age-dependent detection and severity modifiers

``` math
f_{d,i} = 1 - \frac{1-f_{d0}}{1+(a_i/a_{d0})^{\gamma_d}},
\qquad
f_{v,i} = 1 - \frac{1-f_{v0}}{1+(a_i/a_v)^{\gamma_v}} .
```

Note that $`q`$ takes **no** offset, matching the IBM.

**Maternal immunity** is algebraic, not a state. It is inherited from
mothers in the model age group containing 20 years (the IBM selects
$`\lfloor \text{age}/365 \rfloor = 20`$) and decays with age from birth:

``` math
\overline{I_{CA}}_{20,j} = \sum_i m_i I_{CA,ij},
\qquad
I_{CM,ij} = P_M\,\overline{I_{CA}}_{20,j}\,\kappa^{C}_i,
\qquad
I_{VM,ij} = P_{VM}\,\overline{I_{VA}}_{20,j}\,\kappa^{V}_i,
```

where $`m_i`$ is the 20-year indicator and $`\kappa^{C}_i`$ the exact
mean of $`e^{-a/d_m}`$ over age band $`i`$,

``` math
\kappa^{C}_i = \frac{d_m}{w_i}\left(e^{-a^-_i/d_m} - e^{-a^+_i/d_m}\right),
```

and $`\kappa^V_i`$ the same with $`d_{vm}`$.

#### B.4 Disease-state equations

Writing $`\Sigma_{ij} = S_{ij}+A_{ij}+U_{ij}`$ for the at-risk pool,
$`\text{SPC}(t)`$ for the slow-parasite-clearance fraction of the
treated (§4.7), and
$`f_t^{\text{eff}}(t) = f_t(t)\,\text{eff}(t)\,(1-\text{ETF}(t))`$ for
the effective treated fraction:

``` math
\dot{S}_{ij} = r_{i-1}S_{i-1,j} + r_U U_{ij} + r_P P_{ij} + r_{P_c} P_{c,ij}
- \Lambda_{ij}S_{ij} - \rho_i S_{ij}
```

``` math
\dot{D}_{ij} = r_{i-1}D_{i-1,j} + \left(1 - f_t^{\text{eff}}\right)h^c_{ij}\Sigma_{ij}
- r_D D_{ij} - \rho_i D_{ij}
```

``` math
\dot{T}_{ij} = r_{i-1}T_{i-1,j} + (1-\text{SPC})\,f_t^{\text{eff}}\,h^c_{ij}\Sigma_{ij}
- r_T T_{ij} - \rho_i T_{ij}
```

``` math
\dot{T}_{s,ij} = r_{i-1}T_{s,i-1,j} + \text{SPC}\,f_t^{\text{eff}}\,h^c_{ij}\Sigma_{ij}
- r_T^{\text{slow}} T_{s,ij} - \rho_i T_{s,ij}
```

``` math
\dot{A}_{ij} = r_{i-1}A_{i-1,j} + h^a_{ij}\left(S_{ij}+U_{ij}\right)
- h^c_{ij}A_{ij} + r_D D_{ij} - r_A A_{ij} - \rho_i A_{ij}
```

``` math
\dot{U}_{ij} = r_{i-1}U_{i-1,j} + r_A A_{ij} - \Lambda_{ij}U_{ij} - r_U U_{ij} - \rho_i U_{ij}
```

``` math
\dot{P}_{ij} = r_{i-1}P_{i-1,j} + r_T T_{ij} + r_T^{\text{slow}} T_{s,ij}
- r_P P_{ij} - \rho_i P_{ij}
```

``` math
\dot{P}_{c,ij} = r_{i-1}P_{c,i-1,j} - r_{P_c}P_{c,ij} - \rho_i P_{c,ij}
```

Three implementation points:

- **Two treated compartments.** The IBM assigns each successfully
  treated individual to standard or slow parasite clearance by a
  Bernoulli draw, giving a *mixture of two exponentials*, not one
  exponential at the blended mean. $`T`$ and $`T_s`$ reproduce that
  split structurally (but see the sojourn caveat below for
  $`r_T^{\text{slow}}`$).
- **Whole-day sojourns.** The IBM advances states once per day with exit
  probability $`1 - e^{-1/d}`$, so its realised mean dwell is
  $`1/(1-e^{-1/d})`$ (5.52 d for $`d_D = 5`$, not 5). `blink` uses
  $`r_X = 1 - e^{-1/d_X}`$ for $`X \in \{A, D, U, T\}`$ so the dwells
  match. Two caveats. It is *not* applied to $`r_P`$, because the IBM
  models prophylaxis as a hazard multiplier rather than a compartment
  with a daily census. It is also *not* currently applied to
  $`r_T^{\text{slow}}`$, which is set to $`1/\text{dt_slow}`$ — whereas
  the IBM resolves slow clearance through the same per-day
  `rate_to_prob` path as ordinary treatment, giving it a realised dwell
  of $`1/(1-e^{-1/\text{dt_slow}})`$. Slow-clearance dwells therefore
  run about a day short of the IBM’s.
- **$`P_c`$ has no ODE inflow.** It is filled only by the
  chemoprevention pulses of §E, and carries neither infection risk nor
  infectivity.

#### B.5 Immunity equations

The IBM boosts immunity on a per-day *event* with an integer refractory
window: `boost_immunity()` fires only if `timestep - last_boosted >= u`,
and the timestep is an integer, so the effective wait is
$`\lceil u\rceil`$ days, after which the person boosts on the first
subsequent day carrying an event (geometric, mean $`1/q`$). The mean
inter-boost gap is therefore $`(\lceil u\rceil - 1) + 1/q`$, i.e. a
boosting rate

``` math
\beta(q, u) = \frac{q}{q\,u_{\text{eff}} + 1},
\qquad u_{\text{eff}} = \lceil u \rceil - 1 .
```

The event probability $`q`$ is the **deduplicated per-day** one, and
differs between $`I_B`$ and the rest:

``` math
q^b_{ij} = 1 - e^{-\varepsilon_{ij}} \quad\text{(bitten at least once)},
\qquad
q^f_{ij} = 1 - e^{-\Lambda_{ij}} \quad\text{(infected that day)}.
```

$`I_{CA}`$, $`I_D`$ and $`I_{VA}`$ are boosted only for individuals
eligible to be infected (the IBM restricts its source set to $`S`$,
$`A`$, $`U`$), so their boost is scaled by the at-risk fraction
$`\alpha_{ij} = \Sigma_{ij}/N_{ij}`$. $`I_B`$ is boosted for everyone
bitten, whatever their state, so it is not scaled:

``` math
\dot{I}_{B,ij} = \beta(q^b_{ij}, u_b) - \frac{I_{B,ij}}{d_{I_B}} + \rho_i\left(I_{B,i-1,j} - I_{B,ij}\right)
```

``` math
\dot{I}_{X,ij} = \alpha_{ij}\,\beta(q^f_{ij}, u_X) - \frac{I_{X,ij}}{d_{I_X}} + \rho_i\left(I_{X,i-1,j} - I_{X,ij}\right),
\qquad X \in \{CA, D, VA\}
```

with $`I_{X,0,j} \equiv 0`$ (immunity enters at zero at birth).

### C. Human → mosquito infectivity

Onward infectivity is the state-weighted, biting-weighted human
infectiousness, passed through the gametocyte-delay Erlang chain
$`X^{(F)}`$. TBV enters here as per-state, per-age multipliers
$`\tau^U_i, \tau^A_i, \tau^D_i, \tau^T_i`$:

``` math
\text{inf}_{ij} = c_D D_{ij}\tau^D_i + c_{A,ij} A_{ij}\tau^A_i + c_U U_{ij}\tau^U_i
+ c_T(t)\left(T_{ij}+T_{s,ij}\right)\tau^T_i
```

``` math
\overline{\text{inf}} = \frac{1}{\bar\psi}\sum_{i,j}\zeta_j\psi_i\,\text{inf}_{ij},
\qquad
\dot{X}^{(F)}_k = \tfrac{n_F}{\tau_l}\left(X^{(F)}_{k-1} - X^{(F)}_k\right),
\qquad
\Lambda^M_s = a_s(t)\,X^{(F)}_{n_F}.
```

### D. The mosquito model

Per species $`s`$, with adult total
$`M_s = S_{M,s} + E_{M,s} + I_{M,s}`$ and larval total
$`n_{L,s} = E_s + L_s`$:

``` math
\dot{E}_s = \beta_s M_s - \frac{E_s}{d_E} - \mu_E E_s\left(1 + \frac{n_{L,s}}{K_s(t)}\right)
```

``` math
\dot{L}_s = \frac{E_s}{d_E} - \frac{L_s}{d_L} - \mu_L L_s\left(1 + \gamma\frac{n_{L,s}}{K_s(t)}\right)
```

``` math
\dot{P}_{L,s} = \frac{L_s}{d_L} - \frac{P_{L,s}}{d_P} - \mu_P P_{L,s}
```

``` math
\dot{S}_{M,s} = \tfrac{1}{2}\frac{P_{L,s}}{d_P} - S_{M,s}\Lambda^M_s - \mu_s(t)S_{M,s}
```

**Extrinsic incubation.** `malariasimulation` uses a *fixed*
$`\tau_M`$-day delay and applies the incubation survival
$`e^{-\mu\tau_M}`$ at the exit, with the current $`\mu`$. `blink`
substitutes an Erlang chain for the delay, but the chain must be
**loss-free**: putting death in every stage would give a
through-survival of
$`\left(\frac{n_P/\tau_M}{n_P/\tau_M + \mu}\right)^{n_P}`$ instead of
$`e^{-\mu\tau_M}`$ — about 4% too high at baseline mortality, and
worsening as vector control raises $`\mu`$. So the chain is a pure delay
of $`S_M\Lambda^M`$, survival is applied once at the exit, and the
incubating stock follows the IBM’s own $`E`$ equation:

``` math
\dot{X}^{(P)}_{s,1} = \tfrac{n_P}{\tau_M}\left(S_{M,s}\Lambda^M_s - X^{(P)}_{s,1}\right),
\qquad
\dot{X}^{(P)}_{s,k} = \tfrac{n_P}{\tau_M}\left(X^{(P)}_{s,k-1} - X^{(P)}_{s,k}\right)
```

``` math
m_s = X^{(P)}_{s,n_P}\,e^{-\mu_s(t)\tau_M}
```

``` math
\dot{E}_{M,s} = S_{M,s}\Lambda^M_s - m_s - \mu_s(t)E_{M,s},
\qquad
\dot{I}_{M,s} = m_s - \mu_s(t)I_{M,s}.
```

**Vector control** enters only through $`a_s(t)`$ and $`\mu_s(t)`$. For
each species and time, `blink` reuses `malariasimulation`’s own net /
IRS decay and combination algebra, population-averaged over the
net-using and sprayed fractions, to obtain a mean repellency $`\bar{r}`$
and a mean survival $`\bar{p}`$, then

``` math
W = (1-Q_0) + Q_0\,\bar{p},
\qquad
Z = Q_0\,\bar{r},
\qquad
f = \frac{1}{\dfrac{t_f}{1-Z} + \left(\dfrac{1}{f_0} - t_f\right)},
```

``` math
a_s = \left(1 - \frac{1-Q_0}{W}\right) f,
\qquad
\mu_s = -f\log\!\left(\frac{p_1^{(0)}W}{1 - Zp_1^{(0)}}\;e^{-\mu^{(0)}(1/f_0 - t_f)}\right),
\qquad p_1^{(0)} = e^{-\mu^{(0)}t_f}.
```

Oviposition is unaffected: the IBM’s `eggs_laid(beta, mu, f)` is
algebraically identical to $`\beta`$ for all $`\mu, f`$, so $`\beta_s`$
is constant and vector control acts purely through $`\mu`$ and $`a`$.

$`a_s(t)`$ and $`\mu_s(t)`$ are **linearly** interpolated (not step),
because the IBM recomputes them every timestep, so IRS logistic decay
within a spray round is a ramp rather than a step.

### E. Discrete events: chemoprevention

MDA, SMC and PMC are applied as instantaneous state jumps between ODE
segments. The integration is split at every scheduled pulse day; at each
pulse the target age band $`[\ell, u)`$ is mapped onto the model age
grid by **fractional overlap**

``` math
o_i = \frac{\max\!\left(0,\ \min(u, a^+_i) - \max(\ell, a^-_i)\right)}{a^+_i - a^-_i},
```

so narrow bands (e.g. PMC dose ages) are captured rather than dropped,
and coarse groups are not over-treated. A fraction

``` math
\text{frac}_i = o_i \times \text{coverage} \times \text{drug efficacy} \times (1 - \text{ETF})
```

of $`S, U, A, D, T`$ in each affected cell is cleared and moved to
$`P_c`$, which then decays to $`S`$ at $`r_{P_c}`$ — the
coverage-weighted Weibull mean of the chemoprevention drugs, distinct
from the clinical $`r_P`$. Pulses take effect the day **after** their
scheduled timestep. PMC, which the IBM delivers on an age trigger, is
approximated as ~monthly pulses over each dose-age band.

### F. Initial conditions

1.  The `malariasimulation` list is back-translated to
    `malariaEquilibrium` parameter names (`R/translate_params.R`), then
    `parameters$eq_params` is merged over the top and takes precedence.
    Note that `set_equilibrium()` *always* stores its own
    back-translation there, even when called with `eq_params = NULL`, so
    after any `set_equilibrium()` call it is that stored set — not
    `R/translate_params.R` — supplying every shared constant.
2.  For each heterogeneity node $`j`$,
    [`malariaEquilibrium::human_equilibrium_no_het()`](https://rdrr.io/pkg/malariaEquilibrium/man/human_equilibrium_no_het.html)
    is solved at $`\text{EIR} = \text{init_EIR}\times\zeta_j`$ and the
    effective treated fraction $`f_t\cdot\text{eff}`$, on the model age
    grid. Its $`\Lambda`$, $`\phi`$, `prop`, $`r`$, `IB`, `ICA`, `ID`,
    `IVA` and `cA` columns are reused directly. The $`S/T/D/A/U/P`$
    block is **re-solved** with a corrected prophylaxis aging recursion,
    $`b_P = (r_T b_T + r_{i-1}P_{i-1})/\beta_P`$ — upstream divides only
    the $`r_{i-1}P_{i-1}`$ term by $`\beta_P`$, leaving $`r_T b_T`$
    undivided. The re-solve makes the seed the fixed point of the
    prophylaxis ODE actually implemented here, so the model holds flat
    even under treatment.
3.  The equilibrium age structure is recomputed under $`\mu_i(0)`$ —
    which may be the custom-demography baseline row — and the
    constant-hazard seed is rescaled onto it.
4.  $`T_0`$ is split between $`T`$ and $`T_s`$ by the $`t = 0`$
    slow-clearance fraction.
5.  Mosquito compartments are seeded from `initial_mosquito_counts()` at
    the FOIM implied by the seeded human infectivity, with `total_M`
    chosen so that $`\sum_s a_s I_{M,s} = \text{init_EIR}/365`$ exactly.
    Lag chains are seeded at their equilibrium values, and $`K_s(0)`$ at
    the corresponding carrying capacity (floored at $`10^{-9}`$ of the
    maximum for species with proportion 0, which would otherwise give
    $`0/0`$ in the larval term).

Because `malariaEquilibrium` encodes the *simplified* forms — a linear
force of infection, refractory boosting on the raw rate with the
un-rounded $`u`$, no at-risk restriction on boosting, and
continuous-time state exit rates — while `blink` reproduces the IBM’s
per-day semantics, the seed is a close **approximation** to `blink`’s
true fixed point rather than the fixed point itself: an undisturbed run
relaxes by up to ~0.5% over the first years and then holds.
`malariasimulation` relaxes off its own `malariaEquilibrium` seed for
the same reason — though the two seeds are not identical: the IBM passes
the *raw* treated fraction $`f_t`$ rather than $`f_t\cdot\text{eff}`$,
solves on a fixed 0–99.9 y grid with heterogeneity handled inside
`human_equilibrium()`, and does not correct the prophylaxis recursion.
Set `parameters$bite_dedup = 0` to reduce the drift to well under 1% (no
setting makes it exactly zero — see §B.2).

### G. Outputs

Per age group (aggregated to output bands in R), all as fractions of the
total human population and multiplied by `human_population` on the way
out — except `p_detect_lm_*`, which is a proportion (the LM-positive
count divided by the band population) and is not scaled:

| Output | Definition |
|----|----|
| `n_age_*` | $`\sum_j N_{ij}`$ |
| `n_detect_lm_*` | $`\sum_j \left(D + T + T_s + q_{ij}A\right)`$ |
| `p_detect_lm_*` | `n_detect_lm_*` / `n_age_*` for the same band — LM prevalence |
| `n_detect_pcr_*` | $`\sum_j \left(D + T + T_s + A + U\right)`$ — the IBM convention, **not** `malariaEquilibrium`’s sub-patent-weighted `pos_PCR` |
| `n_inc_clinical_*` | $`\sum_j h^c_{ij}\Sigma_{ij}`$ — counted with the *clinical* hazard, so a day’s integral is exactly the IBM’s $`\phi p N`$ |
| `n_inc_severe_*` | $`\sum_j \theta_{ij}\Lambda_{ij}\Sigma_{ij}`$ — severe is drawn from the same infected set as clinical in the IBM, so it keeps the total infection hazard |
| `n_inc_*` | $`\sum_j \Lambda_{ij}\Sigma_{ij}`$ (all new infections) |
| `EIR` | $`365\,\bar\varepsilon`$, per adult per year |
| `FOIM`, `ft` | $`\Lambda^M_1`$ and $`f_t(t)`$ |
| `S_count`, `D_count`, `A_count`, `U_count`, `Tr_count`, `Ph_count` | Population totals by infection state (diagnostic) |

**Six state columns for eight compartments.** `Tr_count` sums
$`T + T_s`$ and `Ph_count` sums $`P + P_c`$, so neither split is
observable in the output. That also makes
`S_count + D_count + A_count + U_count + Tr_count + Ph_count` a complete
population-conservation check.

**Why each family gets its own bands.** Each output family is emitted
**only** over its own `*_rendering_ages` band list, exactly as
`malariasimulation` does, with `n_age_*` over the union. `postie` treats
every output column as an independent age stratum, so if a family were
emitted over the union of all bands, the same ages would appear in two
columns and any person-day-weighted aggregate would count them twice.
This was a real bug (see `NEWS.md`): it inflated all-age clinical
incidence by ~1.21×.

### H. Numerical solution

The model is compiled from odin2 to C++ (`src/malaria_ode.cpp`) and
integrated by dust2’s adaptive solver. Defaults `atol = rtol = 1e-8` and
`step_size_max = 1` preserve the flat equilibrium exactly. The step cap
binds only near equilibrium: there the solution is nearly flat, so the
adaptive stepper proposes a very long step and can jump clean over an
interpolation knot — the day an intervention changes — missing it
entirely. For long dynamic projections, `atol = rtol = 1e-6` with
`step_size_max = 10` runs several times faster with negligible effect on
aggregate outputs.

The PEV, TBV, carrying-capacity and (constant-demography) mortality
grids are built to extend a year past `timesteps` so the stepper never
extrapolates them; the vector-control grid ends exactly at `timesteps`.
Every intervention series starts at its **baseline** value at $`t = 0`$,
so the equilibrium seed is preserved and interventions act only from
their scheduled timesteps.
