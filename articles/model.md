# The fleet model: formal specification

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

The canonical description of the system of ordinary differential
equations `fleet` integrates: §1 for the one design constraint
everything follows from, §2 for what the state space does and does not
carry, and the appendix for the equations as implemented. Read it
alongside `inst/odin/malaria_ode.R` and the input builders in
`R/build_inputs.R` / `R/interventions.R`. Everything below describes
`malariasimulation` **v3.0.0**.

For a task-oriented introduction with runnable code, see [Getting
started with
`fleet`](https://pwinskill.github.io/fleet/articles/fleet.md); for every
`malariasimulation` argument one by one,
[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md);
for what to do differently because this is an ODE and not the IBM,
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md).

## 1. Scope and design

`fleet` is a deterministic mean-field (ODE) counterpart to the
`malariasimulation` individual-based model (IBM) of *Plasmodium
falciparum*.

*Mean field* means the model tracks the **distribution of the population
over states** rather than the people themselves: **the state must not
carry per-individual attributes**. Any quantity that depends on two
attributes belonging to the *same* person — owning a net *and* being
vaccinated, say — is therefore unavailable by construction. Almost
everything in §2.2 and §2.3 follows from that one fact.

`fleet` targets the IBM’s *code*, not the published equations.
`malariasimulation`’s implementation departs from the analytic
Griffin-style model in four places: it deduplicates bites within a day,
rounds immunity-boosting refractory windows to whole days, holds people
in a state for whole days, and applies extrinsic-incubation survival at
a fixed delay. Wherever the two differ, `fleet` reproduces what the IBM
does.

The IBM tracks, for each person, their net, their spray status, their
vaccine doses and antibody titres, their drug and the date they took it,
and the date each of their four immunity functions last boosted. A
mean-field model cannot afford a state dimension for each. `fleet`
therefore keeps only two human structuring dimensions (age and biting
heterogeneity) and represents everything else either as a **time-varying
coefficient** in the ODE, or as a **discrete state jump** between
integration segments. §2.1–2.3 set out exactly which is which.

### Abbreviations

|  |  |
|----|----|
| **EIR** | Entomological inoculation rate: infectious bites per person per year. `fleet` reports it per *adult* per year |
| **FOI** / **FOIM** | Force of infection (on humans) / force of infection on mosquitoes |
| **EIP** | Extrinsic incubation period, the delay in the mosquito between an infectious blood meal and becoming infectious |
| **LM** / **PCR** | Light microscopy / polymerase chain reaction (the two prevalence diagnostics reported) |
| **ETF** / **SPC** | Early treatment failure / slow parasite clearance, the two antimalarial-resistance arms `malariasimulation` implements |
| **PEV** / **TBV** | Pre-erythrocytic vaccine (RTS,S, R21) / transmission-blocking vaccine |
| **TRA** / **TBA** | Transmission-reducing activity (an antibody-level quantity) / transmission-blocking activity (the resulting per-state reduction in infectivity) |
| **MDA / SMC / PMC** | Mass drug administration / seasonal malaria chemoprevention / perennial malaria chemoprevention |
| **EPI** | The Expanded Programme on Immunization routine schedule. Unrelated to `set_epi_outputs()`, where “epi” means *epidemiological* |
| **IRS** | Indoor residual spraying |
| **Erlang chain** | A linear chain of $`n`$ sequential exponential stages. Its total has mean $`\tau`$ and variance $`\tau^2/n`$, so larger $`n`$ sharpens it toward the IBM’s fixed delay |
| **Gauss–Hermite quadrature** | A deterministic stand-in for Monte Carlo: it replaces random per-person draws with a small fixed set of representative values and weights that reproduce the distribution’s moments. “Tensor” means the product grid over several such axes |
| **knots** | The times at which an interpolated series is allowed to change value |
| **seed** | Throughout this article, the *initial state vector* — never a random-number seed. `fleet` has no random component |

## 2. State space

### 2.1 What has a state dimension

Every human compartment is a two-dimensional array over age group $`i`$
and biting-heterogeneity node $`j`$; every mosquito compartment is
indexed by species $`s`$. Human compartments are **per-capita
densities** (fractions of the total human population), which is why run
time is independent of `human_population`.

| Dimension | Symbol | Default size | Set by | What it indexes |
|----|----|----|----|----|
| Age group | $`i = 1,\dots,n_a`$ | 52 | `age_lower =` / [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md) | Graded age grid: monthly to 1 y, quarterly to 5 y, yearly to 15 y, 5-yearly to 80 y, plus an absorbing top group. Individuals age by a linear-chain (McKendrick) flow $`r_i`$ between adjacent groups. |
| Biting heterogeneity | $`j = 1,\dots,n_z`$ | 5 | `parameters$n_heterogeneity_groups`, `enable_heterogeneity` | Gauss–Hermite nodes $`\zeta_j`$ and weights $`w_j`$ of the log-normal relative-biting-rate distribution with variance $`\sigma^2`$ (`sigma_squared`); $`\zeta_j = \exp(z_j\sigma - \sigma^2/2)`$. Fixed at birth, never mixed: a person stays in their stratum for life, as in the IBM. |
| Vector species | $`s = 1,\dots,n_v`$ | 1 | `set_species()` | Independent aquatic + adult mosquito populations, each with its own biting rate, mortality, carrying capacity and vector-control response, coupled only through the shared human population. |
| EIR lag stage | $`1,\dots,n_E`$ | 10 (`n_eir`) | `ode_tuning(n_eir =)` | Erlang chain approximating the sporozoite-to-infection delay $`\tau_E`$ (`de`). |
| FOIM lag stage | $`1,\dots,n_F`$ | 10 (`n_foim`) | `ode_tuning(n_foim =)` | Erlang chain approximating the gametocyte delay $`\tau_l`$ (`delay_gam`). |
| EIP stage | $`1,\dots,n_P`$ | 20 (`n_eip`) | `ode_tuning(n_eip =)` | Loss-free Erlang chain approximating the fixed extrinsic incubation period $`\tau_M`$ (`dem`). |

Human compartments, all $`[n_a \times n_z]`$. The last column gives the
column name in the output table, which is where you will meet these in
practice. Eight compartments are reported in six columns.

| Compartment | odin | Output column | Meaning |
|----|----|----|----|
| $`S`$ | `S` | `S_count` | Susceptible |
| $`D`$ | `D` | `D_count` | Clinical disease, untreated |
| $`T`$ | `Tr` | `Tr_count` (with $`T_s`$) | Successfully treated, standard clearance |
| $`T_s`$ | `Tr_slow` | folded into `Tr_count` | Successfully treated, **slow parasite clearance**. A separate compartment, so the two clearance speeds form a genuine mixture of exponentials rather than one exponential at the blended mean (§B.4) |
| $`A`$ | `A` | `A_count` | Asymptomatic patent infection |
| $`U`$ | `U` | `U_count` | Sub-patent infection |
| $`P_{1..k_P}`$ | `Ph` | `Ph_count` (with $`P_c`$) | Post-treatment prophylaxis, an **Erlang chain** of $`k_P`$ stages so that the exit-time distribution of the whole $`T + P`$ sojourn approximates the IBM’s Weibull protection rather than an exponential at its mean (§B.4): $`k_P`$ = 16 for SP-AQ, 20 for DHA-PQP, 1 for AL (whose 10-day protection is less variable than $`T`$ itself); 1 with no clinical drugs |
| $`P_{c,1..k_{P_c}}`$ | `Ph_c` | folded into `Ph_count` | Chemoprevention prophylaxis (MDA/SMC/PMC), its own chain because its drug differs from the clinical first line, and hence so do its mean duration and shape; $`k_{P_c} = \text{round}(1/\text{CV}_W^2)`$: 14 for SP-AQ, 15 for DHA-PQP |
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
| $`E`$, $`L`$, $`P_L`$ | Aquatic early larval, late larval, pupal (odin names `ME`, `ML`, `MP`). Pupae emerge at $`P_L/d_{P_L}`$, of which the **half that are female** enter $`S_M`$; males are not modelled, so the factor of $`\tfrac12`$ on emergence is a sex ratio and not a loss |
| $`S_M`$ | Susceptible adult females |
| $`X^{(P)}_{s,1..n_P}`$ | EIP delay chain (odin `Xi`): a pure delay line, carries no mortality |
| $`E_M`$ | Incubating (exposed) adults (odin `Em_inc`) |
| $`I_M`$ | Infectious adults |

**Total state size.** With the defaults ($`n_a = 52`$, $`n_z = 5`$,
$`n_v = 1`$, $`n_E = n_F = 10`$, $`n_P = 20`$) and no drugs, so that
$`k_P = k_{P_c} = 1`$:

``` math
\underbrace{(10 + k_P + k_{P_c})\,n_a n_z}_{3120\ \text{human}} \;+\; \underbrace{n_E + n_F}_{20\ \text{lags}} \;+\; \underbrace{n_v(6 + n_P)}_{26\ \text{mosquito}} \;=\; 3166 .
```

Every extra prophylaxis stage adds $`n_a n_z = 260`$ states: an SP-AQ
SMC scenario ($`k_{P_c} = 14`$) has 6546.

### 2.2 What is captured *without* a dimension

These are the mechanisms the IBM carries as per-individual attributes
and `fleet` carries as **time-varying (and sometimes age-varying)
coefficients**. Each row is a dimension that was *not* added.

| Mechanism | IBM representation | `fleet` representation | Cost of the collapse |
|----|----|----|----|
| **Clinical treatment coverage** | Per-person Bernoulli draw at the moment of a clinical episode | Scalar $`f_t(t)`$, step-interpolated over the union of all drugs’ `timesteps`; splits the clinical inflow between $`T/T_s`$ and $`D`$ | None at the flow level |
| **Which drug a person received** | Drug index stored per treated individual, driving their efficacy, infectivity and prophylaxis duration | Three scalar series $`\text{eff}(t)`$, $`c_T(t)`$, $`r_P(t)`$, each a coverage-share-weighted blend across the active clinical drugs, step-interpolated over the coverage change times | A first-line **switch** is modelled (the blend moves in time); a genuinely *mixed* first line is represented by its mean, not its mixture |
| **Post-treatment prophylaxis clock** | Weibull survival $`W(t - t_{\text{drug}})`$ multiplying each treated person’s infection probability, running from the dose | An Erlang chain entered after $`T`$, its mean and its stage count moment-matched to $`W`$ (§B.4) | The chain reproduces the integrated protection exactly and the survival curve to within a few points (§B.4) |
| **Early treatment failure** | Per-person draw diverting a treated case back to clinical | Folded into the effective treated fraction $`f_t^{\text{eff}}(t)`$ (§B.4) | None |
| **Bed net ownership, net age, retention, repeat rounds** | Per-person net-receipt timestep; efficacy decays from that date; nets lost stochastically | Two scalars per species, $`a_s(t)`$ and $`\mu_s(t)`$, computed from the full mean-field mixture over *all past distributions* (each weighted by most-recent-receipt probability $`\times`$ retention survival, each decaying its own $`d_{n0}`$ / $`r_n`$) | Protection is population-averaged rather than correlated within individuals across bites, so nets and IRS tend to **under**-suppress transmission relative to the IBM |
| **IRS spray status and spray age** | Per-person house-spray timestep | Folded into the same $`a_s(t)`$, $`\mu_s(t)`$; mixture over all past rounds weighted by most-recent-spray probability, with no retention factor (sprayed protection never expires in the IBM) | As above |
| **Vaccination status, dose count, booster stratum, antibody titre** | Per-person dose history and last-vaccination date. Antibody parameters are *not* stored: the IBM re-draws all four from their profile distributions every timestep, for every vaccinated person | One $`[n_a \times n_t]`$ multiplier $`v_i(t)`$ on the infection hazard. The vaccinated are partitioned analytically by the most recent booster each has reached; each stratum’s efficacy is the **average of the Hill efficacy over the IBM’s own 4-D antibody distribution** ([`create_pev_profile()`](https://pwinskill.github.io/fleet/articles/parameters.html#create_pev_profile)) | The antibody average is exact, since the IBM’s draws are independent. What is lost is the dose *history*: vaccinated status is uncorrelated with anything else (nets, treatment), `min_wait` re-vaccination exclusion is not applied ([`set_mass_pev()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_mass_pev)), and seasonal boosters become a fixed days-since-primary schedule (warned) |
| **TBV antibody titre** | Per-person, mapped to transmission-reducing activity | Four $`[n_a \times n_t]`$ infectivity multipliers, one per infectious state ($`U`$, $`A`$, $`D`$, $`T`$), because the IBM’s TRA-to-TBA map is state-specific (§C) | Mean-field only |
| **Immunity boosting refractory clock** | Per-person `last_boosted` timestep; a boost fires only if $`\ge u`$ days have passed | An analytic renewal rate on the *deduplicated per-day* event probability, with the wait rounded as the IBM rounds it (§B.5) | Exact for the mean inter-boost interval |
| **Age-specific mortality / demographic transition** | Per-person death hazard by age band and year | $`\mu_i(t)`$, a constant-interpolated $`[n_a \times n_t]`$ series; its $`t = 0`$ row also sets the seed’s age structure (§F) | Ages above the top `deathrate_agegroups` take the top rate (the IBM removes them); a warning fires |
| **Maternal immunity** | Sampled from a mother drawn at birth | Purely algebraic, from the mean acquired immunity of 20-year-old mothers (§B.3) — **no state variable at all** | None beyond the mean-field average over mothers |
| **Seasonality** | Daily rainfall drives larval carrying capacity | $`K_s(t)`$, linearly interpolated on a daily grid, from the same truncated Fourier series | None (identical functional form) |
| **Repeated bites in one day** | Bitten individuals are collected in a set, so a person bitten several times counts once and can be infected at most once per day | Saturating hazard instead of the linear $`b\varepsilon`$ (§B.2) | None at the daily-flow level; it is the *linear* form that would be wrong |

Mass drug administration is carried without a dimension too, but not as
a coefficient: it cannot be written as a smooth rate. It is a short,
near-instantaneous clearance of a target age band, so `fleet` splits the
integration at every scheduled MDA/SMC/PMC timestep, applies a **state
jump** (§E), then resumes, keeping chemoprevention out of the ODE
right-hand side entirely.

### 2.3 What is not captured at all

| Feature | Why | Behaviour |
|----|----|----|
| *P. vivax* | A different model (hypnozoites, relapse) | Error at input |
| Individual mosquitoes (`individual_mosquitoes = TRUE`) | `fleet` is always compartmental | Ignored |
| Intervention correlation (`get_correlation_parameters()`, `run_simulation(correlations =)`) | Purely individual-level; the mean field assumes independence between interventions | Warned and ignored |
| Per-individual stochastic variation | Deterministic model | By construction |
| Partner-drug resistance, late clinical / parasitological failure, reinfection during prophylaxis | `malariasimulation` itself rejects non-zero values for these upstream in `set_antimalarial_resistance()` | Cannot reach the model |

## A. Notation

Abbreviations are at the end of §1.

### Symbols

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

## B. The human model

### B.1 Demography

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

### B.2 Exposure and the infection hazard

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
at most once. `fleet` reproduces that cap; the daily infection
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
and diverge at seasonal peaks and in high-$`\zeta`$ strata.
`bite_dedup = 0` does **not** make the seed an exact fixed point: it
gates only the hazard, while the boosting terms of §B.5 keep the
deduplicated event probability, the integer refractory window
$`\lceil u\rceil - 1`$ and the at-risk restriction, and the state exit
rates keep the whole-day form $`1 - e^{-1/d}`$. What it does is roughly
halve the residual drift. Measured as the largest excursion from the
seeded value over an undisturbed 15-year run at EIR 20, PfPR(2-10)
departs by at most 0.14% with `bite_dedup = 0` against 0.28% at the
default. The *endpoint* tells the opposite story (the default happens to
return close to its seed, ending -0.05% off it against +0.14% with `0`),
so the metric has to be stated.

**Hazard splitting.** In the IBM each person resolves at most one
infection outcome per day, so the daily clinical count is $`\phi p N`$
for susceptibles. (For people already in $`A`$ or $`U`$ the infection
outcome competes with disease progression inside the *same* hazard
resolution, which scales their realised infection probability by
$`\frac{\Lambda}{\Lambda + r}\cdot\frac{1 - e^{-(\Lambda+r)}}{1 - e^{-\Lambda}}`$.
`fleet`’s ODE reproduces that competition automatically, since
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

### B.3 Immunity-dependent probabilities

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

$`q`$ takes **no** offset, matching the IBM.

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

### Why severe incidence comes out low

Of the three acquired-immunity Hill functions, $`\theta`$ has the
*shallowest* exponent, not the steepest:

|            | exponent        | scale              |
|------------|-----------------|--------------------|
| $`b`$      | $`k_b = 2.160`$ | $`I_{B0} = 43.9`$  |
| $`\phi`$   | $`k_c = 2.369`$ | $`I_{C0} = 18.02`$ |
| $`\theta`$ | $`k_v = 2.000`$ | $`I_{V0} = 1.096`$ |

What singles $`\theta`$ out is the **scale**, not the exponent. A Hill
function is flat well below its scale parameter and only approaches its
full log-log slope well above it, so the local slope is a matter of
position. Writing $`X =
f_{v,i}\left(\left(I_{VA}+\delta+I_{VM}\right)/I_{V0}\right)^{k_v}`$ for
the denominator term,

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

Position separates them. $`I_{V0}`$ is sixteen times smaller than
$`I_{C0}`$ while the two acquired immunities take nearly the same values
at every age, so they differ almost entirely in what they are divided
by. At `fleet`’s own seed at the reference EIR of 20:

| age | $`I_{VA}/I_{V0}`$ | $`\partial\log\theta/\partial\log I_{VA}`$ | $`I_{CA}/I_{C0}`$ | $`\partial\log\phi/\partial\log I_{CA}`$ |
|----|----|----|----|----|
| 1 | 1.80 | $`-0.65`$ | 0.12 | $`-0.02`$ |
| 3 | 6.62 | $`-1.80`$ | 0.44 | $`-0.30`$ |
| 5 | 12.20 | $`-1.95`$ | 0.82 | $`-0.91`$ |
| 10 | 26.31 | $`-1.87`$ | 1.79 | $`-1.88`$ |

$`\theta`$ is within 5% of the steepest slope it can have from age 3.4
to age 10.5. $`\phi`$ at age 3 is less than a fifth as steep, and does
not overtake $`\theta`$ until age 10.

Sitting out on that tail makes $`\theta`$**convex** in $`I_{VA}`$ across
essentially the whole population, where $`\phi`$ is concave in
$`I_{CA}`$ through early childhood at this EIR. That fixes the sign of
the mean-field error. `fleet` evaluates $`\theta`$ once per (age group
$`\times`$ heterogeneity group) at that stratum’s mean $`I_{VA}`$; the
IBM averages $`\theta`$ over individuals whose bite histories differ
*inside* the stratum. By Jensen’s inequality
$`\mathbb{E}[\theta(I_{VA})] \geq
\theta(\mathbb{E}[I_{VA}])`$, so the mean field returns the smaller
number and severe incidence comes out low. `fleet` resolves the
log-normal biting strata explicitly (§2.1), so only the spread within an
(age, $`\zeta`$) cell is lost and the realised gap is a few percent
rather than tens of percent, but the direction is the same.
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
says what to do about it.

### B.4 Disease-state equations

Writing $`\Sigma_{ij} = S_{ij}+A_{ij}+U_{ij}`$ for the at-risk pool,
$`\text{SPC}(t)`$ for the slow-parasite-clearance fraction of the
treated
([`set_antimalarial_resistance()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_antimalarial_resistance)),
and $`f_t^{\text{eff}}(t) = f_t(t)\,\text{eff}(t)\,(1-\text{ETF}(t))`$
for the effective treated fraction:

``` math
\dot{S}_{ij} = r_{i-1}S_{i-1,j} + r_U U_{ij} + k_P r_P P_{k_P,ij} + k_{P_c} r_{P_c} P_{c,k_{P_c},ij}
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
\dot{P}_{1,ij} = r_{i-1}P_{1,i-1,j} + r_T T_{ij} + r_T^{\text{slow}} T_{s,ij}
- k_P r_P P_{1,ij} - \rho_i P_{1,ij}
```

``` math
\dot{P}_{m,ij} = r_{i-1}P_{m,i-1,j} + k_P r_P\left(P_{m-1,ij} - P_{m,ij}\right) - \rho_i P_{m,ij},
\qquad m = 2,\dots,k_P
```

``` math
\dot{P}_{c,1,ij} = r_{i-1}P_{c,1,i-1,j} - k_{P_c} r_{P_c} P_{c,1,ij} - \rho_i P_{c,1,ij}
```

``` math
\dot{P}_{c,m,ij} = r_{i-1}P_{c,m,i-1,j} + k_{P_c} r_{P_c}\left(P_{c,m-1,ij} - P_{c,m,ij}\right) - \rho_i P_{c,m,ij},
\qquad m = 2,\dots,k_{P_c}
```

Four implementation points:

- **Prophylaxis chains.** The IBM does not have a prophylaxis *state*:
  it multiplies each treated person’s infection probability by the
  Weibull survival $`W(t - t_{\text{drug}})`$, so a treated cohort’s
  mean protection at lag $`t`$ is exactly $`W(t)`$. A single compartment
  reproduces only the mean of $`W`$ and decays exponentially, keeping
  just 42% of SP-AQ recipients protected at day 30 against the Weibull’s
  70%. The $`k`$-stage chains above, each stage left at $`k\,r`$, keep
  the mean $`1/r`$ and are moment-matched to $`W`$. For $`P_c`$ that is
  simply $`k_{P_c} = \text{round}(1/\text{CV}_W^2)`$ (14 for SP-AQ),
  which puts the chain’s survival within a few points of $`W`$ (0.67 at
  day 30). $`P`$ is entered only after the exponential $`T`$ sojourn,
  whereas the IBM’s clock runs from the dose in parallel with $`T`$: a
  treated person there is unprotected at lag $`t`$ with probability
  $`F_T(t)\,(1 - W(t))`$. Matching the integrated protection gives the
  chain mean $`\bar{d}_W - \int_0^\infty e^{-r_T t} W(t)\,dt`$ (which
  reduces to $`\bar{d}_W - 1/r_T`$ when $`W \approx 1`$ throughout
  $`T`$), and matching the variance of the whole $`T + P`$ sojourn gives
  $`k_P`$: 16 for SP-AQ, 20 for DHA-PQP, and 1 for AL, whose 10-day
  protection is already less variable than $`T`$ itself, so no chain
  length can sharpen it (its output is the same to five decimals at 1
  and 20 stages). While in a chain, people are fully protected from
  $`\Lambda`$ and are not infectious.

- **Two treated compartments.** The IBM assigns each successfully
  treated individual to standard or slow parasite clearance by a
  Bernoulli draw, giving a *mixture of two exponentials*, not one
  exponential at the blended mean. $`T`$ and $`T_s`$ reproduce that
  split structurally (but see the sojourn caveat below for
  $`r_T^{\text{slow}}`$).

- **Whole-day sojourns.** The IBM advances states once per day with exit
  probability $`1 - e^{-1/d}`$, so its realised mean dwell is
  $`1/(1-e^{-1/d})`$ (5.52 d for $`d_D = 5`$, not 5). `fleet` uses
  $`r_X = 1 - e^{-1/d_X}`$ for $`X \in \{A, D, U, T, T_{\text{slow}}\}`$
  so the dwells match (slow clearance is resolved through the same
  per-day `rate_to_prob` path as ordinary treatment in the IBM). It is
  *not* applied to $`r_P`$ or $`r_{P_c}`$, because the IBM models
  prophylaxis as a hazard multiplier rather than a compartment with a
  daily census.

- **The $`P_c`$ chain has no ODE inflow.** It is filled, at stage 1,
  only by the chemoprevention pulses of §E, and carries neither
  infection risk nor infectivity.

### B.5 Immunity equations

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

## C. Human → mosquito infectivity

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

## D. The mosquito model

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
$`e^{-\mu\tau_M}`$ at the exit, with the current $`\mu`$. `fleet`
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
each species and time, `fleet` reuses `malariasimulation`’s own net /
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

## E. Discrete events: chemoprevention

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

of **every** state in each affected cell — $`S, U, A, D, T, T_s`$ and
the existing $`P`$ and $`P_c`$ stages, because the IBM resets
`drug_time` for everyone it successfully treats, renewing the protection
of those already covered — is cleared and moved to the first stage of
$`P_c`$, which then runs down its chain to $`S`$ at $`k_{P_c} r_{P_c}`$
per stage: $`r_{P_c}`$ is 1 / the coverage-weighted Weibull mean of the
chemoprevention drugs and $`k_{P_c}`$ their coverage-weighted
$`1/\text{CV}_W^2`$, both distinct from the clinical chain. Pulses take
effect the day **after** their scheduled timestep. PMC, which the IBM
delivers on an age trigger, is approximated as ~monthly pulses over each
dose-age band.

## F. Initial conditions

1.  The `malariasimulation` list is back-translated to
    `malariaEquilibrium` parameter names (`R/translate_params.R`), then
    `parameters$eq_params` is merged over the top and takes precedence.
    `set_equilibrium()` *always* stores its own back-translation there,
    even when called with `eq_params = NULL`, so after any
    `set_equilibrium()` call it is that stored set, not
    `R/translate_params.R`, supplying every shared constant.
2.  For each heterogeneity node $`j`$,
    [`malariaEquilibrium::human_equilibrium_no_het()`](https://rdrr.io/pkg/malariaEquilibrium/man/human_equilibrium_no_het.html)
    is solved at $`\text{EIR} = \text{init_EIR}\times\zeta_j`$ and the
    effective treated fraction $`f_t\cdot\text{eff}`$, on the model age
    grid. Its $`\Lambda`$, $`\phi`$, `prop`, $`r`$, `IB`, `ICA`, `ID`,
    `IVA` and `cA` columns are reused directly. The $`S/T/D/A/U/P`$
    block is **re-solved** with a corrected prophylaxis aging recursion,
    $`b_P = (r_T b_T + r_{i-1}P_{i-1})/\beta_P`$ — upstream divides only
    the $`r_{i-1}P_{i-1}`$ term by $`\beta_P`$, leaving $`r_T b_T`$
    undivided — and generalised to the $`k_P`$-stage chain (stage 1 fed
    by $`r_T T`$, stage $`m`$ by $`k_P r_P P_{m-1}`$, every stage aging
    at $`r_{i-1}`$). The re-solve makes the seed the fixed point of the
    prophylaxis ODE actually implemented here, so the model holds flat
    even under treatment.
3.  The equilibrium age structure is recomputed under $`\mu_i(0)`$
    (which may be the custom-demography baseline row), and the
    constant-hazard seed is rescaled onto it.
4.  $`T_0`$ is split between $`T`$ and $`T_s`$ by the $`t = 0`$
    slow-clearance fraction.
5.  Mosquito compartments are seeded from `initial_mosquito_counts()` at
    the FOIM implied by the seeded human infectivity, with `total_M`
    chosen so that $`\sum_s a_s I_{M,s} = \text{EIR}_{\text{seed}}/365`$
    exactly. Under the default demography
    $`\text{EIR}_{\text{seed}} = \text{init_EIR}`$; under a custom
    demography `total_M` is fixed to the IBM’s (`ibm_total_M()`,
    [`set_equilibrium()`](https://pwinskill.github.io/fleet/articles/parameters.html#set_equilibrium))
    and $`\text{EIR}_{\text{seed}}`$ is the root of
    `total_M(EIR) = total_M_IBM`, found by
    [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) on
    $`\log \text{EIR}`$ (steps 2–4 are re-solved at each trial EIR). Lag
    chains are seeded at their equilibrium values, and $`K_s(0)`$ at the
    corresponding carrying capacity (floored at $`10^{-9}`$ of the
    maximum for species with proportion 0, which would otherwise give
    $`0/0`$ in the larval term).

Because `malariaEquilibrium` encodes the *simplified* forms — a linear
force of infection, refractory boosting on the raw rate with the
un-rounded $`u`$, no at-risk restriction on boosting, and
continuous-time state exit rates — while `fleet` reproduces the IBM’s
per-day semantics, the seed is a close **approximation** to `fleet`’s
true fixed point rather than the fixed point itself: an undisturbed run
relaxes by up to ~0.5% over the first years and then holds.
`malariasimulation` relaxes off its own `malariaEquilibrium` seed for
the same reason, though the two seeds are not identical: the IBM passes
the *raw* treated fraction $`f_t`$ rather than $`f_t\cdot\text{eff}`$,
solves on a fixed 0–99.9 y grid with heterogeneity handled inside
`human_equilibrium()`, and does not correct the prophylaxis recursion.
Set `parameters$bite_dedup = 0` to reduce the drift to well under 1% (no
setting makes it exactly zero; see §B.2).

## G. Outputs

Per age group (aggregated to output bands in R), all as fractions of the
total human population and multiplied by `human_population` on the way
out, except `p_detect_lm_*`, which is a proportion (the LM-positive
count divided by the band population) and is not scaled:

| Output | Definition |
|----|----|
| `n_age_*` | $`\sum_j N_{ij}`$ |
| `n_detect_lm_*` | $`\sum_j \left(D + T + T_s + q_{ij}A\right)`$ |
| `p_detect_lm_*` | `n_detect_lm_*` / `n_age_*` for the same band (LM prevalence) |
| `n_detect_pcr_*` | $`\sum_j \left(D + T + T_s + A + U\right)`$: the IBM convention, **not** `malariaEquilibrium`’s sub-patent-weighted `pos_PCR` |
| `n_inc_clinical_*` | $`\sum_j h^c_{ij}\Sigma_{ij}`$, counted with the *clinical* hazard, so a day’s integral is exactly the IBM’s $`\phi p N`$ |
| `n_inc_severe_*` | $`\sum_j \theta_{ij}\left[\Lambda_{ij}\left(S_{ij}+U_{ij}\right) + p_{ij}A_{ij}\right]`$: $`\theta`$ times the all-infection count below, because the IBM draws severe from the same *deduplicated* infected set as clinical |
| `n_inc_*` | $`\sum_j \left[\Lambda_{ij}\left(S_{ij}+U_{ij}\right) + p_{ij}A_{ij}\right]`$ (all new infections) |
| `EIR` | $`365\,\bar\varepsilon`$, per adult per year |
| `FOIM`, `ft` | $`\Lambda^M_1`$ and $`f_t(t)`$ |
| `S_count`, `D_count`, `A_count`, `U_count`, `Tr_count`, `Ph_count` | Population totals by infection state (diagnostic) |

**How age groups map onto bands.** By **exact overlap**: a model age
group contributes to a rendering band the fraction of its own width that
lies inside $`[\ell, u)`$, the same weighting the chemoprevention pulses
use (§E). A group straddling a band boundary is split between the two
bands rather than handed whole to one of them, so over any set of bands
that partition a span of the age axis the weights sum to exactly 1 per
group inside that span: no family double-counts, and none loses anyone
it renders. The absorbing top group is the one case the invariant cannot
cover, since it is not a finite interval. It is captured whole by any
band reaching above its lower edge, and `fleet` reports the population
that no band captured.

**Why the two infection counts split the pool.** These outputs are
*rates*, which R integrates over a day, so the right rate depends on how
fast each compartment drains. $`S`$ and $`U`$ drain at the full
$`\Lambda`$, and
$`\int_0^1 \Lambda X e^{-\Lambda t}\,\mathrm{d}t = X\left(1 - e^{-\Lambda}\right) = Xp`$
— the IBM’s deduplicated count, exactly. $`A`$ drains only at $`h^c`$,
because a sub-clinical re-infection leaves an $`A`$ in $`A`$; counting
it at $`\Lambda`$ would over-count by $`\Lambda/p`$, so $`A`$ takes
$`p_{ij}`$ directly. Clinical incidence is unaffected: it is counted
with $`h^c`$, and
$`\int h^c A\,\mathrm{d}t = A\left(1-e^{-h^c}\right) = A\phi p`$
already.

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

## H. Numerical solution

The model is compiled from odin2 to C++ (`src/malaria_ode.cpp`) and
integrated by dust2’s adaptive solver. Defaults `atol = rtol = 1e-8` and
`step_size_max = 1` preserve the flat equilibrium exactly. The step cap
binds only near equilibrium: there the solution is nearly flat, so the
adaptive stepper proposes a very long step and can jump clean over an
interpolation knot (the day an intervention changes), missing it
entirely. See
[`?ode_tuning`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
for when to loosen these settings and what it costs.

The PEV, TBV, carrying-capacity and (constant-demography) mortality
grids are built to extend a year past `timesteps` so the stepper never
extrapolates them; the vector-control grid ends exactly at `timesteps`.
Every intervention series starts at its **baseline** value at $`t = 0`$,
so the equilibrium seed is preserved and interventions act only from
their scheduled timesteps.

That last property is a property of the **knots**, not of the series’
values. The vector-control, carrying-capacity, PEV and TBV series are
interpolated linearly, which is right for the within-round decay each of
them carries but means a value that changes on a scheduled day ramps in
from whatever knot precedes it. Each of these grids therefore carries a
knot at $`t_k - 1`$ as well as at every scheduled day $`t_k`$, pinning
the pre-deployment value one day out and confining the ramp to the
single step onto the scheduled day. Without it a deployment inherits the
spacing of the background grid: off the 10-day vector-control grid, nets
scheduled for day 100 reach 66% of their effect on EIR by day 99. With
it, a round takes effect within a day of its schedule, the same
resolution the chemoprevention pulses (§E) already have.
