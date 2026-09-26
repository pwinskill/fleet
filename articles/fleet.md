# Getting started with fleet

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

`fleet` is a fast, deterministic **mean-field twin** of the
[malariasimulation](https://github.com/mrc-ide/malariasimulation)
individual-based model of *Plasmodium falciparum* and *P. vivax*
malaria. It reproduces the same age- and biting-heterogeneity-structured
human model and compartmental mosquito model, advances them one day at a
time in the IBM’s own order, takes the **same parameter list** as the
IBM, and is seeded at the
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium)
solution (for vivax, `malariaEquilibriumVivax`’s). A 30-year falciparum
run takes under two seconds and a vivax run under half a minute,
independent of population size. That speed is the point: the intended
uses are calibration, sweeps and quick scenario work, where a stochastic
IBM would be too slow. It is not ready for those jobs yet (see the note
above).

Reach for the IBM (not this) when you need stochastic variation or
individual heterogeneity beyond the mean field.

## Install and load

``` r

# install.packages("remotes")
remotes::install_github("pwinskill/fleet")
```

``` r

library(fleet)
```

`malariasimulation` (for parameter lists and the `set_*` builders) and
`postie` (for post-processing) are used throughout; both are listed
under Suggests.

## A basic run

Build a `malariasimulation` parameter list and run it, exactly as you
would run the IBM. You supply a target adult EIR (infectious bites per
adult per year), and the model is seeded at the corresponding
equilibrium. The list asks for the age bands to report, as it does of
the IBM: here incidence in under-5s and at all ages, beside the default
2-10 prevalence band.

``` r

p <- malariasimulation::get_parameters(list(
  human_population = 1000,
  clinical_incidence_rendering_min_ages = c(0, 0),
  clinical_incidence_rendering_max_ages = c(5, 100) * 365 - 1,
  severe_incidence_rendering_min_ages = c(0, 0),
  severe_incidence_rendering_max_ages = c(5, 100) * 365 - 1))

# 15 years at a daily step, seeded at the equilibrium for EIR = 20
out <- run_simulation_ode(timesteps = 15 * 365, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))

out[1:3, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650",
           "n_inc_clinical_0_1824", "EIR")]
#>   timestep n_age_730_3650 n_detect_lm_730_3650 n_inc_clinical_0_1824 EIR
#> 1        1       288.1847             158.3112             0.8722185  20
#> 2        2       288.1847             158.3001             0.8724063  20
#> 3        3       288.1847             158.2896             0.8725619  20
```

The return value is `malariasimulation`’s daily table, with the columns
the IBM would give this parameter list, under the same names and with
the same meanings. Column tags carry the age band **in days**, both ends
included: `730_3650` is the 2-10y band, and `0_1824` the under-5s. Every
`n_*` column is a count: `n_detect_lm_*` the LM-positive and
`n_inc_clinical_*` the day’s new clinical cases in the band.
(`p_detect_lm_*` is the IBM’s expected-count partner of `n_detect_lm_*`,
not a prevalence.) `EIR` is per adult per year. Because the model is
mean-field, rates are population-independent; `human_population` only
rescales the counts.

With no interventions the run settles rather than holding perfectly
still:

``` r

pfpr <- out$n_detect_lm_730_3650 / out$n_age_730_3650   # LM prevalence, 2-10y
range(out$EIR)   # ~20 throughout
#> [1] 19.92115 20.00000
range(pfpr)      # a ~0.4% relaxation over the first years
#> [1] 0.5473898 0.5493395
```

Row 1 is day 1, whose state is the seed. Over these fifteen years LM
prevalence in 2-10 year olds runs between 0.5495 (the seed) and 0.5475
(its low point, a little over a year in), a spread of 0.35%, and is flat
to 0.008% over the last five years. EIR stays within 0.4% of its target
throughout. Incidence relaxes further off the seed than prevalence does,
a few per cent over the first years at this EIR and more at higher ones.

That relaxation is expected, not a numerical error. The
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium) seed
solves the model’s continuous-time form, so it is a close
**approximation** to `fleet`’s true fixed point on the daily clock
rather than the fixed point itself; `malariasimulation` is seeded from
the same solution and drifts off it too.
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
lists the differences. **Burn in before calibrating, and read a settled
window rather than the first days.**

## Reading outputs with postie

Because the table is `malariasimulation`’s, it goes straight into
[postie](https://github.com/mrc-ide/postie): the same two calls you
would make on an IBM run, with no fleet-specific wrapper in between.

``` r

# prevalence: one <diagnostic>_prevalence_<lo>_<hi> column per age band (years)
prevalence <- postie::get_prevalence(out, diagnostic = "lm")
utils::tail(prevalence["lm_prevalence_2_10"], 3)
#>      lm_prevalence_2_10
#> 5473          0.5478469
#> 5474          0.5478469
#> 5475          0.5478469

# rates: clinical / severe incidence, mortality and DALYs by age band
rates <- postie::get_rates(out)
#> Warning in postie::get_rates(out): required column `ft` not found, assuming ft
#> = 0
#> Warning in postie::get_rates(out): required column `ft_sev` (probability
#> hospitalisation | severe case) not found, assuming ft_sev = 0.8
utils::head(rates[, c("time", "age_lower", "age_upper",
                      "clinical", "severe", "dalys")])
#> # A tibble: 6 × 6
#>    time age_lower age_upper clinical    severe    dalys
#>   <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 2000          0      5.00  0.00412 0.000138  0.00131 
#> 2 2000          0    100.0   0.00184 0.0000327 0.000176
#> 3 2000.         0      5.00  0.00412 0.000138  0.00131 
#> 4 2000.         0    100.0   0.00184 0.0000327 0.000176
#> 5 2000.         0      5.00  0.00412 0.000138  0.00131 
#> 6 2000.         0    100.0   0.00184 0.0000327 0.000176
```

Use `diagnostic = "pcr"` for PCR prevalence, and postie’s own arguments
(`scaler`, `treatment_scaler`, `life_expectancy`, …) as you normally
would. This means a post-processing pipeline written for the IBM works
on a `fleet` run unchanged.

### Choosing age bands

The bands come from the `malariasimulation` rendering fields on the
parameter list, the same fields the IBM reads, and each family is
reported over its own list only: `prevalence_rendering_*` for detection,
`clinical_incidence_*`, `severe_incidence_*` and `incidence_rendering_*`
for the three incidence families, and `age_group_rendering_*` for
population alone. A band need not line up with the model’s age groups:

``` r

p_bands <- p
p_bands$prevalence_rendering_min_ages <- c(0,  2, 10) * 365
p_bands$prevalence_rendering_max_ages <- c(2, 10, 100) * 365 - 1
p_bands$clinical_incidence_rendering_min_ages <- c(0, 2) * 365
p_bands$clinical_incidence_rendering_max_ages <- c(2, 10) * 365 - 1

out_bands <- run_simulation_ode(10 * 365, malariasimulation::set_equilibrium(p_bands, init_EIR = 20))
grep("^n_detect_lm_|^n_inc_clinical_", names(out_bands), value = TRUE)
#> [1] "n_inc_clinical_0_729"    "n_inc_clinical_730_3649"
#> [3] "n_detect_lm_0_729"       "n_detect_lm_730_3649"   
#> [5] "n_detect_lm_3650_36499"
```

## Layering interventions

Interventions are configured with the ordinary
`malariasimulation::set_*` builders and applied automatically from the
list. There are no extra arguments to
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md).
As in the IBM, settings read on the day (treatment coverage, the drug
mix, resistance, death rates) change on their own timestep, and
deployments (nets, spraying, vaccination, chemoprevention) act from the
day after theirs. A few representative examples follow.

### Bed nets

``` r

p_nets <- malariasimulation::set_bednets(
  malariasimulation::get_parameters(list(human_population = 1000)),
  timesteps = 5 * 365, coverages = 0.8, retention = 5 * 365,
  dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
  rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)

nets <- run_simulation_ode(12 * 365, malariasimulation::set_equilibrium(p_nets, init_EIR = 20))
pfpr_nets <- nets$n_detect_lm_730_3650 / nets$n_age_730_3650
c(baseline = pfpr_nets[1], trough = min(pfpr_nets))
#>  baseline    trough 
#> 0.5493395 0.1810560
```

Prevalence falls after deployment and then relaxes as the nets decay.
Both exponential (log-uniform) and logistic net retention are supported,
and repeated distributions accumulate correctly (net usage is the full
mean-field mixture over all past distributions).

### Seasonal SMC (needs a drug)

Chemoprevention and clinical treatment require a drug to be registered
first with `set_drugs()`. Here, four monthly SMC rounds in under-5s:

``` r

p_smc <- malariasimulation::set_drugs(
  malariasimulation::get_parameters(list(human_population = 1000)),
  list(malariasimulation::SP_AQ_params))
p_smc <- malariasimulation::set_smc(
  p_smc, drug = 1, timesteps = c(365, 395, 425, 455), coverages = rep(0.9, 4),
  min_ages = rep(round(0.25 * 365), 4), max_ages = rep(round(5 * 365), 4))
p_smc$prevalence_rendering_min_ages <- round(0.25 * 365)
p_smc$prevalence_rendering_max_ages <- round(5 * 365)

smc <- run_simulation_ode(3 * 365, malariasimulation::set_equilibrium(p_smc, init_EIR = 20))
pfpr_smc <- smc$n_detect_lm_91_1825 / smc$n_age_91_1825   # the 0.25-5y target band
c(baseline = pfpr_smc[1], trough = min(pfpr_smc))
#>   baseline     trough 
#> 0.47085538 0.02299911
```

SMC, MDA and PMC rounds are applied between days, as the IBM applies
them: of the target band, a fraction (coverage x drug efficacy) is
treated. The clinical and LM-detectable asymptomatic among them pass
through the treated state first, still detectable and infectious for a
few days, and everyone treated is then protected by a prophylaxis chain
matched to the drug’s Weibull protection curve.

### PEV vaccine (RTS,S via EPI, with a booster)

``` r

p_pev <- malariasimulation::set_pev_epi(
  malariasimulation::get_parameters(list(human_population = 1000)),
  profile = malariasimulation::rtss_profile,
  timesteps = 1, coverages = 0.9, min_wait = 0, age = round(6 * 30),
  booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
  booster_profile = list(malariasimulation::rtss_booster_profile))

pev <- run_simulation_ode(12 * 365, malariasimulation::set_equilibrium(p_pev, init_EIR = 20))
pfpr_pev <- pev$n_detect_lm_730_3650 / pev$n_age_730_3650
c(baseline = pfpr_pev[1], year12 = pfpr_pev[nrow(pev)])
#>  baseline    year12 
#> 0.5493395 0.5109142
```

PEV reduces the infection hazard by a per-age, per-time factor. Efficacy
is averaged over the per-individual antibody distribution the IBM
samples (a 4-D Gauss-Hermite rule over `cs`/`rho`/`ds`/`dl`), not
evaluated at the profile median; primary and booster doses, EPI and mass
campaigns, and time-varying EPI coverage are all modelled, and a
campaign’s protection ages with the cohort it vaccinated.
Transmission-blocking vaccines are available via `set_tbv()`.

## Demography

`set_demography()` supplies age-specific mortality, which reshapes the
equilibrium age structure the model is seeded at:

``` r

dr <- c(0.048, 0.007, 0.003, 0.004, 0.008, 0.020, 0.050, 0.120) / 365
ag <- round(c(1, 5, 10, 20, 40, 60, 80, 100) * 365)

p_dem <- malariasimulation::set_demography(
  malariasimulation::get_parameters(list(human_population = 1000)),
  agegroups = ag, timesteps = 0, deathrates = matrix(dr, nrow = 1))
p_dem$prevalence_rendering_min_ages <- c(0, 0)
p_dem$prevalence_rendering_max_ages <- c(5, 100) * 365

dem <- run_simulation_ode(8 * 365, malariasimulation::set_equilibrium(p_dem, init_EIR = 20))
under5_fraction <- dem$n_age_0_1825[nrow(dem)] / dem$n_age_0_36500[nrow(dem)]
under5_fraction     # ~0.09 here, vs ~0.21 under the default constant hazard
#> [1] 0.08659583
```

Custom demography is **time-varying**: `mu_age(t)` is interpolated over
`deathrate_timesteps`, so a demographic transition is modelled rather
than frozen. The baseline (`t = 0`) row additionally seeds the
equilibrium age structure. Match `default_age_lower(max_age =)` to the
top `deathrate_agegroups`; a warning fires if the model age grid runs
past them.

## Seasonality and burn-in

Turn on seasonality with the Fourier rainfall parameters. Unlike an
aseasonal run, a seasonal run **oscillates** around a limit cycle rather
than holding flat. The state is seeded at the annual-mean equilibrium,
so the first ~10 years are a transient onto the cycle, and the seasonal
annual-mean EIR sits about 7% below the aseasonal `init_EIR` target
(nonlinear averaging). Run long and read a **burned-in** year for
calibration or comparison:

``` r

p_seas <- malariasimulation::get_parameters(list(
  human_population = 1000, model_seasonality = TRUE,
  g0 = 0.285, g = c(-0.33, -0.13, 0.052), h = c(-0.35, 0.020, 0.10)))

seas <- run_simulation_ode(20 * 365, malariasimulation::set_equilibrium(p_seas, init_EIR = 20))
final_year <- seas[(nrow(seas) - 365 + 1):nrow(seas), ]
range(final_year$EIR)   # seasonal swing in the settled cycle
#> [1]  0.04853132 57.94258994
mean(final_year$EIR)    # about 7% below the aseasonal target of 20
#> [1] 18.54437
```

## P. vivax

A `get_parameters(parasite = "vivax")` list runs the *P. vivax* model:
hypnozoite batches and relapse, radical cure with liver-stage
protection, and no severe disease. The output is the IBM’s vivax table,
with relapses and hypnozoite carriage beside the familiar columns:

``` r

pv <- malariasimulation::get_parameters(list(
  human_population = 1000,
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * 365 - 1),
  parasite = "vivax")
pv <- malariasimulation::set_drugs(pv, list(malariasimulation::CQ_PQ_params_vivax))
pv <- malariasimulation::set_clinical_treatment(pv, drug = 1, timesteps = 1, coverages = 0.4)

out_pv <- run_simulation_ode(5 * 365, malariasimulation::set_equilibrium(pv, init_EIR = 3))
tail(out_pv[, c("timestep", "n_inc_clinical_0_36499", "n_relapses",
                "n_with_hypnozoites", "iaa_mean")], 3)
#>      timestep n_inc_clinical_0_36499 n_relapses n_with_hypnozoites iaa_mean
#> 1823     1823              0.8002892   7.645336           186.6042 9.237516
#> 1824     1824              0.8004321   7.645060           186.5968 9.236775
#> 1825     1825              0.8005750   7.644786           186.5895 9.236034
```

Radical cure starts on day 1 here, from a seed without it, so these rows
are still settling towards the treated state; a level worth comparing
wants a burn-in first
([`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)).
MDA, SMC and PMC are refused under vivax, because `malariasimulation`
itself fails on them.

## Discretisation settings

Every run above uses the default numerics, which live in
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)‘s
fourth argument, `tuning` — an
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
object, or a plain named list of just the fields you want to change.
Nothing epidemiological sits there: model parameters, the target EIR
included, stay on the parameter list. There are only three settings —
the age grid, the prophylaxis chains’ stage counts and the mosquito
sub-steps per day — and every default is the validated choice. The one
worth reaching for is the age grid. The default, 118 groups, is the
smallest on which every falciparum claim in `fleetcheck` passes; a grid
under half as fine runs over twice as fast, and is adequate where a
result does not rest on the level of severe incidence or on young
children at high transmission:

``` r

seas_coarse <- run_simulation_ode(
  20 * 365, malariasimulation::set_equilibrium(p_seas, init_EIR = 20),
  tuning = list(age_lower = default_age_lower(n_group = 53)))

c(default = mean(final_year$EIR),
  coarse  = mean(seas_coarse$EIR[(nrow(seas_coarse) - 365 + 1):nrow(seas_coarse)]))
#>  default   coarse 
#> 18.54437 18.53354
```

[`?ode_tuning`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
documents every field.

## Where to go next

Four things to read before you trust a number from this model.

|  |  |
|----|----|
| **[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)** | Where the mean field departs and what to do about it: burn-in, age grids, output bands, which settings to leave alone, which results to treat with caution, and what a run costs. |
| **[fleetcheck](https://pwinskill.github.io/fleetcheck/)** | How well it matches the IBM, as a separate project: a register of claims, each with the criterion that decides it, the value measured against it and a verdict. Not completed validation — the evidence that exists so far, with the open discrepancies named and any failing claim left failing. |
| **[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)** | The formal specification: the daily update, and what each state dimension carries as against what is captured *without* one. |
| **[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md)** | Every `malariasimulation` `set_*()` function, argument by argument. |

Scope: `fleet` runs whichever parasite the list names, as
`malariasimulation` does, except for chemoprevention under vivax, which
`malariasimulation` itself cannot run. The model is always
compartmental, so the individual-mosquito code path does not apply.

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] fleet_0.0.0.9003
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0                dplyr_1.2.1                  
#>  [3] compiler_4.6.1                tidyselect_1.2.1             
#>  [5] Rcpp_1.1.2                    stringr_1.6.0                
#>  [7] tidyr_1.3.2                   jquerylib_0.1.4              
#>  [9] systemfonts_1.3.2             textshaping_1.0.5            
#> [11] yaml_2.3.12                   fastmap_1.2.0                
#> [13] statmod_1.5.2                 malariaEquilibrium_1.0.1     
#> [15] R6_2.6.1                      generics_0.1.4               
#> [17] postie_1.1.0                  knitr_1.52                   
#> [19] MASS_7.3-65                   tibble_3.3.1                 
#> [21] desc_1.4.3                    monty_0.4.14                 
#> [23] bslib_0.12.0                  pillar_1.11.1                
#> [25] rlang_1.3.0                   stringi_1.8.9                
#> [27] cachem_1.1.0                  malariasimulation_3.0.0      
#> [29] dust2_0.3.28                  xfun_0.61                    
#> [31] fs_2.1.0                      sass_0.4.10                  
#> [33] otel_0.2.0                    cli_3.6.6                    
#> [35] withr_3.0.3                   pkgdown_2.2.1                
#> [37] magrittr_2.0.5                malariaEquilibriumVivax_1.0.1
#> [39] digest_0.6.39                 lifecycle_1.0.5              
#> [41] vctrs_0.7.3                   evaluate_1.0.5               
#> [43] glue_1.8.1                    ragg_1.5.2                   
#> [45] purrr_1.2.2                   rmarkdown_2.32               
#> [47] tools_4.6.1                   pkgconfig_2.0.3              
#> [49] htmltools_0.5.9
```
