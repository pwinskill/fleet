# Getting started with fleet

`fleet` is not yet validated for real use; see the [status
note](https://pwinskill.github.io/fleet/#not-ready-for-real-use).

This article runs the model on a `malariasimulation` parameter list,
reads the output with `postie`, and adds interventions, custom
demography, seasonality, *P. vivax* and the discretisation settings. It
assumes familiarity with `malariasimulation`.

``` r

# install.packages("remotes")
remotes::install_github("pwinskill/fleet")
```

``` r

library(fleet)
```

`malariasimulation` builds the parameter lists and `postie`
post-processes the output. Both are in Suggests and are installed
separately (see the
[README](https://pwinskill.github.io/fleet/#installation)).

## Running the model

Build a parameter list and set the target EIR, in infectious bites per
adult per year, with `set_equilibrium()`, as for the IBM; the run is
seeded at the corresponding equilibrium. The output age bands are also
set on the list: here clinical and severe incidence in under-5s and at
all ages, beside the default 2-10 prevalence band.

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

The result is `malariasimulation`’s daily table: the columns the IBM
would return for this parameter list, with the same names and meanings.
Column tags give the age band in days, both ends included, so `730_3650`
is 2-10 years and `0_1824` is under 5. `n_*` columns are counts:
`n_detect_lm_*` is the number LM-positive and `n_inc_clinical_*` the
day’s new clinical cases in the band. `EIR` is per adult per year.
Because the model is a mean field, rates do not depend on
`human_population`, which only scales the counts.

With no interventions the run settles rather than holding exactly still:

``` r

pfpr <- out$n_detect_lm_730_3650 / out$n_age_730_3650   # LM prevalence, 2-10y
range(out$EIR)   # ~20 throughout
#> [1] 19.92115 20.00000
range(pfpr)      # a ~0.4% relaxation over the first years
#> [1] 0.5473898 0.5493395
```

Row 1 is day 1, whose state is the seed. LM prevalence in 2-10 year olds
falls from the seed by 0.35% to a low a little over a year in, and is
flat to within 0.008% over the last five years. EIR stays within 0.4% of
its target. Incidence moves further from the seed than prevalence does:
a few per cent over the first years at this EIR, and more at higher
ones.

This relaxation is expected. The
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium) seed
solves the model’s continuous-time form, so it approximates `fleet`’s
fixed point on the daily clock rather than being that fixed point.
`malariasimulation` is seeded from the same solution and drifts from it
too;
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
lists the differences. Burn in before calibrating, and read a settled
window rather than the first days.

## Output

The table goes into [postie](https://github.com/mrc-ide/postie) with the
same calls as an IBM run:

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

Use `diagnostic = "pcr"` for PCR prevalence, and `postie`’s other
arguments (`scaler`, `treatment_scaler`, `life_expectancy` and so on) as
for an IBM run.

### Age bands

Bands come from the `malariasimulation` rendering fields on the
parameter list, and each output family is reported over its own fields
only: `prevalence_rendering_*` for detection;
`clinical_incidence_rendering_*`, `severe_incidence_rendering_*` and
`incidence_rendering_*` for the three incidence families; and
`age_group_rendering_*` for population alone. A band need not align with
the model’s age groups:

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

## Interventions

Interventions are set with the usual `malariasimulation::set_*()`
functions and read from the parameter list;
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
has no extra arguments for them. As in the IBM, settings read on the day
(treatment coverage, the drug mix, resistance, death rates) change on
their own timestep, and deployments (nets, spraying, vaccination,
chemoprevention) act from the day after theirs.
[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md)
covers every argument.

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
Exponential (log-uniform) and logistic net retention are both supported.
Net usage is the mixture over all past distributions, so repeated
distributions accumulate.

### Seasonal malaria chemoprevention

Chemoprevention and clinical treatment need a drug registered first with
`set_drugs()`. Here, four monthly SMC rounds in children aged 3 months
to 5 years:

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

SMC, MDA and PMC rounds are applied between days, as in the IBM. A
fraction of the target band equal to coverage times drug efficacy is
treated. The clinical and LM-detectable asymptomatic among them pass
through the treated state first, still detectable and infectious for a
few days, and everyone treated is then protected by a prophylaxis chain
matched to the drug’s Weibull protection curve.

### Pre-erythrocytic vaccination

RTS,S through the EPI schedule, with a booster:

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

PEV scales the infection hazard by an age- and time-dependent factor.
Efficacy is averaged over the per-individual antibody distribution the
IBM samples (a 4-D Gauss-Hermite rule over `cs`, `rho`, `ds` and `dl`)
rather than evaluated at the profile median. Primary and booster doses,
EPI and mass campaigns, and time-varying EPI coverage are modelled, and
a campaign’s protection ages with the cohort it vaccinated.
Transmission-blocking vaccines are set with `set_tbv()`.

## Demography

`set_demography()` sets age-specific mortality, which also reshapes the
equilibrium age structure the run is seeded at:

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

Custom demography is time-varying: `mu_age(t)` is interpolated over
`deathrate_timesteps`, so a demographic transition is modelled rather
than frozen. The `t = 0` row also sets the seeded age structure. Match
`default_age_lower(max_age =)` to the top of `deathrate_agegroups`; a
warning is given if the model’s age grid runs past it.

## Seasonality

Seasonality is set with the Fourier rainfall parameters. A seasonal run
oscillates around a limit cycle rather than settling to a fixed point.
It is seeded at the annual-mean equilibrium, so the first ~10 years are
a transient onto the cycle, and the annual-mean EIR on the cycle is
about 7% below the aseasonal `init_EIR` target, through nonlinear
averaging. Read a burned-in year for calibration or comparison:

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

## *P. vivax*

A list from `get_parameters(parasite = "vivax")` runs the *P. vivax*
model: hypnozoite batches and relapse, radical cure with liver-stage
protection, and no severe disease. The output is the IBM’s vivax table,
which adds relapses and hypnozoite carriage to the usual columns:

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
#> 1823     1823              0.8015565   7.658682           186.7889 9.239522
#> 1824     1824              0.8016998   7.658410           186.7815 9.238782
#> 1825     1825              0.8018432   7.658138           186.7742 9.238043
```

Radical cure starts on day 1 here, from a seed without it, so these rows
are still settling towards the treated state; burn in before comparing
levels
([`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)).
MDA, SMC and PMC are refused under *P. vivax*, because
`malariasimulation` itself fails on them.

## Discretisation settings

Numerical settings are the fourth argument of
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md),
`tuning`: an
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
object, or a named list of the fields to change. They are the age grid,
the stage counts of the prophylaxis chains, the mosquito sub-steps per
day and, for *P. vivax*, whether hypnozoite batch counts share levels.
Model parameters, the target EIR included, stay on the parameter list,
and every default is the validated choice.

The age grid is the setting most worth changing. The default of 118
groups is the smallest on which every falciparum claim in `fleetcheck`
passes. A grid under half as fine runs over twice as fast, and is
adequate where a result does not rest on the level of severe incidence
or on young children at high transmission:

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

## Next steps

Before relying on a result, read
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
for where the mean field departs from the IBM and what to do about it,
and check the relevant claims in
[fleetcheck](https://pwinskill.github.io/fleetcheck/).

## Session information

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
#> [1] fleet_0.0.0.9006
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
