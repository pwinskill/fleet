# Getting started with fleet

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

`fleet` is a fast, deterministic **mean-field (ODE) twin** of the
[malariasimulation](https://github.com/mrc-ide/malariasimulation)
individual-based model of *Plasmodium falciparum* malaria. It reproduces
the same age- and biting-heterogeneity-structured human model and
compartmental mosquito model, takes the **same parameter list** as the
IBM, and is seeded at the
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium)
fixed point. A multi-decade run finishes in a couple of seconds,
independent of population size. That speed is the point: the intended
uses are calibration, sweeps and quick scenario work, where a stochastic
IBM would be too slow. It is not ready for those jobs yet (see the note
above).

Reach for the IBM (not this) when you need stochastic variation,
individual heterogeneity beyond the mean field, or *P. vivax*.

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

Build an unmodified `malariasimulation` parameter list and run it. You
supply a target adult EIR (infectious bites per adult per year); the
model seeds itself at the corresponding equilibrium, so no separate
burn-in is needed for an aseasonal run.

``` r

p <- malariasimulation::get_parameters(list(human_population = 1000))

# 15 years at a daily step, seeded at the equilibrium for EIR = 20
out <- run_simulation_ode(timesteps = 15 * 365, parameters = malariasimulation::set_equilibrium(p, init_EIR = 20))

out[1:3, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650",
           "p_detect_lm_730_3650", "EIR", "ft")]
#>   timestep n_age_730_3650 n_detect_lm_730_3650 p_detect_lm_730_3650 EIR ft
#> 1        0       288.0659             158.1277            0.5489289  20  0
#> 2        1       288.0659             158.1177            0.5488942  20  0
#> 3        2       288.0659             158.1083            0.5488614  20  0
```

The return value is a wide, `malariasimulation`-style daily count table.
Column tags carry the age band **in days**: `730_3650` is the canonical
2-10y band (always included), and `0_36500` is all ages. `p_detect_lm_*`
is LM-detectable prevalence, `EIR` is per adult per year, and `ft` is
the treated fraction. Because the model is mean-field, prevalence and
EIR are population-independent; `human_population` only rescales the
count columns.

With no interventions the run settles rather than holding perfectly
still:

``` r

range(out$EIR)                    # ~20 throughout
#> [1] 19.93045 20.00984
range(out$p_detect_lm_730_3650)   # a ~0.3% relaxation over the first years
#> [1] 0.5473664 0.5489289
```

Over these fifteen years LM prevalence in 2-10 year olds runs between
0.5485 (the seed, at day 0) and 0.5470 (its low point near the end of
the first year), a spread of 0.28%, and is flat to 0.02% over the last
five years. EIR stays within 0.5% of its target throughout.

That relaxation is expected, not a numerical error. The
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium) seed
solves a *simplified* form of the model, so it is a close
**approximation** to `fleet`’s true fixed point rather than the fixed
point itself; `malariasimulation` is seeded from the same solution and
drifts off it the same way.
[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)
lists the four simplifications. **Burn in before calibrating, and read a
settled window rather than day 0.**

## Reading outputs with postie

Because the count table is `malariasimulation`-shaped, it goes straight
into [postie](https://github.com/mrc-ide/postie): the same two calls you
would make on an IBM run, with no fleet-specific wrapper in between.

``` r

# prevalence: one <diagnostic>_prevalence_<lo>_<hi> column per age band (years)
prevalence <- postie::get_prevalence(out, diagnostic = "lm")
utils::tail(prevalence["lm_prevalence_2_10"], 3)
#>      lm_prevalence_2_10
#> 5474          0.5483602
#> 5475          0.5483602
#> 5476          0.5483602

# rates: clinical / severe incidence, mortality and DALYs by age band
rates <- postie::get_rates(out)
#> Warning in postie::get_rates(out): required column `ft_sev` (probability
#> hospitalisation | severe case) not found, assuming ft_sev = 0.8
utils::head(rates[, c("time", "age_lower", "age_upper",
                      "clinical", "severe", "dalys")])
#> # A tibble: 6 × 6
#>    time age_lower age_upper clinical    severe    dalys
#>   <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 2000.         2        10  0.00338 0.0000440 0.000406
#> 2 2000.         0       100  0.00183 0.0000315 0.000168
#> 3 2000          2        10  0.00338 0.0000440 0.000408
#> 4 2000          0       100  0.00183 0.0000315 0.000169
#> 5 2000.         2        10  0.00338 0.0000440 0.000408
#> 6 2000.         0       100  0.00183 0.0000315 0.000169
```

Use `diagnostic = "pcr"` for PCR prevalence, and postie’s own arguments
(`scaler`, `treatment_scaler`, `life_expectancy`, …) as you normally
would. This means a post-processing pipeline written for the IBM works
on a `fleet` run unchanged.

### Getting finer age resolution

By default you get the 2-10y and all-ages bands. To emit more, set the
`malariasimulation` rendering fields on the parameter list before
running (the same fields the IBM uses):

``` r

p_bands <- p
p_bands$prevalence_rendering_min_ages <- c(0,  2, 10) * 365
p_bands$prevalence_rendering_max_ages <- c(2, 10, 100) * 365
p_bands$clinical_incidence_rendering_min_ages <- c(0, 2) * 365
p_bands$clinical_incidence_rendering_max_ages <- c(2, 10) * 365

out_bands <- run_simulation_ode(10 * 365, malariasimulation::set_equilibrium(p_bands, init_EIR = 20))
grep("^p_detect_lm_", names(out_bands), value = TRUE)
#> [1] "p_detect_lm_0_730"      "p_detect_lm_730_3650"   "p_detect_lm_3650_36500"
```

## Layering interventions

Interventions are configured with the ordinary
`malariasimulation::set_*` builders and applied automatically from the
list. There are no extra arguments to
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md).
Every module starts at its baseline at t = 0, so the equilibrium seed is
preserved and each intervention acts from its scheduled timestep. A few
representative examples follow.

### Bed nets

``` r

p_nets <- malariasimulation::set_bednets(
  malariasimulation::get_parameters(list(human_population = 1000)),
  timesteps = 5 * 365, coverages = 0.8, retention = 5 * 365,
  dn0 = matrix(0.387, 1, 1), rn = matrix(0.563, 1, 1),
  rnm = matrix(0.24, 1, 1), gamman = 2.64 * 365)

nets <- run_simulation_ode(12 * 365, malariasimulation::set_equilibrium(p_nets, init_EIR = 20))
c(baseline = nets$p_detect_lm_730_3650[1],
  trough   = min(nets$p_detect_lm_730_3650))
#>  baseline    trough 
#> 0.5489289 0.1810685
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
band <- "p_detect_lm_91_1825"        # the 0.25-5y target band, tag in days
c(baseline = smc[[band]][1], trough = min(smc[[band]]))
#>   baseline     trough 
#> 0.46982293 0.01038361
```

SMC/MDA/PMC are applied as pulsed mass drug administration between
integration segments: a fraction (coverage x drug efficacy) of the
target age band is cleared and moved to the first stage of a
chemoprevention prophylaxis chain matched to the drug’s Weibull
protection curve.

### PEV vaccine (RTS,S via EPI, with a booster)

``` r

p_pev <- malariasimulation::set_pev_epi(
  malariasimulation::get_parameters(list(human_population = 1000)),
  profile = malariasimulation::rtss_profile,
  timesteps = 1, coverages = 0.9, min_wait = 0, age = round(6 * 30),
  booster_spacing = round(12 * 30), booster_coverage = matrix(0.8, 1, 1),
  booster_profile = list(malariasimulation::rtss_booster_profile))

pev <- run_simulation_ode(12 * 365, malariasimulation::set_equilibrium(p_pev, init_EIR = 20))
c(baseline = pev$p_detect_lm_730_3650[1],
  year12   = pev$p_detect_lm_730_3650[nrow(pev)])
#>  baseline    year12 
#> 0.5489289 0.5117402
```

PEV reduces the infection hazard by a per-age, per-time factor. Efficacy
is averaged over the per-individual antibody distribution the IBM
samples (a 4-D Gauss-Hermite rule over `cs`/`rho`/`ds`/`dl`), not
evaluated at the profile median; primary and booster doses, EPI and mass
campaigns, and time-varying EPI coverage are all modelled.
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
annual-mean EIR sits a few percent below the aseasonal `init_EIR` target
(nonlinear averaging). Run long and read a **burned-in** year for
calibration or comparison:

``` r

p_seas <- malariasimulation::get_parameters(list(
  human_population = 1000, model_seasonality = TRUE,
  g0 = 0.285, g = c(-0.33, -0.13, 0.052), h = c(-0.35, 0.020, 0.10)))

seas <- run_simulation_ode(20 * 365, malariasimulation::set_equilibrium(p_seas, init_EIR = 20))
final_year <- seas[(nrow(seas) - 365 + 1):nrow(seas), ]
range(final_year$EIR)   # seasonal swing in the settled cycle
#> [1]  0.04947592 58.08644424
mean(final_year$EIR)    # a few % below the aseasonal target of 20
#> [1] 18.64028
```

## Solver settings

Every run above uses the default numerics, which live in
[`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)’s
fourth argument, `tuning` — an
[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
object, or a plain named list of just the fields you want to change.
Nothing epidemiological sits there: model parameters, the target EIR
included, stay on the parameter list. Every default is the validated
choice, so leave `tuning` alone except on a long projection like the
20-year run above, where the solver rather than the daily output grid
sets the step count. There, loosening the relative tolerance runs about
1.4–1.7x faster on a seasonal projection and moves the aggregate outputs
far less than anything you would report.

``` r

seas_fast <- run_simulation_ode(
  20 * 365, malariasimulation::set_equilibrium(p_seas, init_EIR = 20),
  tuning = list(rtol = 1e-6, step_size_max = 10))

c(default = mean(final_year$EIR),
  faster  = mean(seas_fast$EIR[(nrow(seas_fast) - 365 + 1):nrow(seas_fast)]))
#>  default   faster 
#> 18.64028 18.64028
```

[`?ode_tuning`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
documents every field, including why `atol` should stay at `1e-8`
whatever else you change, and why `step_size_max` is a safety rail
rather than a speed control.

## Where to go next

Four things to read before you trust a number from this model.

|  |  |
|----|----|
| **[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)** | Where the mean field departs and what to do about it: burn-in, age grids, output bands, which settings to leave alone, which results to treat with caution, and what a run costs. |
| **[`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)** | How well it matches the IBM: the same parameter list through both models across eighteen scenarios, plus a 63-country site-file comparison. Not completed validation — the evidence that exists so far, with the open discrepancies named. |
| **[`vignette("model")`](https://pwinskill.github.io/fleet/articles/model.md)** | The formal specification: the ODE system, and what each state dimension carries as against what is captured *without* one. |
| **[`vignette("parameters")`](https://pwinskill.github.io/fleet/articles/parameters.md)** | Every `malariasimulation` `set_*()` function, argument by argument. |

Scope: `fleet` is **P. falciparum only**. *P. vivax* parameter lists are
rejected at input, and the model is always compartmental, so the
individual-mosquito code path does not apply.

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
#> [1] fleet_0.0.0.9002
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0           dplyr_1.2.1              compiler_4.6.1          
#>  [4] tidyselect_1.2.1         Rcpp_1.1.2               stringr_1.6.0           
#>  [7] tidyr_1.3.2              jquerylib_0.1.4          systemfonts_1.3.2       
#> [10] textshaping_1.0.5        yaml_2.3.12              fastmap_1.2.0           
#> [13] statmod_1.5.2            malariaEquilibrium_1.0.1 R6_2.6.1                
#> [16] generics_0.1.4           postie_1.1.0             knitr_1.52              
#> [19] MASS_7.3-65              tibble_3.3.1             desc_1.4.3              
#> [22] monty_0.4.14             bslib_0.12.0             pillar_1.11.1           
#> [25] rlang_1.3.0              stringi_1.8.9            cachem_1.1.0            
#> [28] malariasimulation_3.0.0  dust2_0.3.28             xfun_0.61               
#> [31] fs_2.1.0                 sass_0.4.10              otel_0.2.0              
#> [34] cli_3.6.6                withr_3.0.3              pkgdown_2.2.1           
#> [37] magrittr_2.0.5           digest_0.6.39            lifecycle_1.0.5         
#> [40] vctrs_0.7.3              evaluate_1.0.5           glue_1.8.1              
#> [43] ragg_1.5.2               purrr_1.2.2              rmarkdown_2.32          
#> [46] tools_4.6.1              pkgconfig_2.0.3          htmltools_0.5.9
```
