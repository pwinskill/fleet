# Getting started with blink

`blink` is a fast, deterministic **mean-field (ODE) twin** of the
[malariasimulation](https://github.com/mrc-ide/malariasimulation)
individual-based model of *Plasmodium falciparum* malaria. It reproduces
the same age- and biting-heterogeneity-structured human model and
compartmental mosquito model, takes the **same parameter list** as the
IBM, and is seeded at the
[malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium)
fixed point. A multi-decade run finishes in a couple of seconds,
independent of population size, which makes it well suited to
calibration, sweeps and quick scenario work where a stochastic IBM would
be too slow.

Reach for the IBM (not this) when you need stochastic variation,
individual heterogeneity beyond the mean field, or *P. vivax* — see
**Scope** at the end.

## Install and load

``` r

# install.packages("remotes")
remotes::install_github("pwinskill/blink")
```

``` r

library(blink)
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
out <- run_simulation_ode(timesteps = 15 * 365, parameters = p, init_EIR = 20)

out[1:3, c("timestep", "n_age_730_3650", "n_detect_lm_730_3650",
           "p_detect_lm_730_3650", "EIR", "ft")]
#>   timestep n_age_730_3650 n_detect_lm_730_3650 p_detect_lm_730_3650 EIR ft
#> 1        0       284.2213             155.9074            0.5485424  20  0
#> 2        1       284.2213             155.8974            0.5485073  20  0
#> 3        2       284.2213             155.8880            0.5484741  20  0
```

The return value is a wide, `malariasimulation`-style daily count table.
Column tags carry the age band **in days**: `730_3650` is the canonical
2-10y band (always included), and `0_36500` is all ages. `p_detect_lm_*`
is LM-detectable prevalence, `EIR` is per adult per year, and `ft` is
the treated fraction. Because the model is mean-field, prevalence and
EIR are population-independent; `human_population` only rescales the
count columns.

With no interventions the model holds flat at the seeded equilibrium, so
the first and last rows agree:

``` r

range(out$EIR)                    # ~20 throughout
#> [1] 19.93161 20.02801
range(out$p_detect_lm_730_3650)   # constant to machine precision
#> [1] 0.5469938 0.5485424
```

## Reading outputs with postie

Because the count table is `malariasimulation`-shaped, it goes straight
into [postie](https://github.com/mrc-ide/postie) — the same two calls
you would make on an IBM run, with no blink-specific wrapper in between.

``` r

# prevalence: one <diagnostic>_prevalence_<lo>_<hi> column per age band (years)
prevalence <- postie::get_prevalence(out, diagnostic = "lm")
utils::tail(prevalence["lm_prevalence_2_10"], 3)
#>      lm_prevalence_2_10
#> 5474          0.5482791
#> 5475          0.5482791
#> 5476          0.5482791

# rates: clinical / severe incidence, mortality and DALYs by age band
rates <- postie::get_rates(out)
#> Warning in postie::get_rates(out): required column `ft_sev` (probability
#> hospitalisation | severe case) not found, assuming ft_sev = 0.8
utils::head(rates[, c("time", "age_lower", "age_upper",
                      "clinical", "severe", "dalys")])
#> # A tibble: 6 × 6
#>    time age_lower age_upper clinical    severe    dalys
#>   <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
#> 1 2000.         2        10  0.00337 0.0000444 0.000409
#> 2 2000.         0       100  0.00180 0.0000314 0.000168
#> 3 2000          2        10  0.00337 0.0000444 0.000412
#> 4 2000          0       100  0.00180 0.0000314 0.000169
#> 5 2000.         2        10  0.00337 0.0000444 0.000412
#> 6 2000.         0       100  0.00180 0.0000314 0.000169
```

Use `diagnostic = "pcr"` for PCR prevalence, and postie’s own arguments
(`scaler`, `treatment_scaler`, `life_expectancy`, …) as you normally
would. This means a post-processing pipeline written for the IBM works
on a `blink` run unchanged.

### Getting finer age resolution

By default you get the 2-10y and all-ages bands. To emit more, set the
`malariasimulation` rendering fields on the parameter list before
running — the same fields the IBM uses:

``` r

p_bands <- p
p_bands$prevalence_rendering_min_ages <- c(0,  2, 10) * 365
p_bands$prevalence_rendering_max_ages <- c(2, 10, 100) * 365
p_bands$clinical_incidence_rendering_min_ages <- c(0, 2) * 365
p_bands$clinical_incidence_rendering_max_ages <- c(2, 10) * 365

out_bands <- run_simulation_ode(10 * 365, p_bands, init_EIR = 20)
grep("^p_detect_lm_", names(out_bands), value = TRUE)
#> [1] "p_detect_lm_0_730"      "p_detect_lm_730_3650"   "p_detect_lm_3650_36500"
```

## Layering interventions

Interventions are configured with the ordinary
`malariasimulation::set_*` builders and applied automatically from the
list — there are no extra arguments to
[`run_simulation_ode()`](https://pwinskill.github.io/blink/reference/run_simulation_ode.md).
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

nets <- run_simulation_ode(12 * 365, p_nets, init_EIR = 20)
c(baseline = nets$p_detect_lm_730_3650[1],
  trough   = min(nets$p_detect_lm_730_3650))
#>  baseline    trough 
#> 0.5485424 0.1800651
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

smc <- run_simulation_ode(3 * 365, p_smc, init_EIR = 20)
band <- "p_detect_lm_91_1825"        # the 0.25-5y target band, tag in days
c(baseline = smc[[band]][1], trough = min(smc[[band]]))
#>    baseline      trough 
#> 0.470577730 0.009556474
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

pev <- run_simulation_ode(12 * 365, p_pev, init_EIR = 20)
c(baseline = pev$p_detect_lm_730_3650[1],
  year12   = pev$p_detect_lm_730_3650[nrow(pev)])
#>  baseline    year12 
#> 0.5485424 0.5113764
```

PEV reduces the infection hazard by a per-age, per-time factor. Efficacy
is averaged over the per-individual antibody distribution the IBM
samples (a 4-D Gauss-Hermite rule over `cs`/`rho`/`ds`/`dl`), not
evaluated at the profile median; primary and booster doses, EPI and mass
campaigns, and time-varying EPI coverage are all modelled.
Transmission-blocking vaccines are available via `set_tbv()`.

### Custom demography

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

dem <- run_simulation_ode(8 * 365, p_dem, init_EIR = 20)
under5_fraction <- dem$n_age_0_1825[nrow(dem)] / dem$n_age_0_36500[nrow(dem)]
under5_fraction     # ~0.09 here, vs ~0.21 under the default constant hazard
#> [1] 0.08787847
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

seas <- run_simulation_ode(20 * 365, p_seas, init_EIR = 20)
final_year <- seas[(nrow(seas) - 365 + 1):nrow(seas), ]
range(final_year$EIR)   # seasonal swing in the settled cycle
#> [1]  0.04953752 58.12753634
mean(final_year$EIR)    # a few % below the aseasonal target of 20
#> [1] 18.65741
```

## The validation story

`blink` is checked two ways. **First, self-consistency:** with no
interventions it holds flat at the `malariaEquilibrium` seed to machine
precision (and still holds flat with treatment, `ft > 0`, thanks to a
corrected prophylaxis-aging recursion). **Second, agreement with the
IBM:** the same parameter list is passed to both models —
`malariasimulation` burned in 30 years at a population of 10,000, this
model seeded at equilibrium — across six scenarios (`comparison/`).
Representative agreement:

- Equilibrium PfPR(2-10) vs EIR tracks the IBM closely (e.g. at EIR 20:
  IBM 0.553 vs ODE 0.547; at EIR 120: 0.788 vs 0.786).
- The age-prevalence profile at EIR 20 agrees to
  `max |IBM - ODE| = 0.008` (mean 0.004) across bands from infancy to
  85y.
- Custom demography reproduces the IBM’s under-5 population fraction
  (IBM 0.090 vs ODE 0.088).

Bed nets, seasonality and SMC scenarios are compared as well (including
clinical and severe incidence time series); see `comparison/plots/` for
the figures. LM prevalence and clinical incidence match the IBM to
within ~1%. **Severe incidence** (and DALYs derived from it) is more
sensitive to the immunity model and can differ by ~5-20% (largest at low
EIR); treat it as indicative.

## Scope and approximations

- **P. falciparum only.** *P. vivax* is rejected; the model is always
  compartmental (the individual-mosquito code path does not apply).
- **Mean-field caveats.** Vector control is population-averaged
  (slightly under-suppresses at deep troughs vs the IBM); drug
  prophylaxis is an Erlang chain matched to the drug’s Weibull
  protection (its mean and variance, not the Weibull itself). PCR
  prevalence follows the IBM convention (all D/Tr/A/U).
- **Full detail.**
  [`vignette("model")`](https://pwinskill.github.io/blink/articles/model.md)
  is the canonical specification: the ODE system, what each state
  dimension carries (and what is captured *without* a dimension), and an
  argument-by-argument table of every `malariasimulation` `set_*`
  function.

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
#> [1] blink_0.0.0.9000
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
#> [28] malariasimulation_3.0.0  dust2_0.3.28             xfun_0.60               
#> [31] fs_2.1.0                 sass_0.4.10              otel_0.2.0              
#> [34] cli_3.6.6                withr_3.0.3              pkgdown_2.2.1           
#> [37] magrittr_2.0.5           digest_0.6.39            lifecycle_1.0.5         
#> [40] vctrs_0.7.3              evaluate_1.0.5           glue_1.8.1              
#> [43] ragg_1.5.2               purrr_1.2.2              rmarkdown_2.32          
#> [46] tools_4.6.1              pkgconfig_2.0.3          htmltools_0.5.9
```
