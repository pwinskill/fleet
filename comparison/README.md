# Validation harness — blink vs malariasimulation

Scripts that run each scenario through **both** models on the *same* parameter
list and render IBM-vs-ODE figures. The IBM (`malariasimulation`) is burned in
~30 years at a population of 10,000 to reach its stochastic steady state; the ODE
(`blink`) is seeded at the `malariaEquilibrium` fixed point. Prevalence /
incidence age bands are set once on the shared parameter list, so both models
emit identical `n_detect_lm_*` / `n_age_*` / `n_inc_*` columns.

## Files

| File | What it does |
| --- | --- |
| `run_comparison.R` | Six scenarios: **A** PfPR–EIR curve, **B** age-prevalence profile, **C** bed-net campaign, **D** seasonal cycle, **E** seasonal SMC, **F** custom demography. Writes `data/*.csv` and `plots/*.png`. |
| `run_incidence.R` | Clinical & severe **incidence** time series (under-5, monthly-annualised) for the nets / seasonal / SMC scenarios. Writes `data/inc_*.csv` and `plots/inc_*.png`. |
| `plot_helpers.R` | Shared figure builders + constants (population, burn-in, age bands, colour/shape/linetype, incidence specs). Sourced by both runners **and** `render_plots.R`, so the live and re-render paths never drift. |
| `render_plots.R` | Re-render every PNG from the saved CSVs (no model runs) — cheap iteration on figure design. |

## Reproduce

```bash
Rscript comparison/run_comparison.R    # ~15 min (six 30-year IBM burn-ins)
Rscript comparison/run_incidence.R     # ~4 min (three 30-year IBM burn-ins)
Rscript comparison/render_plots.R      # seconds (re-plot from CSVs)
```

(The scripts pin the Windows-arm64 R library path at the top; adjust `.libPaths()`
and the `ROOT` path for another machine.)

## Design notes

- **Pooled rates.** Prevalence and incidence are pooled as `sum(cases) / sum(person-time)`
  over the observation window, not a mean of per-day ratios.
- **Fair EIR axis.** The IBM's rendered `EIR_<species>` is total infectious bites;
  dividing by population and ×365 recovers the ODE's per-adult-per-year convention.
- **Greyscale/CVD-safe figures.** The IBM/ODE distinction is carried by colour,
  line type **and** point shape, so the panels read in monochrome and under
  colour-vision deficiency. The noisier (rare-event) IBM severe series is drawn at
  reduced opacity so it doesn't obscure the ODE mean.
