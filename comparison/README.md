# Comparison harness — blink vs malariasimulation

Scripts that run every comparison scenario through **both** models on the *same*
parameter list and render the figures used in `vignette("comparison")` and the
package README. The IBM (`malariasimulation`) is run as `N_REP` stochastic
replicates per scenario (10,000 people, 30-year burn-in) in parallel; `blink` is
run once, seeded at the `malariaEquilibrium` fixed point. Rendering age bands are
set once on the shared parameter list so both models emit identical
`n_detect_lm_*` / `n_age_*` / `n_inc_*` columns, and one summariser reduces both.

## Files

| File | What it does |
| --- | --- |
| `theme.R` | House style (palette, theme, series scales, envelope helper, `save_fig()`) **and** the scenario constants (`BURN_Y`, `POP`, `N_REP`, `EIR_GRID`, age bands, seasonality, intervention labels). Sourced by every other script so they cannot drift. |
| `run_replicates.R` | Builds the scenarios (EIR grid, seasonal, custom demography, five interventions), runs blink and the IBM replicates on a PSOCK cluster, writes `data/rep_{eq,age,monthly,doy,timing}.csv`. `CMP_SMOKE=1` gives a 4-year, 1-replicate end-to-end check into `data/smoke/`; `CMP_ONLY=a,b` re-runs a subset and merges it into the existing CSVs. |
| `render_figures.R` | Draws every `cmp_*.png` from the CSVs (no model runs) into `man/figures/` and `vignettes/`. The 63-country panel reads the `blink2_validate` results directly. |
| `summary_tables.R` | The numbers quoted in the article, as markdown tables in `data/tables.md`. |

## Reproduce

```bash
Rscript comparison/run_replicates.R     # ~25 min on 10 workers
Rscript comparison/render_figures.R     # seconds
Rscript comparison/summary_tables.R     # seconds
```

(The scripts pin the Windows-arm64 R library path at the top; adjust `.libPaths()`
and `ROOT` for another machine. `N_WORKERS` is set in `run_replicates.R`.)

## Design notes

- **Replicates, not a realisation.** The IBM is always drawn as the median of its
  replicates with a 10–90% band; a single noisy run would misstate both the
  agreement and the disagreement.
- **Pooled rates.** Prevalence and incidence are pooled as `sum(cases) /
  sum(person-time)` over each window, never as a mean of per-day ratios.
- **Fair EIR axis.** The IBM's `EIR_<species>` is total infectious bites per day;
  `/ population * 365` recovers blink's per-adult-per-year convention.
- **Rendering cost.** The IBM's per-band output rendering dominates its run time,
  so only the reference (EIR 20) and demography scenarios carry the 12-band age
  profile; every other scenario renders just the default 2–10 y prevalence and
  0–5 y incidence bands.
- **Greyscale/CVD-safe figures.** Colour, line type and point shape all encode
  the model; the palette (indigo/coral, the site colours) passes the dataviz
  validator for CVD separation and contrast. Text never wears a series colour.
