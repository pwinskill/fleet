# Comparison harness: fleet vs malariasimulation

Scripts that run every comparison scenario through **both** models on the *same*
parameter list and render the figures used in `vignette("comparison")` and the
package README. The IBM (`malariasimulation`) is run as `N_REP` stochastic
replicates per scenario (10,000 people, 30-year burn-in) in parallel; `fleet` is
run once, seeded at the `malariaEquilibrium` fixed point. Rendering age bands are
set once on the shared parameter list so both models emit identical
`n_detect_lm_*` / `n_age_*` / `n_inc_*` columns, and one summariser reduces both.

## Files

| File | What it does |
| --- | --- |
| `constants.R` | Paths (`ROOT`, `VDIR`) and the scenario constants (`BURN_Y`, `POP`, `N_REP`, `EIR_GRID`, age bands, seasonality, intervention labels). Sourced by every other script, so they cannot drift. Split out of `theme.R` so the two scripts that draw nothing — the IBM run and the drift check — need no plotting stack. |
| `theme.R` | House style: palette, theme, series scales, envelope helper, `save_fig()`. Sources `constants.R`. |
| `scenarios.R` | The scenario definitions, the shared per-run summariser, and `run_fleet()`. Sourced by both `run_replicates.R` and `check_drift.R`, so a drift check cannot silently test different scenarios from the ones the reference was built on. |
| `check_drift.R` | **The one to run often.** Re-runs fleet only (~2 min) against the committed IBM rows. Reports movement and agreement separately; exits non-zero on drift. |
| `run_replicates.R` | Runs fleet and the IBM replicates on a PSOCK cluster, writes `data/rep_{eq,age,monthly,doy,timing}.csv` and `data/ibm_reference.json`. `CMP_SMOKE=1` gives a 4-year, 1-replicate end-to-end check into `data/smoke/`; `CMP_ONLY=a,b` re-runs a subset and merges it into the existing CSVs. |
| `render_figures.R` | Draws every `cmp_*.png` from the CSVs (no model runs) into `man/figures/` and `vignettes/`. The 63-country hex panel is a **snapshot** and is not redrawn: it needs a ~7-hour validation run in a separate checkout. `CMP_REFRESH_SITES=1` re-takes it. |
| `summary_tables.R` | The numbers quoted in the article, as markdown tables in `data/tables.md`. The site-file statistics come from the committed `data/site_snapshot.json`, not from a live run, so they survive on a machine without the validation checkout. |

## Checking a change: don't re-run the IBM

The IBM does not depend on fleet, so its committed rows stay valid for any
fleet-side change. Re-run fleet alone and compare:

```bash
Rscript comparison/check_drift.R          # ~2 min, exits 1 on drift
CMP_ONLY=eir_20,smc Rscript comparison/check_drift.R   # a subset, ~20 s
CMP_STRICT=1 Rscript comparison/check_drift.R          # also fail if ANY value moved
```

It answers two questions separately, because they mean different things.
**Did anything move?** fleet now against fleet's committed rows: a moved number
is not automatically wrong (a deliberate model fix moves numbers) but it must be
seen, and silent movement is how a regression ships. **Is the match still good?**
fleet now against the committed IBM medians and 10–90% bands, at the thresholds
the article's claims rest on. Only the second one fails the run by default.

If the movement was intended, refresh the committed fleet rows without touching
the IBM:

```bash
CMP_FLEET_ONLY=1 Rscript comparison/run_replicates.R   # ~1 min
```

### The unit-test suite pins fleet's numbers too — refresh both

These committed rows are not the only baseline of fleet's own output.
`tests/testthat/test-reference.R` holds a second one: absolute levels for one
quantity from every output family, at three EIRs, pinned to 1e-6 in
`tests/testthat/reference-values.csv` and checked on every `devtools::test()` and
every `R CMD check` run. It answers a different question — *did anything move at
all*, rather than *is the match to the IBM still good* — but it is a snapshot of
the same model, so a deliberate model change makes both stale at once. Regenerate
both, in the commit that makes the change:

```bash
FLEET_REGENERATE_REFERENCE=1 Rscript -e 'devtools::test(filter = "reference")'
CMP_FLEET_ONLY=1 Rscript comparison/run_replicates.R
```

Do one and not the other and the one left behind stops carrying information:
stale comparison rows make `check_drift.R` report movement that was reviewed and
accepted weeks ago, and a stale CSV turns the test suite red for the same reason
— which is the pressure that gets a reference regenerated to make a red test
green. Read both diffs; between them they are the record of what the change did.
`.github/workflows/CI.md` has the same note from the CI side.

### When the IBM *does* need re-running

Only when the reference itself goes stale, which `check_drift.R` tells you about
rather than leaving you to remember. `data/ibm_reference.json` records the date,
the malariasimulation version, `N_REP`/`POP`/burn-in, and a digest of the
scenario definitions (every parameter that differs from a bare
`get_parameters()`, plus each EIR and horizon). The check warns when the
installed malariasimulation or the scenario digest no longer matches.

## Full re-run

```bash
Rscript comparison/run_replicates.R     # ~25 min on 10 workers; re-stamps ibm_reference.json
Rscript comparison/render_figures.R     # seconds
Rscript comparison/summary_tables.R     # seconds
```

(No paths are hardcoded. Each script finds the checkout root by walking up to the
`DESCRIPTION`, so it runs from any working directory, via `Rscript` or `source()`.
Two optional environment variables: `FLEET_LIB` is prepended to the library path
(the libraries already on it are kept), for installations that do not pick up
`R_LIBS_USER`; `FLEET_VALIDATE` points at the site-file validation results, which
default to `../fleet_validate` and are skipped when absent. `N_WORKERS` is set
in `run_replicates.R`.)

## Design notes

- **Replicates, not a realisation.** The IBM is always drawn as the median of its
  replicates with a 10–90% band; a single noisy run would misstate both the
  agreement and the disagreement.
- **Pooled rates.** Prevalence and incidence are pooled as `sum(cases) /
  sum(person-time)` over each window, never as a mean of per-day ratios.
- **Fair EIR axis.** The IBM's `EIR_<species>` is total infectious bites per day;
  `/ population * 365` recovers fleet's per-adult-per-year convention.
- **Rendering cost.** The IBM's per-band output rendering dominates its run time,
  so only the reference (EIR 20) and demography scenarios carry the 12-band age
  profile; every other scenario renders just the default 2–10 y prevalence and
  0–5 y incidence bands.
- **Greyscale/CVD-safe figures.** Colour, line type and point shape all encode
  the model; the palette (indigo/coral, the site colours) passes the dataviz
  validator for CVD separation and contrast. Text never wears a series colour.
