# Contributing to fleet

Thanks for your interest. `fleet` is a deterministic mean-field twin of the
[malariasimulation](https://github.com/mrc-ide/malariasimulation) individual-based
model, advanced one day at a time on the IBM's own clock. That framing drives almost
every rule below: where the two models can be made to agree, they should agree *by
construction*, not by fitting.

> `fleet` is a work in progress and is not ready for real use. See the warning at the
> top of the [README](https://github.com/pwinskill/fleet#not-ready-for-real-use).

## The one thing that will catch you out

**`inst/odin/malaria_daily.R` is the model. `src/malaria_daily.cpp` is what actually runs.**

The C++ is generated from the odin source and is committed to the repository. If you
edit the odin file and do not regenerate, everything passes: `R CMD check` compiles
the old model, the tests exercise the old model, and the drift check reports "nothing
moved", which reads as *my change was numerically inert* rather than *my change was
never compiled*.

After any edit to `inst/odin/malaria_daily.R`:

```r
odin2::odin_package(".")
```

then commit the regenerated files alongside your change:

```
src/malaria_daily.cpp   src/cpp11.cpp   R/dust.R   R/cpp11.R   inst/dust/malaria_daily.cpp
```

CI enforces this: the `generated-code` job in `.github/workflows/R-CMD-check.yaml`
regenerates and fails on any diff.

**`devtools::load_all()` and `devtools::test()` compile the model with debug flags**
(`-O0`), which runs it about ten times slower, and they leave those object files in
`src/`. `R CMD INSTALL` on the source directory then reuses them, so an install made
after a `load_all()` is silently the slow build. Install with
`R CMD INSTALL --preclean .` (or delete `src/*.o` first) before timing anything.

## Setup

```r
install.packages("remotes")
remotes::install_deps(dependencies = TRUE)
```

None of the six modelling dependencies (`odin2`, `dust2`, `monty`, `malariaEquilibrium`,
`malariasimulation`, `postie`) are on CRAN. `DESCRIPTION`'s `Remotes` points all six at
GitHub, so that is where `install_deps()` fetches them from. You need a C++ toolchain
(Rtools on Windows).

## Running the tests

```r
devtools::test()
```

The suite requires `malariasimulation`. It **errors** rather than skips if that package
is missing, so a runner that cannot build the IBM goes red instead of reporting
"OK, 0 errors" having tested almost nothing. Set `FLEET_ALLOW_SKIP=1` to skip instead,
locally only.

## Changing the model

Two committed baselines pin `fleet`'s numbers: this repository's reference
values, and `fleetcheck`'s committed fleet rows. A deliberate model change means
refreshing **both**, and reading both diffs.

| baseline | what it pins | regenerate with |
|---|---|---|
| `tests/testthat/reference-values.csv`<br>`tests/testthat/reference-interventions.csv` | absolute output levels, to 1e-6 | `FLEET_REGENERATE_REFERENCE=1 Rscript -e 'testthat::test_file("tests/testthat/test-reference.R", package = "fleet", load_package = "installed")'`, against a fresh `R CMD INSTALL` |

Agreement with the IBM is a separate repository,
[fleetcheck](https://github.com/pwinskill/fleetcheck). A change that moves model
output moves its numbers too, so refresh it in the same change: clone it beside
this one and run `CMP_FLEET_ONLY=1 Rscript validations/02-scenarios/run.R`,
then `validations/02-scenarios/render.R`. The IBM rows do not need re-running —
nothing in `fleet` can affect them — which is why that takes under a minute rather
than the hours a full sweep costs.

**Never regenerate a baseline to make a red test green.** The diff *is* the record of
what your change did to the model. Read it, and put it in the pull request.

### House rules for the model itself

- **Replicate malariasimulation exactly.** Where the IBM has a known mechanism, mirror
  it, including its per-day discretisation and its arithmetic. Mean-matched
  approximations and fudge factors are not accepted, even when they improve agreement.
  If the mean field genuinely cannot carry a mechanism, document the departure and
  measure it rather than tuning around it.
- **Cite the source.** Most mean-field choices here name the malariasimulation file and
  line they mirror. Keep that up: it is what makes a divergence auditable later.
- **Deliberate departures carry their evidence.** `bite_dedup` and
  `acquired_immunity_offset` both exist because an A/B measurement justified them, and
  both record it.
- **`set_equilibrium()` must be the last call on a parameter list.** It stores
  `eq_params`, which freezes the translated biological constants; edits made afterwards
  are ignored. `fleet` warns when it detects this.

## Pull requests

- Branch from `main`; keep changes focused.
- Say what moved. If output numbers changed, include the baseline diff and explain it.
- New behaviour needs a test that fails without the change. For anything touching model
  output, that means an absolute assertion, not only a relation: a test that checks
  "intervention-on is below intervention-off" survives multiplying the whole output by a
  constant.
- Documentation counts as part of the change. If you falsify a claim in the README or a
  vignette, fix it in the same pull request.

## Reporting a problem

Open an issue with your `sessionInfo()`, the `fleet` commit, and a reproducible example.
For anything numerical, say which parameter list you used and how you seeded it.
