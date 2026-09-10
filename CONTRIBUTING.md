# Contributing to blink

Thanks for your interest. `blink` is a deterministic mean-field (ODE) twin of the
[malariasimulation](https://github.com/mrc-ide/malariasimulation) individual-based
model. That framing drives almost every rule below: where the two models can be made
to agree, they should agree *by construction*, not by fitting.

> `blink` is a work in progress and is not ready for real use. See the warning at the
> top of the [README](https://github.com/pwinskill/blink#not-ready-for-real-use).

## The one thing that will catch you out

**`inst/odin/malaria_ode.R` is the model. `src/malaria_ode.cpp` is what actually runs.**

The C++ is generated from the odin source and is committed to the repository. If you
edit the odin file and do not regenerate, everything passes: `R CMD check` compiles
the old model, the tests exercise the old model, and the drift check reports "nothing
moved", which reads as *my change was numerically inert* rather than *my change was
never compiled*.

After any edit to `inst/odin/malaria_ode.R`:

```r
odin2::odin_package(".")
```

then commit the regenerated files alongside your change:

```
src/malaria_ode.cpp   src/cpp11.cpp   R/dust.R   R/cpp11.R   inst/dust/malaria_ode.cpp
```

CI enforces this: the `generated-code` job in `.github/workflows/R-CMD-check.yaml`
regenerates and fails on any diff.

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
"OK, 0 errors" having tested almost nothing. Set `BLINK_ALLOW_SKIP=1` to skip instead,
locally only.

## Changing the model

Two committed baselines pin `blink`'s numbers. A deliberate model change means
refreshing **both**, and reading both diffs.

| baseline | what it pins | regenerate with |
|---|---|---|
| `tests/testthat/reference-values.csv`<br>`tests/testthat/reference-interventions.csv` | absolute output levels, to 1e-6 | `BLINK_REGENERATE_REFERENCE=1 Rscript -e 'devtools::test(filter = "reference")'` |
| `comparison/data/rep_*.csv` | agreement with the IBM | `CMP_BLINK_ONLY=1 Rscript comparison/run_replicates.R` |

The IBM rows do not need re-running: nothing in `blink` can affect them. See
[comparison/README.md](https://github.com/pwinskill/blink/blob/main/comparison/README.md).

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
  are ignored. `blink` warns when it detects this.

## Pull requests

- Branch from `main`; keep changes focused.
- Say what moved. If output numbers changed, include the baseline diff and explain it.
- New behaviour needs a test that fails without the change. For anything touching model
  output, that means an absolute assertion, not only a relation: a test that checks
  "intervention-on is below intervention-off" survives multiplying the whole output by a
  constant, and two real mutations once passed the entire suite that way.
- Documentation counts as part of the change. If you falsify a claim in the README or a
  vignette, fix it in the same pull request.

## Reporting a problem

Open an issue with your `sessionInfo()`, the `blink` commit, and a reproducible example.
For anything numerical, say which parameter list you used and how you seeded it.
