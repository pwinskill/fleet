# Contributing to fleet

`fleet` is a deterministic mean-field counterpart to the
[malariasimulation](https://github.com/mrc-ide/malariasimulation)
individual-based model (IBM), advanced one day at a time on the IBM’s
clock. Most of the rules below follow from that: where the two models
can be made to agree, they should agree by construction, not by fitting.

`fleet` is a work in progress and not ready for real use; see the
[README](https://github.com/pwinskill/fleet#not-ready-for-real-use).

## Generated code

`inst/odin/malaria_daily.R` is the model source; `src/malaria_daily.cpp`
is what runs. The C++ is generated from the odin source and committed to
the repository. If you edit the odin file and do not regenerate,
everything still passes: `R CMD check` compiles the old model, the tests
exercise the old model, and the drift check reports that nothing moved,
which reads as a numerically inert change rather than one that was never
compiled.

After any edit to `inst/odin/malaria_daily.R`, regenerate:

``` r

odin2::odin_package(".")
```

and commit the regenerated files with the change:

    src/malaria_daily.cpp   src/cpp11.cpp   R/dust.R   R/cpp11.R   inst/dust/malaria_daily.cpp

The `generated-code` job in `.github/workflows/R-CMD-check.yaml`
regenerates the code on Linux and fails on any difference. odin2 prints
a fractional literal to 17 digits, which differs between platforms
unless the number is a short binary fraction: 0.5 is safe, 1e-12 is not.
Write any other constant as a quotient of whole numbers, `1 / 1e12`,
which prints the same everywhere; `test-odin-literals.R` checks this
locally.

## Setup

``` r

install.packages("remotes")
remotes::install_deps(dependencies = TRUE)
```

None of the seven modelling dependencies (`odin2`, `dust2`, `monty`,
`malariaEquilibrium`, `malariaEquilibriumVivax`, `malariasimulation`,
`postie`) is on CRAN. `Remotes` in `DESCRIPTION` points them at GitHub,
which is where `install_deps()` fetches them from. A C++ toolchain is
required (Rtools on Windows).

## Tests

``` r

devtools::test()
```

The suite requires `malariasimulation`. It errors rather than skips if
that package is missing, so that a runner that cannot build the IBM
fails instead of reporting success having tested almost nothing. Set
`FLEET_ALLOW_SKIP=1` to skip instead, locally only.

`devtools::load_all()` and `devtools::test()` compile the model with
debug flags (`-O0`), which makes it about ten times slower, and leave
those object files in `src/`. `R CMD INSTALL` on the source directory
then reuses them, so an install made after `load_all()` is the slow
build without saying so. Install with `R CMD INSTALL --preclean .`, or
delete `src/*.o` first, before timing anything.

## Changing the model

Two committed baselines pin `fleet`’s output, so a deliberate model
change means refreshing both and reading both diffs.

The first is this repository’s reference values, which pin absolute
output levels to 1e-6: `tests/testthat/reference-values.csv`,
`reference-interventions.csv`, `reference-vivax.csv` and
`reference-vivax-interventions.csv`. Regenerate them against a fresh
`R CMD INSTALL`:

    FLEET_REGENERATE_REFERENCE=1 Rscript -e 'testthat::test_file("tests/testthat/test-reference.R", package = "fleet", load_package = "installed")'

The second is the committed `fleet` rows in
[fleetcheck](https://github.com/pwinskill/fleetcheck), the separate
repository that compares `fleet` with the IBM. A change that moves model
output moves these too, so refresh them in the same change: clone
`fleetcheck` beside this repository and run
`CMP_FLEET_ONLY=1 Rscript validations/02-scenarios/run.R`, then
`validations/02-scenarios/render.R`. For a change that reaches the vivax
block, repeat with `CMP_PARASITE=pv`, then `render_pv.R`. The IBM rows
do not need re-running, since nothing in `fleet` can affect them, so
this takes minutes rather than the hours of a full sweep.

Never regenerate a baseline to make a failing test pass. The diff is the
record of what your change did to the model: read it, and include it in
the pull request.

### Model conventions

- Replicate `malariasimulation` exactly. Where the IBM has a known
  mechanism, mirror it, including its per-day discretisation and its
  arithmetic. Mean-matched approximations and fudge factors are not
  accepted, even when they improve agreement. If the mean field cannot
  carry a mechanism, document the departure and measure it rather than
  tuning around it.
- Cite the source. Most mean-field choices name the `malariasimulation`
  file and line they mirror; keep doing so, as it is what makes a
  divergence auditable later.
- Deliberate departures carry their evidence. `bite_dedup` and
  `acquired_immunity_offset` both exist because an A/B measurement
  justified them, and both record it.
- `set_equilibrium()` must be the last call on a parameter list. It
  stores `eq_params`, which freezes the translated biological constants,
  so edits made afterwards are ignored. `fleet` warns when it detects
  this.

## Pull requests

- Branch from `main` and keep changes focused.
- State what moved. If output changed, include the baseline diff and
  explain it.
- New behaviour needs a test that fails without the change. For anything
  that touches model output, that means an absolute assertion, not only
  a relation: a test that checks that output with an intervention is
  below output without it still passes if the whole output is multiplied
  by a constant.
- Documentation is part of the change. If a change falsifies a claim in
  the README or a vignette, fix it in the same pull request.

## Reporting a problem

Open an issue with your
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), the `fleet`
commit, and a reproducible example. For anything numerical, say which
parameter list you used and how you seeded it.
