<!-- Two things about where this file lives, both learned the hard way.

     NOT named README.md: GitHub resolves a repository's front page in the order
     .github/README.md, README.md, docs/README.md, so a file at .github/README.md
     silently REPLACES the package README on the repo home page. It did, briefly.

     NOT in .github/ directly, but one level down in .github/workflows/: pkgdown's
     package_mds() globs every *.md in the package root AND in .github/, with only
     a hardcoded skip-list, and renders each one as a site page. This file was
     being published at /CI.html -- linked from nothing, but in the sitemap and
     the site search. There is no config key to exclude it. The glob is not
     recursive, so a subdirectory is out of reach, and next to the workflows it
     documents is where it belonged anyway. -->  

# CI: what runs, and why it is built this way

All four workflows are live. Each also keeps its `workflow_dispatch` trigger, so
any of them can be run on demand from the Actions tab (→ the workflow → *Run
workflow*) without waiting for a matching push.

| Workflow | When it runs | What it does |
| --- | --- | --- |
| `R-CMD-check.yaml` | every push and PR, except docs-only commits; superseded runs cancelled | a 5-runner matrix, plus one Ubuntu job checking the generated model code is current |
| `pkgdown.yaml` | every push, PR and release | builds and deploys the site to `gh-pages` |

`dependabot.yml` is active and needs no switch.

The Monday run is the one with no substitute: it is the only thing that catches
`malariasimulation` changing underneath us, which is the single staleness cause
that never appears in our own commits. GitHub emails the repository owner when a
scheduled run fails.

## Why it is built this way

**The compiled model can go stale while everything stays green.**
`inst/odin/malaria_ode.R` is the model's source of truth, but the artifact that
actually gets compiled is the committed, generated `src/malaria_ode.cpp` — along
with `R/dust.R`, `R/cpp11.R`, `src/cpp11.cpp` and `inst/dust/`. No ordinary build
regenerates any of them. So: edit the odin model, forget to run
`odin2::odin_package(".")`, and all five check runners compile the *old* model
and pass, the tests pass, and `check_drift.R` passes most emphatically of all —
"nothing moved" reads as *my change was numerically inert* when what it really
means is *my change was never compiled*. The `generated code matches inst/odin`
job in `R-CMD-check.yaml` closes that gap: it regenerates from `inst/odin` and
fails if the working tree comes out dirty. It runs once, on Ubuntu, rather than
as a sixth matrix row, because the generated code comes out the same whatever
platform writes it and the package itself is never built. It also goes red when
odin2 or dust2 themselves move, since the generator's version is stamped into the
header of every file it writes and the `Remotes:` pins are to HEAD rather than to
a tag; the fix is the same either way — regenerate and commit.

**And its dependency set is hard dependencies only.** The job calls one
generator, so it installs `Imports`/`LinkingTo` plus `odin2` and the few
generator helpers named in the workflow — not `dependencies: "all"`, which would
drag in the whole `Suggests` tree: `malariasimulation` and `postie`, both GitHub
remotes that have to compile from source, plus `knitr`, `rmarkdown` and
`testthat`. This job loads none of them. With them installed, a
`malariasimulation` build failure on a runner turns the staleness guard red for a
reason that has nothing to do with staleness — which is exactly the failure this
guard cannot afford, because its whole value is that a red means one specific
thing.

**Two entries in that list look wrong and are not.** `odin2` is named as the
GitHub ref `mrc-ide/odin2`, not as `any::odin2`, because `any::` is a
CRAN-style ref: pkgdepends resolves it from repository metadata alone and never
consults GitHub or `Remotes:`. odin2 is not on CRAN, so `any::odin2` fails at
resolution ("Can't find package called any::odin2") and the job dies during the
dependency install. `withr` is on CRAN, so `any::` is right for it, but nothing
visible asks for it: `odin_package()` ends in `cpp11::cpp_register()`, whose
`get_call_entries()` wraps its scan in `withr::with_collate("C", ...)` to make
the generated `CallEntries` table locale-independent. withr is only a `Suggests`
of cpp11 and is not in cpp11's own `Config/Needs/cpp11/cpp_register` list, so
nothing in this job's hard-dependency closure installs it, and without it the
job writes `cpp11.R` and then fails with "there is no package called 'withr'".

**And the guard reads `git status`, not `git diff`.** `git diff` compares
against the index, so a generated file that has never been committed is
untracked and invisible to it: the guard would pass on a tree that is missing a
file the build needs. `git status --porcelain` sees both. It stays scoped to
`src`, `R/dust.R`, `R/cpp11.R` and `inst/dust`, because `setup-r-dependencies`
leaves an untracked `.github/pkg.lock` in the tree and an unscoped check would
go red on every run.

**A missing `malariasimulation` now fails the check instead of quietly emptying
it.** Nearly every `test_that()` block opens with
`skip_if_not_installed("malariasimulation")` — it is a Suggests dependency and
the source of every parameter list the model is driven with. So a runner where it
was missing or unloadable used to run the handful of dependency-free tests, skip
the rest, and report `OK, 0 errors`: a green tick meaning "we tested three
things". `tests/testthat/setup.R` now `stop()`s in that case, which changes how
the matrix fails. Before, a runner that could not build or load the IBM from the
`Remotes:` pin passed alongside the four that could, and nothing said the
difference. Now that row goes red at the test step, naming the missing package,
while the others stay green — so a one-platform dependency problem reads as one,
instead of as a model that passed everywhere. Nothing else in the check would
have said so:
`vignettes/fleet.Rmd` gates every chunk on `malariasimulation` and `postie` being
installed, so it knits to prose and passes with no code run at all, and the other
two vignettes are `eval = FALSE` throughout. The test suite is the only part of
`R CMD check` that notices.

**`paths-ignore`, not `paths`, wherever the question is "could this break the
model".** A deny-list of things that provably cannot (markdown, PNGs, the pkgdown
site) stays correct when the package grows a new source directory. An allow-list
would quietly stop checking it, and nothing would say so.

**Dependabot handles Actions deprecations.** Runner images and actions are retired
on GitHub's schedule, not ours; the failure mode without it is a red build one
morning for a reason that has nothing to do with the package. Monthly and grouped
into a single PR, because a PR that gets ignored is worse than one that gets read.

**It is hard dependencies only as well.** The two render scripts read committed
CSVs and JSON and draw pictures. Neither loads `malariasimulation` or `postie`,
so on the action's default of `"all"` this job built both from source for
nothing: a long install, and two more ways for it to go red without a figure
having moved. Neither loads fleet either, which is why the job never installs
the package under test.

## Two committed baselines of fleet's own numbers

There are two, they are regenerated by different commands, and a deliberate model
change moves **both**.

| Baseline | Pinned in | Checked by | Regenerate with |
| --- | --- | --- | --- |
| Unit-test reference | `tests/testthat/reference-values.csv` | `tests/testthat/test-reference.R`, in every check run | `FLEET_REGENERATE_REFERENCE=1 Rscript -e 'devtools::test(filter = "reference")'` |

## Environment variables

The two the test suite reads. Neither should ever be set in CI.

| Variable | Read by | What it does |
| --- | --- | --- |
| `FLEET_ALLOW_SKIP` | `tests/testthat/setup.R` | Any non-empty value restores the old behaviour: a missing `malariasimulation` skips the suite instead of erroring. For a deliberate local run on a machine where the IBM genuinely cannot be built, where you want the handful of dependency-free tests. Setting it in CI re-opens exactly the hole the `stop()` closes. |
| `FLEET_REGENERATE_REFERENCE` | `tests/testthat/test-reference.R` | Any non-empty value rewrites `tests/testthat/reference-values.csv` and `reference-interventions.csv` from the current model, then passes trivially. Only for a model change that was intended and reviewed: read the diff, it is the record of what the change did. Never to turn a red test green. |

## A note on Actions minutes

Private repositories consume the account's Actions quota; public ones get
unlimited standard-runner minutes. The full `R-CMD-check` matrix costs roughly 87
billed minutes per run — macOS is charged at 10x and is 56% of that, despite being
the fastest job on the wall clock. The generated-code job, and the drift and
figures workflows, are Ubuntu-only and cost a few minutes each. The generated-code
job is cheap only because of its `dependencies: '"hard"'`, and the figures job
for the same reason: what they spend is a dependency install, so widening either
set is the one edit that would make them expensive as well as flaky. The drift
check is the opposite case and its default `'"all"'` is deliberate: `Suggests` is
where `malariasimulation` lives, and narrowing that job would leave
`check_drift.R` with nothing to compare against.

If this repository stays private, the cheapest further saving by far is moving the
macOS row off every-push and onto the weekly schedule. If it goes public, none of
this matters and the matrix can stay as it is.
