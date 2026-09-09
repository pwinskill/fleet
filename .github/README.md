# CI: what runs, what does not yet, and how to turn it on

Two of the four workflows are **prepared but inert**. They carry only a
`workflow_dispatch` trigger, so they never fire on their own; their real triggers
sit directly above, commented, one uncomment away. That way they can be tested on
demand (Actions tab → the workflow → *Run workflow*) before anything becomes
automatic.

| Workflow | Runs now | When enabled |
| --- | --- | --- |
| `R-CMD-check.yaml` | every push and PR, 5-runner matrix | same, minus docs-only commits, and superseded runs cancelled |
| `pkgdown.yaml` | every push, PR and release | unchanged |
| `comparison.yaml` | **manual only** | every push touching the model, plus Mondays |
| `figures.yaml` | **manual only** | pushes touching the comparison data or renderers |

## Going live

Three edits, all uncommenting a block that is already written and marked
`TO ENABLE`:

1. **`comparison.yaml`** — uncomment `push`, `pull_request` and `schedule`.
   This is the one that matters: it means a change to the model cannot land
   without the blink-vs-IBM match being re-checked.
2. **`figures.yaml`** — uncomment `push` and `pull_request`.
3. **`R-CMD-check.yaml`** — uncomment the two `paths-ignore` blocks and the
   `concurrency` block.

Nothing else needs configuring. `dependabot.yml` is already active on merge and
needs no switch.

## Why it is built this way

**`paths-ignore`, not `paths`, wherever the question is "could this break the
model".** A deny-list of things that provably cannot (markdown, PNGs, the pkgdown
site) stays correct when the package grows a new source directory. An allow-list
would quietly stop checking it, and nothing would say so. `figures.yaml` is the
exception and uses `paths`, because there the question is the narrow one of
whether a specific set of inputs changed.

**The drift check does not re-run the IBM.** The IBM does not depend on blink, so
its committed rows stay valid for any blink-side change. `check_drift.R` re-runs
blink alone against that frozen reference: about two minutes, against
twenty-five for a full IBM sweep. That is the whole reason this can run on every
push.

**Its dependency list is deliberately short.** `check_drift.R` sources
`comparison/constants.R` rather than `theme.R`, so it needs no plotting stack:
blink, malariasimulation, digest, jsonlite. Every package in that list is one
more thing that can break the check for a reason unrelated to the model.

**The weekly run is the only thing that catches malariasimulation moving.** It is
the one staleness cause that never appears in our own commits. `check_drift.R`
compares the installed version and a digest of the scenario definitions against
`comparison/data/ibm_reference.json`, and fails when the scenarios have changed —
because at that point it would be comparing blink on the new scenarios against
the IBM on the old ones, which is not an answer to the question no matter what it
prints. GitHub emails the repo owner when a scheduled run fails.

**Dependabot handles Actions deprecations.** Runner images and actions are retired
on GitHub's schedule, not ours; the failure mode without it is a red build one
morning for a reason that has nothing to do with the package. Monthly and grouped
into a single PR, because a PR that gets ignored is worse than one that gets read.

**`figures.yaml` fails on the tables, warns on the images.** A PNG can differ
byte-for-byte between machines from font hinting alone, so failing on images would
produce noise that trains you to ignore it. `tables.md` is plain text computed
from the same CSVs, so a genuinely stale figure almost always shows up there too.

## A note on Actions minutes

Private repositories consume the account's Actions quota; public ones get
unlimited standard-runner minutes. The full `R-CMD-check` matrix costs roughly 87
billed minutes per run — macOS is charged at 10x and is 56% of that, despite being
the fastest job on the wall clock. The drift and figures workflows are Ubuntu-only
and cost a few minutes each.

If this repository stays private, the cheapest further saving by far is moving the
macOS row off every-push and onto the weekly schedule. If it goes public, none of
this matters and the matrix can stay as it is.
