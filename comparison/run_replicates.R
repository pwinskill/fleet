# Run every comparison scenario through BOTH models and save tidy CSVs.
#
#   Rscript comparison/run_replicates.R          # ~25 min on 10 workers (see cost note)
#   Rscript comparison/render_figures.R          # seconds: figures from the CSVs
#   CMP_ONLY=nets,smc Rscript comparison/run_replicates.R   # re-run a subset, merge into the CSVs
#
# The IBM (malariasimulation) is run N_REP times per scenario with different
# seeds, in parallel on a PSOCK cluster, and summarised per replicate; blink is
# run once (it is deterministic). Both get the SAME parameter list, with the
# rendering bands set once so their output columns line up exactly.
#
# Every scenario is burned in BURN_Y years in the IBM so it reaches its own
# stochastic steady state; blink is seeded at the malariaEquilibrium fixed point
# and integrated over the same horizon so both are on a common clock.
#
# Cost note: the IBM's per-band rendering dominates its run time at 10k people
# (~2.5 s per simulated year with the three default ranges, ~5 s with the 12-band
# age profile added to every family), so only the reference scenario carries the
# age-profile bands. Jobs are load-balanced, longest first.

## No absolute paths anywhere in here. BLINK_LIB prepends an R library, for
## installations that do not pick up R_LIBS_USER (the Windows-arm64 setup this was
## developed on); leave it unset and your normal library is used. ROOT is found by
## walking up to the DESCRIPTION, so these scripts run from any working directory
## and on anyone's checkout, whether via Rscript or source().
if (nzchar(.l <- Sys.getenv("BLINK_LIB"))) .libPaths(.l)
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(dirname(.f), "/") else getwd()
while (!file.exists(file.path(ROOT, "DESCRIPTION")) && dirname(ROOT) != ROOT)
  ROOT <- dirname(ROOT)
if (!file.exists(file.path(ROOT, "DESCRIPTION")))
  stop("run this from inside the blink checkout (no DESCRIPTION found above ", getwd(), ")")
suppressMessages(library(malariasimulation))
source(file.path(ROOT, "comparison", "theme.R"))     # scenario constants
DDIR <- file.path(ROOT, "comparison", "data"); dir.create(DDIR, showWarnings = FALSE)
N_WORKERS <- 10L
## CMP_SMOKE=1 -> a few-minute end-to-end check: 4-year horizon, one replicate,
## interventions at year 1 so every builder actually fires, output to data/smoke/.
SMOKE <- nzchar(Sys.getenv("CMP_SMOKE"))
if (SMOKE) { N_REP <- 1L; N_WORKERS <- 3L; BURN_Y <- 1L
  DDIR <- file.path(DDIR, "smoke"); dir.create(DDIR, showWarnings = FALSE) }
PROG <- file.path(DDIR, "progress.log"); if (file.exists(PROG)) invisible(file.remove(PROG))
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

## Scenario definitions and the shared summariser. Kept in their own file so the
## drift check can reach them without starting a run.
source(file.path(ROOT, "comparison", "scenarios.R"))
## CMP_BLINK_ONLY=1 -> re-run only blink (the IBM rows are kept): for a blink-side
## model change, when the IBM results are unaffected
BLINK_ONLY <- nzchar(Sys.getenv("CMP_BLINK_ONLY"))


## ---- run blink (deterministic, seconds) --------------------------------------
ode <- run_blink()

## ---- run the IBM replicates in parallel --------------------------------------
ibm <- list()
if (!BLINK_ONLY) {
jobs <- expand.grid(scenario = names(scenarios), rep = seq_len(N_REP),
                    stringsAsFactors = FALSE)
## rough cost in default-band sim-years, so the load balancer starts the longest jobs first
jobs$cost <- vapply(jobs$scenario, function(nm) scenarios[[nm]]$years *
  (if (length(scenarios[[nm]]$p$prevalence_rendering_min_ages) > 1) 2 else 1), numeric(1))
jobs <- jobs[order(-jobs$cost, jobs$rep), ]; rownames(jobs) <- NULL
log_msg("IBM: %d runs (%d scenarios x %d reps) on %d workers; ~%.0f min at 2.5 s per sim-year",
        nrow(jobs), length(scenarios), N_REP, N_WORKERS, sum(jobs$cost) * 2.5 / N_WORKERS / 60)
cl <- parallel::makeCluster(N_WORKERS)
## workers do not inherit .libPaths(), so hand them the parent's rather than
## hardcoding one: whatever library this session is using, they use too
.libs <- .libPaths()
parallel::clusterExport(cl, c("jobs", "scenarios", "summarise_run", "tag_parts", "band_lo",
                              "band_hi", "AGE_TAGS", "POP", "BURN_Y", "PROG", ".libs"))
invisible(parallel::clusterEvalQ(cl, {
  .libPaths(.libs)
  suppressMessages(library(malariasimulation))
}))
t0 <- Sys.time()
ibm <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(j) {
  nm <- jobs$scenario[j]; k <- jobs$rep[j]; s <- scenarios[[nm]]
  set.seed(1000L + 7L * k)
  el <- system.time(out <- run_simulation(timesteps = s$years * 365, parameters = s$p))[["elapsed"]]
  r <- summarise_run(out, s$years)
  r$timing <- data.frame(years = s$years, elapsed_s = el)
  cat(sprintf("[%s] %-10s rep %d  %5.1f min\n", format(Sys.time(), "%H:%M:%S"), nm, k, el / 60),
      file = PROG, append = TRUE)
  tag_parts(r, nm, "IBM", k)
}, chunk.size = 1L)
parallel::stopCluster(cl)
log_msg("IBM done in %.1f min", as.numeric(Sys.time() - t0, units = "mins"))
}

## ---- combine and write -------------------------------------------------------
bind <- function(part) do.call(rbind, c(lapply(ode, `[[`, part), lapply(ibm, `[[`, part)))
for (part in c("eq", "age", "monthly", "doy", "timing")) {
  d <- bind(part); f <- file.path(DDIR, paste0("rep_", part, ".csv"))
  if ((length(ONLY) || BLINK_ONLY) && file.exists(f)) {   # partial run: replace just those rows
    old <- read.csv(f, stringsAsFactors = FALSE)
    drop <- old$scenario %in% names(scenarios) & (if (BLINK_ONLY) old$model == "blink" else TRUE)
    old <- old[!drop, ]
    for (nm in setdiff(names(d), names(old))) old[[nm]] <- NA
    for (nm in setdiff(names(old), names(d))) d[[nm]] <- NA
    d <- rbind(old[names(d)], d)
  }
  write.csv(d, f, row.names = FALSE)
  log_msg("wrote rep_%s.csv (%d rows)", part, nrow(d))
}
## the IBM rows just changed, so stamp what produced them. Skipped under
## BLINK_ONLY and CMP_ONLY, where the committed IBM rows are only partly
## refreshed and the existing stamp still describes the rest.
if (!BLINK_ONLY && !length(ONLY) && !SMOKE) {
  jsonlite::write_json(ibm_reference(), file.path(DDIR, "ibm_reference.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  log_msg("wrote ibm_reference.json (scenario digest %s)", scenario_digest())
}
log_msg("ALL DONE")
