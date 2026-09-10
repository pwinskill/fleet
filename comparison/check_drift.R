# Has a change moved blink, and is it still matching the IBM?
#
#   Rscript comparison/check_drift.R              # ~2 min
#   CMP_ONLY=eir_20,smc Rscript comparison/check_drift.R    # a subset, faster
#   CMP_STRICT=1 Rscript comparison/check_drift.R # also fail if ANY number moved
#
# The IBM does not depend on blink, so its committed rows stay valid for any
# blink-side change and there is no reason to re-run a 25-minute IBM sweep to
# find out whether the match still holds. This re-runs blink only, against the
# frozen reference.
#
# Two questions, reported separately because they mean different things:
#
#   1. Did anything MOVE?  blink now against blink's committed rows. A moved
#      number is not automatically wrong -- a deliberate model fix moves numbers
#      -- but it must be SEEN. Silent movement is how a regression ships.
#   2. Is the match still GOOD?  blink now against the committed IBM medians and
#      10-90% bands, at the thresholds the article's claims rest on. This is the
#      gate, and the only thing that fails the run by default.
#
# Exit 0 = still matching. Exit 1 = drift beyond threshold (or, under
# CMP_STRICT, any movement at all).

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
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
source(file.path(ROOT, "comparison", "constants.R"))
source(file.path(ROOT, "comparison", "scenarios.R"))
DDIR <- file.path(ROOT, "comparison", "data")
STRICT <- nzchar(Sys.getenv("CMP_STRICT"))

## ---- thresholds ---------------------------------------------------------------
## `rel` is the largest tolerated |blink / IBM median - 1| over the EIR grid, set
## from what the models achieve today with room to spare (currently 1.6%, 2.6%,
## 2.1% and 6.2% against the limits below).
##
## `out` is how many grid points may fall OUTSIDE the IBM's 10-90% replicate
## range -- a count of failures, not of successes, so CMP_ONLY can check a subset
## without the limit becoming unreachable. These are set AT today's values, with
## no slack: 0, 1, 1 and 2 of 6. That is deliberate rather than an oversight.
## blink is deterministic and the IBM reference is frozen, so nothing here
## fluctuates; a point moving in or out of the band is a real change every time,
## and there is no noise for headroom to absorb. Loosen one only when you have
## decided the new value is correct, and say so in the commit.
##
## Severe is loosest on both because it is the rarest outcome and the noisiest in
## the IBM: see the article's Programmes section for the widths.
LIMITS <- list(
  pfpr_2_10 = list(rel = 0.02, out = 0L, label = "PfPR 2-10"),
  clin_0_5  = list(rel = 0.05, out = 1L, label = "clinical, 0-5"),
  clin_all  = list(rel = 0.05, out = 1L, label = "clinical, all ages"),
  sev_all   = list(rel = 0.10, out = 2L, label = "severe, all ages"))
## anything moving by more than this against the committed blink rows is reported;
## solver output is deterministic, so this is a floating-point floor, not a budget
MOVE_TOL <- 1e-6

fail <- character(); warn <- character()
rule <- function(t) cat("\n", t, "\n", strrep("-", nchar(t)), "\n", sep = "")

## ---- 0. is the reference still describing these scenarios? ---------------------
rule("Reference")
ref_f <- file.path(DDIR, "ibm_reference.json")
if (!file.exists(ref_f)) {
  warn <- c(warn, "no ibm_reference.json: the committed IBM rows have no provenance. Re-run comparison/run_replicates.R to stamp one.")
  cat("  (absent)\n")
} else {
  ref <- jsonlite::read_json(ref_f, simplifyVector = TRUE)
  cat(sprintf("  IBM rows generated %s\n", ref$generated))
  cat(sprintf("  malariasimulation  %s (installed: %s)\n", ref$malariasimulation,
              as.character(utils::packageVersion("malariasimulation"))))
  cat(sprintf("  %d replicates of %s people, %d-year burn-in\n",
              ref$n_rep, format(ref$population, big.mark = ","), ref$burn_in_years))
  ## a FAIL, not a warning, for the same reason as the digest below: the whole
  ## design rests on the IBM rows staying valid while blink changes. A different
  ## malariasimulation is a different IBM, so section 2 would be comparing blink
  ## against a reference the current upstream would no longer produce -- and
  ## `blink` exists to be a twin of one specific version. This is the condition
  ## the weekly run exists to catch, so it must be loud enough to stop the build.
  if (!identical(ref$malariasimulation, as.character(utils::packageVersion("malariasimulation"))))
    fail <- c(fail, sprintf("malariasimulation has changed since the IBM rows were made (%s -> %s), so the agreement below compares blink against a reference the installed IBM would no longer reproduce. Re-run comparison/run_replicates.R (without CMP_BLINK_ONLY) to rebuild the IBM rows.",
                            ref$malariasimulation, utils::packageVersion("malariasimulation")))
  now <- scenario_digest()
  if (!identical(ref$scenario_digest, now)) {
    ## a FAIL, not a warning: with the scenarios changed, section 2 below is
    ## comparing blink-on-new-scenarios against IBM-on-old-scenarios, which is
    ## not a valid answer to "is the match still good" no matter what it prints
    fail <- c(fail, sprintf("the scenario definitions have changed since the IBM rows were made (digest %s -> %s), so the agreement below compares blink on the new scenarios against the IBM on the old ones. Re-run comparison/run_replicates.R.",
                            ref$scenario_digest, now))
    cat("  scenario digest    CHANGED\n")
  } else {
    cat(sprintf("  scenario digest    %s (unchanged)\n", now))
  }
}

## ---- run blink against the frozen reference -----------------------------------
rule("Re-running blink")
new_eq <- do.call(rbind, lapply(run_blink(), `[[`, "eq"))
old <- read.csv(file.path(DDIR, "rep_eq.csv"), stringsAsFactors = FALSE)
MET <- names(LIMITS)

## ---- 1. did anything move? ----------------------------------------------------
rule("1. Movement against the committed blink rows")
ob <- old[old$model == "blink", ]
cmp <- merge(new_eq[, c("scenario", MET)], ob[, c("scenario", MET)],
             by = "scenario", suffixes = c(".new", ".old"))
moved <- do.call(rbind, lapply(MET, function(m) {
  a <- cmp[[paste0(m, ".new")]]; b <- cmp[[paste0(m, ".old")]]
  rel <- ifelse(abs(b) > 0, abs(a / b - 1), NA_real_)
  k <- which(!is.na(rel) & rel > MOVE_TOL)
  if (!length(k)) NULL else
    data.frame(scenario = cmp$scenario[k], outcome = LIMITS[[m]]$label,
               committed = b[k], now = a[k], rel = rel[k])
}))
if (is.null(moved)) {
  cat("  nothing moved: every value reproduces the committed rows to 1e-6.\n")
} else {
  moved <- moved[order(-moved$rel), ]
  cat(sprintf("  %d of %d scenario x outcome values moved:\n\n", nrow(moved),
              nrow(cmp) * length(MET)))
  cat(sprintf("    %-12s %-20s %12s %12s %9s\n", "scenario", "outcome", "committed", "now", "change"))
  for (i in seq_len(min(nrow(moved), 25))) with(moved[i, ], cat(sprintf(
    "    %-12s %-20s %12.5g %12.5g %+8.2f%%\n", scenario, outcome, committed, now, 100 * rel)))
  if (nrow(moved) > 25) cat(sprintf("    ... and %d more\n", nrow(moved) - 25))
  cat("\n  Movement is not automatically wrong. If it was intended, refresh the\n",
      "  committed rows with CMP_BLINK_ONLY=1 Rscript comparison/run_replicates.R\n", sep = "")
  if (STRICT) fail <- c(fail, sprintf("%d values moved (CMP_STRICT)", nrow(moved)))
}

## ---- 2. is the match still good? ----------------------------------------------
rule("2. Agreement with the committed IBM reference")
e <- old[old$model == "IBM" & grepl("^eir_", old$scenario), ]
if (!nrow(e)) {
  warn <- c(warn, "no IBM rows for the EIR grid in rep_eq.csv; cannot check agreement.")
} else {
  cat(sprintf("    %-20s %10s %8s %15s %7s\n", "outcome", "max |rel|", "limit", "outside band", "limit"))
  for (m in MET) {
    lim <- LIMITS[[m]]
    st <- do.call(rbind, lapply(split(e, e$scenario), function(g) data.frame(
      scenario = g$scenario[1], med = median(g[[m]]),
      lo = quantile(g[[m]], .1), hi = quantile(g[[m]], .9))))
    st <- merge(st, new_eq[, c("scenario", m)], by = "scenario")
    v <- st[[m]]
    rel <- max(abs(v / st$med - 1))
    outb <- sum(v < st$lo | v > st$hi)
    bad <- rel > lim$rel || outb > lim$out
    cat(sprintf("    %-20s %9.1f%% %7.1f%% %9d of %-3d %6d  %s\n", lim$label,
                100 * rel, 100 * lim$rel, outb, nrow(st), lim$out,
                if (bad) "<-- DRIFT" else "ok"))
    if (bad) fail <- c(fail, sprintf("%s: max |rel| %.1f%% (limit %.1f%%), outside band %d of %d (limit %d)",
                                     lim$label, 100 * rel, 100 * lim$rel, outb, nrow(st), lim$out))
  }
}

## ---- verdict -------------------------------------------------------------------
rule("Verdict")
for (w in warn) cat("  WARNING: ", w, "\n", sep = "")
if (length(fail)) {
  for (f in fail) cat("  FAIL: ", f, "\n", sep = "")
  cat("\n  Drifted. Either the change is wrong, or the reference needs refreshing.\n")
  quit(status = 1)
}
cat("  Still matching.", if (length(warn)) " (with warnings above)" else "", "\n", sep = "")
