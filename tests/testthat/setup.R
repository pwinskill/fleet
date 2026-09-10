# Fail loudly when the suite cannot actually test anything.
#
# All but a handful of the test_that() blocks here open with skip_if_not_installed
# ("malariasimulation") -- it is a Suggests dependency and the source of every
# parameter list the model is driven with. If it fails to install, every one of
# those blocks skips and the run still reports "OK, 0 errors" having exercised
# almost nothing. A green tick that means "we tested three things" is worse than
# a red one.
#
# So: absent malariasimulation is an ERROR here, not a skip.
#
# Escape hatch: set the environment variable BLINK_ALLOW_SKIP to any non-empty
# value to fall back to the old skip-everything behaviour -- for a deliberate run
# on a machine where malariasimulation genuinely cannot be built, where you want
# the handful of dependency-free tests to run. Never set it in CI.
if (!requireNamespace("malariasimulation", quietly = TRUE) &&
    !nzchar(Sys.getenv("BLINK_ALLOW_SKIP"))) {
  stop("'malariasimulation' is not installed, so nearly every test in this suite ",
       "would skip and the run would report success having tested almost nothing.\n",
       "Install it (Remotes: mrc-ide/malariasimulation), or set BLINK_ALLOW_SKIP=1 ",
       "to accept a near-empty run.", call. = FALSE)
}
