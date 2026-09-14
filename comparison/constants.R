# Constants shared by every comparison script: paths, and the scenario settings
# that must be identical across them.
#
# Split out of theme.R so the scripts that do not draw anything -- the 25-minute
# IBM run and the drift check -- do not have to load ggplot2, patchwork, ragg and
# systemfonts just to read POP and EIR_GRID. That keeps the drift check's CI
# dependency list to fleet, malariasimulation, digest and jsonlite, which is both
# faster and one less thing to break when a plotting package changes.
#
# theme.R sources this, so anything that draws still gets these too.

## ---- paths ------------------------------------------------------------------
## Every comparison script sets ROOT (the checkout root) before sourcing this
## file, so nothing here is machine-specific either. VDIR points at the optional
## site-file validation results, which live in their own checkout beside this
## one; override with FLEET_VALIDATE. The site figure and its table are skipped
## when it is absent.
if (!exists("ROOT") || !file.exists(file.path(ROOT, "DESCRIPTION")))
  stop("source this from a comparison/ script, which sets ROOT to the fleet checkout.")
VDIR <- Sys.getenv("FLEET_VALIDATE", file.path(dirname(ROOT), "fleet_validate"))

## ---- shared scenario constants (must match run_replicates.R) ----------------
BURN_Y   <- 30L                 # IBM burn-in years before observation / intervention
POP      <- 10000L              # IBM population
N_REP    <- 10L                  # IBM replicates per scenario
AGE_EDGES <- c(0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 85)   # age-profile bands, years
SEASON   <- list(g0 = 0.285, g = c(-0.33, -0.13, 0.052), h = c(-0.35, 0.020, 0.10))
EIR_GRID <- c(1, 3, 10, 20, 50, 120)
EIR_REF  <- 20                  # reference transmission for age profiles + interventions

## intervention scenario labels, in display order (highest expected impact last)
INT_LABELS <- c(
  treatment = "Treatment scale-up\n20% \u2192 60% of clinical cases",
  pev       = "RTS,S via EPI\n90% at 5 months, booster",
  smc       = "Seasonal SMC\n4 rounds/yr, ages 0.25\u20135",
  irs       = "Indoor residual spraying\n80% coverage, annual",
  nets      = "Bed-net campaign\n80% coverage, one round")

## long-horizon programme scenarios (ts_*), in tier order: nothing, one thing,
## everything. All EIR 20 seasonal with 20% baseline case management, so each row
## is the row above it with one more thing added.
## kept short on purpose: these are row strips down the left of a 4-column
## figure, and every character of label is width taken from the panels
TS_LABELS <- c(
  ts_none  = "No\ninterventions",
  ts_nets  = "Bed nets\n5 campaigns",
  ts_smc   = "Seasonal\nSMC",
  ts_treat = "Case\nmanagement",
  ts_all   = "All three\ntogether")
TS_YEARS <- 15L                 # years followed past deployment
TS_NET_EVERY <- 3L              # years between net campaigns
