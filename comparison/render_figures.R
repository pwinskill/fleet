# Build every comparison figure from the saved CSVs (no model runs).
#
#   Rscript comparison/render_figures.R                 # figures -> man/figures + vignettes
#   CMP_SMOKE=1 Rscript comparison/render_figures.R     # from data/smoke -> comparison/plots/smoke
#
# Writes cmp_*.png. All visual decisions live in theme.R; this file only shapes
# data and composes panels. Figures:
#   core_eir       four equilibrium relationships vs EIR (2x2)
#   core_age       age profiles at EIR 20: prevalence, clinical, severe
#   core_seasonal  the settled annual cycle: prevalence and clinical incidence
#   core_sites     63-country monthly comparison (from fleet_validate results)
#   int_timeseries five interventions x {prevalence, clinical}, time series
#   programme_ts   five long-horizon programmes x four outcomes, 15 years
#   int_impact     % reduction per intervention, IBM (replicate range) vs fleet

## No absolute paths anywhere in here. FLEET_LIB is prepended to the library
## path, for installations that do not pick up R_LIBS_USER (the Windows-arm64
## setup this was developed on); the libraries already on the path are kept, so a
## FLEET_LIB holding only some of the dependencies still works. Leave it unset and
## your normal library is used. ROOT is found by walking up to the DESCRIPTION, so
## these scripts run from any working directory and on anyone's checkout, whether
## via Rscript or source().
if (nzchar(.l <- Sys.getenv("FLEET_LIB"))) .libPaths(c(.l, .libPaths()))
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(dirname(.f), "/") else getwd()
while (!file.exists(file.path(ROOT, "DESCRIPTION")) && dirname(ROOT) != ROOT)
  ROOT <- dirname(ROOT)
if (!file.exists(file.path(ROOT, "DESCRIPTION")))
  stop("run this from inside the fleet checkout (no DESCRIPTION found above ", getwd(), ")")
suppressMessages({library(dplyr); library(tidyr)})
source(file.path(ROOT, "comparison", "theme.R"))
SMOKE <- nzchar(Sys.getenv("CMP_SMOKE"))
DDIR  <- file.path(ROOT, "comparison", "data")
if (SMOKE) {                                   # smoke data must never overwrite the real figures
  BURN_Y <- 1L; DDIR <- file.path(DDIR, "smoke")
  out_dir <- file.path(ROOT, "comparison", "plots", "smoke")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_fig <- function(g, name, width, height, dpi = 200)
    ggsave(file.path(out_dir, paste0("cmp_", name, ".png")), g, width = width,
           height = height, dpi = dpi, device = ragg::agg_png, bg = SURFACE)
}
rd <- function(part) read.csv(file.path(DDIR, paste0("rep_", part, ".csv")), stringsAsFactors = FALSE)
eq <- rd("eq"); age <- rd("age"); monthly <- rd("monthly"); doy <- rd("doy"); timing <- rd("timing")
n_rep <- max(eq$rep)
ibm_note <- sprintf(paste(
  "IBM: %d stochastic replicates of %s people, %d-year burn-in; line/point = median, band/bar = 10\u201390%% range across replicates.",
  "fleet: one deterministic run seeded at equilibrium."), n_rep, format(POP, big.mark = ","), BURN_Y)
PREV_LAB <- "LM prevalence, ages 2\u201310"
CLIN_LAB <- "clinical episodes per child-year, ages 0\u20135"

## ============================================================================
## 1. core_eir -- four equilibrium relationships vs EIR
## ============================================================================
## no panel titles here: the y axis already names the outcome and its age band,
## so a title on top of it can only restate it or editorialise. The shape of each
## curve is what the panel is for, and the article says what to make of it.
EIR_MET <- list(
  list(key = "pfpr_2_10", lab = PREV_LAB, pct = TRUE),
  list(key = "clin_0_5", lab = CLIN_LAB),
  list(key = "clin_all", lab = "clinical episodes per person-year, all ages"),
  list(key = "sev_all", lab = "severe episodes per 1,000 person-years, all ages"))
e_long <- eq %>% filter(grepl("^eir_", scenario)) %>%
  mutate(init_EIR = as.numeric(sub("eir_", "", scenario))) %>%
  select(init_EIR, model, rep, all_of(vapply(EIR_MET, `[[`, "", "key"))) %>%
  pivot_longer(-c(init_EIR, model, rep), names_to = "metric", values_to = "y")
ibm_e <- e_long %>% filter(model == "IBM") %>% group_by(metric) %>%
  group_modify(~ envelope(.x, by = "init_EIR")) %>% ungroup()
ode_e <- e_long %>% filter(model == "fleet") %>% rename(mid = y)

panel_eir <- function(metric_id, ylab) {
  gi <- filter(ibm_e, metric == metric_id); go <- filter(ode_e, metric == metric_id)
  ggplot() +
    geom_linerange(data = gi, aes(init_EIR, ymin = lo, ymax = hi, colour = model),
                   linewidth = 0.7, alpha = 0.55, show.legend = FALSE) +
    ## dashed over solid, hollow over filled -- see the layer-order rule in theme.R
    geom_line(data = go, aes(init_EIR, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = gi, aes(init_EIR, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_point(data = gi, aes(init_EIR, mid, colour = model, shape = model, fill = model),
               size = 2.6, stroke = 0.5) +
    geom_point(data = go, aes(init_EIR, mid, colour = model, shape = model, fill = model),
               size = 2.6, stroke = 0.5) +
    scale_x_log10(breaks = EIR_GRID, minor_breaks = NULL) +
    scale_models() + guide_models() +
    labs(x = "EIR (infectious bites per adult per year)", y = ylab) +
    theme_cmp()
}
ps <- lapply(seq_along(EIR_MET), function(i) {
  m <- EIR_MET[[i]]
  p <- panel_eir(m$key, m$lab) +
    if (isTRUE(m$pct)) scale_y_continuous(limits = c(0, 1), labels = scales::percent)
    else scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)))
  ## the x axis is the same in all four; name it once, on the bottom row
  if (i <= 2) p + labs(x = NULL) else p
})
## one collected legend for the whole figure, not a legend sitting inside panel 1:
## with no panel titles left to anchor the eye, a legend in one panel pushes that
## panel's plot area down and the top row stops lining up
g <- patchwork::wrap_plots(ps, ncol = 2) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Core transmission relationships at equilibrium",
    subtitle = "The same parameter list through both models, across a 120-fold range of transmission intensity",
    caption = cap("x = the EIR passed to set_equilibrium(); each model's realised EIR is reported in the article.", ibm_note),
    theme = theme_cmp()) &
  theme(legend.position = "top", legend.justification = "left")
save_fig(g, "core_eir", width = 10, height = 9)

## ============================================================================
## 2. core_age -- age profiles at the reference EIR
## ============================================================================
a <- age %>% filter(scenario == paste0("eir_", EIR_REF)) %>%
  pivot_longer(c(prev, clin, sev), names_to = "metric", values_to = "y")
ibm_a <- a %>% filter(model == "IBM") %>% group_by(metric) %>%
  group_modify(~ envelope(.x, by = "age_mid")) %>% ungroup()
ode_a <- a %>% filter(model == "fleet") %>% rename(mid = y)
lab_a <- c(prev = "LM prevalence", clin = "clinical episodes per person-year",
           sev = "severe episodes per 1,000 person-years")
ttl_a <- c(prev = "Prevalence peaks in\nschool-age children",
           clin = "Clinical disease\nconcentrates in the young",
           sev  = "Severe disease is rarer\nand earlier still")
panel_age <- function(m) {
  gi <- filter(ibm_a, metric == m); go <- filter(ode_a, metric == m)
  ggplot() +
    geom_ribbon(data = gi, aes(age_mid, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
    geom_line(data = go, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = gi, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_point(data = gi, aes(age_mid, mid, colour = model, shape = model, fill = model), size = 2, stroke = 0.4) +
    geom_point(data = go, aes(age_mid, mid, colour = model, shape = model, fill = model), size = 2, stroke = 0.4) +
    scale_models() + guide_models() +
    scale_x_continuous(breaks = c(0, 10, 20, 40, 60, 80)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)),
                       labels = if (m == "prev") scales::percent else waiver()) +
    labs(title = ttl_a[[m]], x = "age (years)", y = lab_a[[m]]) +
    theme_cmp() + theme(legend.position = if (m == "prev") "top" else "none")
}
g <- (panel_age("prev") | panel_age("clin") | panel_age("sev")) + plot_annotation(
  title = sprintf("Age structure of infection and disease at EIR %s", EIR_REF),
  subtitle = "Both models share the immunity functions that shape these profiles; fleet tracks them on a 52-group age grid",
  caption = cap("Points at age-band midpoints; bands are finer in childhood. Severe incidence is the most immunity-sensitive output and the noisiest in the IBM.", ibm_note),
  theme = theme_cmp())
save_fig(g, "core_age", width = 10, height = 5.2)

## ============================================================================
## 2b. core_demography -- custom demography: age structure and age-prevalence
## ============================================================================
if ("demography" %in% age$scenario) {
  d <- age %>% filter(scenario == "demography") %>%
    mutate(dens = pop_frac / (age_hi - age_lo)) %>%
    pivot_longer(c(dens, prev), names_to = "metric", values_to = "y")
  ibm_d <- d %>% filter(model == "IBM") %>% group_by(metric) %>%
    group_modify(~ envelope(.x, by = "age_mid")) %>% ungroup()
  ode_d <- d %>% filter(model == "fleet") %>% rename(mid = y)
  panel_dem <- function(m, ylab, title, pct = FALSE) {
    gi <- filter(ibm_d, metric == m); go <- filter(ode_d, metric == m)
    ggplot() +
      geom_ribbon(data = gi, aes(age_mid, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
      geom_line(data = go, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
      geom_line(data = gi, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
      geom_point(data = gi, aes(age_mid, mid, colour = model, shape = model, fill = model), size = 2, stroke = 0.4) +
      geom_point(data = go, aes(age_mid, mid, colour = model, shape = model, fill = model), size = 2, stroke = 0.4) +
      scale_models() + guide_models() +
      scale_x_continuous(breaks = c(0, 10, 20, 40, 60, 80)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)),
                         labels = if (pct) scales::percent else scales::percent) +
      labs(title = title, x = "age (years)", y = ylab) +
      theme_cmp() + theme(legend.position = if (pct) "none" else "top")
  }
  g <- (panel_dem("dens", "share of the population per year of age",
                  "The age pyramid the mortality schedule implies") |
        panel_dem("prev", "LM prevalence", "Age-prevalence under that demography", TRUE)) +
    plot_annotation(
      title = sprintf("Custom demography at EIR %s: high infant and elderly mortality", EIR_REF),
      subtitle = cap("set_demography() with age-specific death rates from 4.8% per year in infancy to 12% per year over 80; fleet derives its equilibrium age structure from the same schedule", width = 115),
      caption = cap("Population shares are per band divided by band width, so bands of different width are comparable; the bands cover ages 0-85.", ibm_note),
      theme = theme_cmp())
  save_fig(g, "core_demography", width = 10, height = 5)
}

## ============================================================================
## 3. core_seasonal -- the settled annual cycle
## ============================================================================
s <- doy %>% filter(scenario == "seasonal") %>%
  pivot_longer(c(pfpr_2_10, clin_0_5), names_to = "metric", values_to = "y")
ibm_s <- s %>% filter(model == "IBM") %>% group_by(metric) %>%
  group_modify(~ envelope(.x, by = "doy")) %>% ungroup()
ode_s <- s %>% filter(model == "fleet") %>% rename(mid = y)
mon_brk <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))[c(1, 4, 7, 10)] + 1
panel_doy <- function(m, ylab, title, pct = FALSE) {
  gi <- filter(ibm_s, metric == m); go <- filter(ode_s, metric == m)
  ggplot() +
    geom_ribbon(data = gi, aes(doy, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
    geom_line(data = go, aes(doy, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = gi, aes(doy, mid, colour = model, linetype = model), linewidth = 0.8) +
    scale_models(shapes = FALSE) +
    scale_x_continuous(breaks = mon_brk, labels = c("Jan", "Apr", "Jul", "Oct")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)),
                       labels = if (pct) scales::percent else waiver()) +
    labs(title = title, x = NULL, y = ylab) + theme_cmp() +
    theme(legend.position = if (pct) "top" else "none")
}
g <- (panel_doy("pfpr_2_10", PREV_LAB, "Prevalence lags the season", TRUE) |
      panel_doy("clin_0_5", CLIN_LAB, "Incidence follows the rains more sharply")) +
  plot_annotation(
    title = sprintf("Seasonal transmission: the settled annual cycle at EIR %s", EIR_REF),
    subtitle = "Final year of the run, weekly bins. Rainfall enters both models through the larval carrying capacity",
    caption = cap("Both models are seeded at the aseasonal equilibrium and converge onto the same limit cycle during the burn-in.", ibm_note),
    theme = theme_cmp())
save_fig(g, "core_seasonal", width = 10, height = 5)

## ============================================================================
## 4. core_sites -- 63-country monthly comparison (a SNAPSHOT)
## ============================================================================
## The committed cmp_core_sites.png is a snapshot and is deliberately NOT redrawn
## on an ordinary render. Its inputs are a ~7-hour validation run against the
## malariaverse site files, living in a separate checkout that is not part of this
## repo, so on any machine without that checkout this block would otherwise be
## skipped silently -- which is fine for the figure (the committed PNG survives)
## but leaves no record of the fact. Re-draw it deliberately, alongside re-taking
## the statistics in summary_tables.R:
##
##   CMP_REFRESH_SITES=1 Rscript comparison/render_figures.R
##   CMP_REFRESH_SITES=1 Rscript comparison/summary_tables.R
vdir <- VDIR
if (!nzchar(Sys.getenv("CMP_REFRESH_SITES"))) {
  message("core_sites: keeping the committed snapshot ",
          "(CMP_REFRESH_SITES=1 to re-draw it from ", vdir, ")")
  vdir <- ""                                   # skips the block below
}
fs <- list.files(file.path(vdir, "results"), pattern = "_compare.rds$", full.names = TRUE)
if (length(fs)) {
  v <- bind_rows(lapply(fs, readRDS))
  ag <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    list(n = length(x), r = cor(x, y), slope = unname(coef(lm(y ~ x))[2]), bias = mean(y - x) / mean(x)) }
  st_c <- ag(v$ms_clinical, v$mo_clinical); st_s <- ag(v$ms_severe, v$mo_severe)
  hexp <- function(x, y, st, unit, title) {
    d <- data.frame(x = x, y = y) %>% filter(is.finite(x), is.finite(y))
    top <- unname(quantile(c(d$x, d$y), 0.999))
    ggplot(d, aes(x, y)) +
      geom_hex(bins = 60) +
      geom_abline(slope = 1, intercept = 0, colour = REF, linewidth = 0.8,
                  linetype = "22") +
      ## The low end of the ramp is near-white on purpose. The counts are wildly
      ## skewed -- 52% of the hex cells carry 0.2% of the sub-site-months, while
      ## the top 5% of cells carry 90% of them -- so a saturated low end spends
      ## most of the plot's ink on almost none of the data and buries the 1:1
      ## ridge it is meant to show. Keeping log10 keeps the sparse cells visible
      ## (they are real, and the scatter is the point); making them faint stops
      ## them out-shouting the ridge.
      scale_fill_gradient(low = "#F4F6FE", high = "#171449", transform = "log10",
                          name = "sub-site\nmonths",
                          breaks = c(1, 10, 100, 1000, 10000),
                          labels = scales::label_comma()) +
      coord_equal(xlim = c(0, top), ylim = c(0, top), expand = FALSE) +
      labs(title = title,
           subtitle = sprintf("r = %.2f \u00b7 slope = %.2f\nfleet \u2212 IBM on average: %+.1f%% of the IBM mean", st$r, st$slope, 100 * st$bias),
           x = sprintf("IBM (%s)", unit), y = sprintf("fleet (%s)", unit)) +
      theme_cmp() + theme(legend.position = "right", legend.justification = "center",
                          legend.title = element_text(size = rel(0.8), colour = INK2),
                          panel.grid.major = element_blank())
  }
  g <- (hexp(v$ms_clinical, v$mo_clinical, st_c, "episodes per person-year", "Monthly clinical incidence") |
        hexp(v$ms_severe, v$mo_severe, st_s, "episodes per person-year", "Monthly severe incidence")) +
    plot_annotation(
      title = sprintf("Country site files: %s sub-site-months across %d countries",
                      format(st_c$n, big.mark = ","), length(unique(v$iso3c))),
      subtitle = sprintf("Every P. falciparum admin-1 \u00d7 urban/rural sub-site in the malariaverse site files, %d\u2013%d, with its full intervention history",
                         min(v$year), max(v$year)),
      caption = cap("Dashed line = perfect agreement. All ages, P. falciparum only on both sides. Cell colour = number of sub-site-months (log scale). IBM values are the site files' own calibration diagnostic runs; fleet was run here from the same site_parameters() lists."),
      theme = theme_cmp())
  save_fig(g, "core_sites", width = 10, height = 5.4)
}

## ============================================================================
## 5. int_timeseries -- five interventions, prevalence + clinical incidence
## ============================================================================
INT <- names(INT_LABELS)
m <- monthly %>% filter(scenario %in% INT, year >= BURN_Y - 3, year < BURN_Y + 6) %>%
  select(scenario, model, rep, year, pfpr_2_10, clin_0_5) %>%
  pivot_longer(c(pfpr_2_10, clin_0_5), names_to = "metric", values_to = "y") %>%
  mutate(scenario = factor(scenario, levels = INT, labels = INT_LABELS),
         metric = factor(metric, levels = c("pfpr_2_10", "clin_0_5"), labels = c(PREV_LAB, CLIN_LAB)))
ibm_m <- m %>% filter(model == "IBM") %>% group_by(scenario, metric) %>%
  group_modify(~ envelope(.x, by = "year")) %>% ungroup()
ode_m <- m %>% filter(model == "fleet") %>% rename(mid = y)
smc_rounds <- data.frame(scenario = factor(INT_LABELS[["smc"]], levels = INT_LABELS),
                         x = as.vector(sapply(0:2, function(y) BURN_Y + y + (c(0, 30, 60, 90) + 200) / 365)))
onset_lab <- data.frame(scenario = factor(INT_LABELS[[INT[1]]], levels = INT_LABELS),
                        metric = factor(PREV_LAB, levels = c(PREV_LAB, CLIN_LAB)))
g <- ggplot() +
  geom_vline(data = smc_rounds, aes(xintercept = x), colour = GRID, linewidth = 0.5) +
  geom_vline(xintercept = BURN_Y, colour = AXIS, linewidth = 0.5, linetype = "22") +
  geom_text(data = onset_lab, aes(x = BURN_Y, y = Inf, label = "deployment"),
            hjust = -0.08, vjust = 1.6, size = 3.1, colour = INK2, family = FONT) +
  geom_ribbon(data = ibm_m, aes(year, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
  geom_line(data = ode_m, aes(year, mid, colour = model, linetype = model), linewidth = 0.75) +
  geom_line(data = ibm_m, aes(year, mid, colour = model, linetype = model), linewidth = 0.75) +
  facet_grid(scenario ~ metric, scales = "free_y", switch = "y") +
  scale_models(shapes = FALSE) +
  scale_x_continuous(breaks = seq(BURN_Y - 3, BURN_Y + 6, 3),
                     labels = function(b) ifelse(b == BURN_Y, "0", sprintf("%+d y", as.integer(b - BURN_Y)))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Intervention impact: the same deployment through both models",
       subtitle = cap(sprintf("Monthly series, three years before to six after deployment, EIR %s. SMC: EIR 15 in a seasonal setting, three years of rounds (marked)", EIR_REF), width = 120),
       x = "years relative to deployment", y = NULL,
       caption = cap("Each row is one intervention layered on the same baseline with the ordinary malariasimulation set_*() builders.", ibm_note)) +
  theme_cmp() + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 1),
                      strip.placement = "outside", panel.spacing.x = unit(2.2, "lines"))
save_fig(g, "int_timeseries", width = 10, height = 12)

## ============================================================================
## 5b. programme_ts -- long-horizon programmes: nothing / one thing / everything
## ============================================================================
## Section 5 isolates one builder over 6 years, which is the shape for attributing
## a difference. This is the shape for seeing what a programme does: 15 years past
## deployment, seasonal, with the repeated-campaign dynamics visible. Four metrics
## get four SEPARATE faceted columns rather than one facet_grid, because
## facet_grid frees y by row and here the scales differ by COLUMN -- prevalence,
## two incidence families and severe share no axis.
TS <- names(TS_LABELS)
if (all(TS %in% monthly$scenario)) {
  TS_MET <- list(
    list(key = "pfpr_2_10", title = "LM prevalence\nages 2–10", pct = TRUE),
    list(key = "clin_0_5",  title = "Clinical incidence\nages 0–5, per child-year"),
    list(key = "clin_all",  title = "Clinical incidence\nall ages, per person-year"),
    list(key = "sev_all",   title = "Severe incidence\nall ages, per 1,000 py"))
  ts_long <- monthly %>%
    filter(scenario %in% TS, year >= BURN_Y - 2, year <= BURN_Y + TS_YEARS) %>%
    select(scenario, model, rep, year, all_of(vapply(TS_MET, `[[`, "", "key"))) %>%
    pivot_longer(-c(scenario, model, rep, year), names_to = "metric", values_to = "y") %>%
    mutate(scenario = factor(scenario, levels = TS, labels = TS_LABELS))
  ibm_t <- ts_long %>% filter(model == "IBM") %>% group_by(scenario, metric) %>%
    group_modify(~ envelope(.x, by = "year")) %>% ungroup()
  ode_t <- ts_long %>% filter(model == "fleet") %>% rename(mid = y)
  ## net distributions marked only in the two rows that have them -- the SMC
  ## pulses are 4 a year for 15 years and would be a picket fence, so those are
  ## left to show themselves in the clinical sawtooth
  net_x <- data.frame(
    scenario = factor(rep(TS_LABELS[c("ts_nets", "ts_all")], each = 5L),
                      levels = TS_LABELS),
    x = rep(BURN_Y + seq(0, by = TS_NET_EVERY, length.out = 5L), times = 2L))
  ## a heavier envelope than the rest of the figure set (ENV_ALPHA = 0.16). This
  ## figure is 183 monthly points in a ~495 px panel -- 2.7 px per month, 32 px
  ## per seasonal cycle -- so the band only opens at the spike tips, and at 0.16
  ## it is invisible there. It has real width to show: over the months carrying
  ## the top quartile of burden the IBM's 10-90 range is 26% of the panel peak
  ## for all-age severe, against 3-6% for prevalence and the two clinical
  ## measures. Kept local to this figure rather than raised in theme.R, so the
  ## four already-reviewed figures are not changed unseen.
  TS_ENV_ALPHA <- 0.40

  ts_col <- function(m, first) {
    gi <- filter(ibm_t, metric == m$key); go <- filter(ode_t, metric == m$key)
    p <- ggplot() +
      geom_vline(data = net_x, aes(xintercept = x), colour = AXIS, linewidth = 0.45) +
      geom_vline(xintercept = BURN_Y, colour = INK2, linewidth = 0.45, linetype = "22") +
      geom_ribbon(data = gi, aes(year, ymin = lo, ymax = hi, fill = model),
                  alpha = TS_ENV_ALPHA) +
      geom_line(data = go, aes(year, mid, colour = model, linetype = model), linewidth = 0.6) +
      geom_line(data = gi, aes(year, mid, colour = model, linetype = model), linewidth = 0.6) +
      facet_grid(scenario ~ ., switch = "y") +
      scale_models(shapes = FALSE) +
      scale_x_continuous(breaks = BURN_Y + seq(0, TS_YEARS, 5),
                         labels = function(b) paste0("+", as.integer(b - BURN_Y), " y")) +
      labs(title = m$title, x = NULL, y = NULL) +
      theme_cmp() + theme(legend.position = "none")
    p <- p + if (isTRUE(m$pct))
      scale_y_continuous(limits = c(0, NA), labels = scales::percent,
                         expand = expansion(mult = c(0, 0.08)))
      else scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)))
    ## scenario strips on the leftmost column only; repeating them four times
    ## would cost a third of the width and say nothing new
    if (first) p + theme(strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 0.5),
                         strip.placement = "outside")
    else p + theme(strip.text.y = element_blank(), strip.background = element_blank())
  }
  ps <- lapply(seq_along(TS_MET), function(i) ts_col(TS_MET[[i]], i == 1L))
  g <- patchwork::wrap_plots(ps, nrow = 1, widths = c(1.12, 1, 1, 1)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Programmes over fifteen years, through both models",
      subtitle = cap(sprintf("Monthly series at EIR %s in a seasonal setting, from two years before deployment. Every row carries 20%% baseline case management; rows 2–4 add one intervention, row 5 adds all three. Grey rules mark the five net distributions.", EIR_REF), width = 128),
      caption = cap("x = years relative to deployment; the dashed rule is deployment. Nets: 80% coverage every 3 years, 5-year mean retention. SMC: 4 monthly rounds a year, ages 3 months to 5 years, 90% coverage. Case management: SP-AQ, coverage of clinical cases raised from 20% to 60%.",
                    "The IBM band is hard to resolve here -- 15 years of monthly points is ~3 px per month -- so its width is given instead: over the months carrying the top quartile of burden the 10-90% replicate range is 26% of the panel peak for severe incidence, and 3-6% for the other three. Severe is the noisiest because a 30-day bin holds only ~20 severe episodes at the seasonal peak in a population of 10,000.",
                    ibm_note),
      theme = theme_cmp()) &
    theme(legend.position = "top", legend.justification = "left")
  save_fig(g, "programme_ts", width = 13.5, height = 11)
} else {
  message("programme_ts skipped: ts_* scenarios not in rep_monthly.csv")
}

## ============================================================================
## 6. int_impact -- % reduction over the first three years, IBM range vs fleet
## ============================================================================
## Four outcomes: the two young-child measures a trial would report, and the two
## all-age measures a programme carries. Severe is the noisiest of them in the
## IBM, which the replicate range shows honestly.
MET_KEY <- c(pfpr = "LM prevalence, ages 2\u201310",
             clin05 = "clinical incidence, ages 0\u20135",
             clinall = "clinical incidence, all ages",
             sevall = "severe incidence, all ages")
MET_R <- unname(MET_KEY)
red <- monthly %>% filter(scenario %in% INT) %>%
  mutate(phase = case_when(year >= BURN_Y - 3 & year < BURN_Y ~ "pre",
                           year >= BURN_Y & year < BURN_Y + 3 ~ "post", TRUE ~ NA_character_)) %>%
  filter(!is.na(phase)) %>%
  group_by(scenario, model, rep, phase) %>%
  summarise(pfpr = mean(pfpr_2_10), clin05 = mean(clin_0_5),
            clinall = mean(clin_all), sevall = mean(sev_all), .groups = "drop") %>%
  pivot_longer(all_of(names(MET_KEY)), names_to = "metric", values_to = "v") %>%
  pivot_wider(names_from = phase, values_from = v) %>%
  mutate(reduction = 1 - post / pre,
         metric = unname(MET_KEY[metric])) %>%
  select(scenario, model, rep, metric, reduction)
ibm_r <- red %>% filter(model == "IBM") %>% group_by(scenario, metric) %>%
  summarise(mid = median(reduction), lo = unname(quantile(reduction, .1)),
            hi = unname(quantile(reduction, .9)), .groups = "drop") %>% mutate(model = "IBM")
ode_r <- red %>% filter(model == "fleet") %>% transmute(scenario, metric, mid = reduction, model = "fleet")
both <- bind_rows(ibm_r, ode_r) %>%
  mutate(scenario = factor(scenario, levels = INT, labels = INT_LABELS),
         metric = factor(metric, levels = MET_R))
## Round before writing. This file is committed and CI byte-compares it against a
## fresh render, but the values come out of quantile() on different hardware, and
## cross-platform floating point is not bit-identical: a Linux runner reproduced
## 0.324024697880216 where this machine wrote ...217, one ulp, and the check
## failed on it. These are summary percentages for an article, so nothing needs
## more than a few significant figures; 10 is far beyond the reporting precision
## and far above the platform noise, which makes the artifact reproducible.
write.csv(dplyr::mutate(both, dplyr::across(where(is.numeric), ~ signif(.x, 10))),
          file.path(DDIR, "int_impact_summary.csv"), row.names = FALSE)
seg  <- both %>% select(scenario, metric, model, mid) %>% pivot_wider(names_from = model, values_from = mid)
XCOL <- c(IBM = 1.08, fleet = 1.22)                 # value columns to the right of the data
vals <- both %>% mutate(x = XCOL[model], lab = scales::percent(mid, accuracy = 1))
hdr  <- data.frame(x = XCOL, lab = names(XCOL))
xmin <- min(-0.04, floor(min(c(both$lo, both$mid), na.rm = TRUE) * 20) / 20 - 0.03)

g <- ggplot(both, aes(y = scenario)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.5) +
  geom_segment(data = seg, aes(x = IBM, xend = fleet, yend = scenario), colour = GRID,
               linewidth = 2.2, lineend = "round") +
  geom_linerange(data = filter(both, model == "IBM"), aes(xmin = lo, xmax = hi, colour = model),
                 linewidth = 0.9, alpha = 0.55) +
  geom_point(aes(x = mid, colour = model, shape = model, fill = model), size = 3.2, stroke = 0.6) +
  geom_text(data = vals, aes(x = x, label = lab), hjust = 0, size = 3.3, colour = INK2, family = FONT) +
  geom_text(data = hdr, aes(x = x, y = Inf, label = lab), hjust = 0, vjust = 1.4, size = 3.3,
            fontface = "bold", colour = INK2, family = FONT) +
  facet_wrap(~metric, nrow = 2) +
  scale_models(lines = FALSE) + guide_models() +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25), limits = c(xmin, 1.36),
                     expand = expansion(0)) +
  scale_y_discrete(limits = rev(unname(INT_LABELS)), expand = expansion(add = c(0.6, 1.3))) +
  coord_cartesian(clip = "off") +
  labs(title = "Intervention impact summarised: reduction over the first three years",
       subtitle = "Relative to the three pre-deployment years of the same run. Circle = IBM median with 10\u201390% replicate range; triangle = fleet",
       x = "reduction relative to baseline", y = NULL,
       caption = cap("Columns give the plotted medians.", ibm_note)) +
  theme_cmp() + theme(panel.grid.major.y = element_blank(), axis.line.x = element_blank(),
                      axis.text.y = element_text(size = rel(0.9), lineheight = 0.95, hjust = 1),
                      panel.spacing.x = unit(1.6, "lines"),
                      panel.spacing.y = unit(2.0, "lines"))
save_fig(g, "int_impact", width = 11, height = 8.6)

## ---- console summary ---------------------------------------------------------
cat("figures written", if (SMOKE) "to comparison/plots/smoke" else "to man/figures and vignettes", "\n\n")
print(both %>% select(scenario, metric, model, mid) %>% mutate(mid = round(mid, 3)) %>%
        pivot_wider(names_from = model, values_from = mid) %>% mutate(scenario = sub("\n.*", "", scenario)), n = 20)
cat("\nrealised EIR (final 3 years):\n")
print(eq %>% filter(grepl("^eir_", scenario)) %>% group_by(scenario, model) %>%
        summarise(eir = round(mean(eir_realised), 2), pfpr = round(mean(pfpr_2_10), 3), .groups = "drop") %>%
        pivot_wider(names_from = model, values_from = c(eir, pfpr)), n = 20)
cat("\nrun time:\n")
print(timing %>% group_by(model) %>%
        summarise(runs = n(), s_per_sim_year = round(sum(elapsed_s) / sum(years), 2),
                  mean_run_s = round(mean(elapsed_s), 1), .groups = "drop"))
