# Build every comparison figure from the saved CSVs (no model runs).
#
#   Rscript comparison/render_figures.R                 # figures -> man/figures + vignettes
#   CMP_SMOKE=1 Rscript comparison/render_figures.R     # from data/smoke -> comparison/plots/smoke
#
# Writes cmp_*.png. All visual decisions live in theme.R; this file only shapes
# data and composes panels. Figures:
#   core_eir       PfPR(2-10) and under-5 clinical incidence vs EIR
#   core_age       age profiles at EIR 20: prevalence, clinical, severe
#   core_seasonal  the settled annual cycle: prevalence and clinical incidence
#   core_sites     63-country monthly comparison (from blink2_validate results)
#   int_timeseries five interventions x {prevalence, clinical}, time series
#   int_impact     % reduction per intervention, IBM (replicate range) vs blink

.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")
suppressMessages({library(dplyr); library(tidyr)})
ROOT <- "C:/Users/pwinskil/Documents/dev/blink2/blink"
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
  "blink: one deterministic run seeded at equilibrium."), n_rep, format(POP, big.mark = ","), BURN_Y)
PREV_LAB <- "LM prevalence, ages 2\u201310"
CLIN_LAB <- "clinical episodes per child-year, ages 0\u20135"

## ============================================================================
## 1. core_eir -- equilibrium PfPR(2-10) and under-5 clinical incidence vs EIR
## ============================================================================
e_long <- eq %>% filter(grepl("^eir_", scenario)) %>%
  mutate(init_EIR = as.numeric(sub("eir_", "", scenario))) %>%
  select(init_EIR, model, rep, pfpr_2_10, clin_0_5) %>%
  pivot_longer(c(pfpr_2_10, clin_0_5), names_to = "metric", values_to = "y")
ibm_e <- e_long %>% filter(model == "IBM") %>% group_by(metric) %>%
  group_modify(~ envelope(.x, by = "init_EIR")) %>% ungroup()
ode_e <- e_long %>% filter(model == "blink") %>% rename(mid = y)

panel_eir <- function(metric_id, ylab, title) {
  gi <- filter(ibm_e, metric == metric_id); go <- filter(ode_e, metric == metric_id)
  ggplot() +
    geom_linerange(data = gi, aes(init_EIR, ymin = lo, ymax = hi, colour = model),
                   linewidth = 0.7, alpha = 0.55, show.legend = FALSE) +
    geom_line(data = gi, aes(init_EIR, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = go, aes(init_EIR, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_point(data = gi, aes(init_EIR, mid, colour = model, shape = model, fill = model),
               size = 2.6, stroke = 0.5) +
    geom_point(data = go, aes(init_EIR, mid, colour = model, shape = model, fill = model),
               size = 2.6, stroke = 0.5) +
    scale_x_log10(breaks = EIR_GRID, minor_breaks = NULL) +
    scale_models() + guide_models() +
    labs(title = title, x = "EIR (infectious bites per adult per year)", y = ylab) +
    theme_cmp() + theme(legend.position = "none")
}
p1 <- panel_eir("pfpr_2_10", PREV_LAB, "Parasite prevalence rises with transmission") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme(legend.position = "top")
p2 <- panel_eir("clin_0_5", CLIN_LAB, "Clinical incidence saturates as immunity builds") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)))
g <- (p1 | p2) + plot_annotation(
  title = "Core transmission relationships at equilibrium",
  subtitle = "The same parameter list through both models, across a 120-fold range of transmission intensity",
  caption = cap("x = the EIR passed to set_equilibrium(); each model's realised EIR is reported in the article.", ibm_note),
  theme = theme_cmp())
save_fig(g, "core_eir", width = 10, height = 5)

## ============================================================================
## 2. core_age -- age profiles at the reference EIR
## ============================================================================
a <- age %>% filter(scenario == paste0("eir_", EIR_REF)) %>%
  pivot_longer(c(prev, clin, sev), names_to = "metric", values_to = "y")
ibm_a <- a %>% filter(model == "IBM") %>% group_by(metric) %>%
  group_modify(~ envelope(.x, by = "age_mid")) %>% ungroup()
ode_a <- a %>% filter(model == "blink") %>% rename(mid = y)
lab_a <- c(prev = "LM prevalence", clin = "clinical episodes per person-year",
           sev = "severe episodes per 1,000 person-years")
ttl_a <- c(prev = "Prevalence peaks in\nschool-age children",
           clin = "Clinical disease\nconcentrates in the young",
           sev  = "Severe disease is rarer\nand earlier still")
panel_age <- function(m) {
  gi <- filter(ibm_a, metric == m); go <- filter(ode_a, metric == m)
  ggplot() +
    geom_ribbon(data = gi, aes(age_mid, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
    geom_line(data = gi, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = go, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
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
  subtitle = "Both models share the immunity functions that shape these profiles; blink tracks them on a 52-group age grid",
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
  ode_d <- d %>% filter(model == "blink") %>% rename(mid = y)
  panel_dem <- function(m, ylab, title, pct = FALSE) {
    gi <- filter(ibm_d, metric == m); go <- filter(ode_d, metric == m)
    ggplot() +
      geom_ribbon(data = gi, aes(age_mid, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
      geom_line(data = gi, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
      geom_line(data = go, aes(age_mid, mid, colour = model, linetype = model), linewidth = 0.8) +
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
      subtitle = cap("set_demography() with age-specific death rates from 4.8% per year in infancy to 12% per year over 80; blink derives its equilibrium age structure from the same schedule", width = 115),
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
ode_s <- s %>% filter(model == "blink") %>% rename(mid = y)
mon_brk <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))[c(1, 4, 7, 10)] + 1
panel_doy <- function(m, ylab, title, pct = FALSE) {
  gi <- filter(ibm_s, metric == m); go <- filter(ode_s, metric == m)
  ggplot() +
    geom_ribbon(data = gi, aes(doy, ymin = lo, ymax = hi, fill = model), alpha = ENV_ALPHA) +
    geom_line(data = gi, aes(doy, mid, colour = model, linetype = model), linewidth = 0.8) +
    geom_line(data = go, aes(doy, mid, colour = model, linetype = model), linewidth = 0.8) +
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
## 4. core_sites -- 63-country monthly comparison (site::site_parameters lists)
## ============================================================================
vdir <- "C:/Users/pwinskil/Documents/dev/blink2/blink2_validate"
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
           subtitle = sprintf("r = %.2f \u00b7 slope = %.2f\nblink \u2212 IBM on average: %+.1f%% of the IBM mean", st$r, st$slope, 100 * st$bias),
           x = sprintf("IBM (%s)", unit), y = sprintf("blink (%s)", unit)) +
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
      caption = cap("Dashed line = perfect agreement. All ages, P. falciparum only on both sides. Cell colour = number of sub-site-months (log scale). IBM values are the site files' own calibration diagnostic runs; blink was run here from the same site_parameters() lists."),
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
ode_m <- m %>% filter(model == "blink") %>% rename(mid = y)
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
  geom_line(data = ibm_m, aes(year, mid, colour = model, linetype = model), linewidth = 0.75) +
  geom_line(data = ode_m, aes(year, mid, colour = model, linetype = model), linewidth = 0.75) +
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
## 6. int_impact -- % reduction over the first three years, IBM range vs blink
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
ode_r <- red %>% filter(model == "blink") %>% transmute(scenario, metric, mid = reduction, model = "blink")
both <- bind_rows(ibm_r, ode_r) %>%
  mutate(scenario = factor(scenario, levels = INT, labels = INT_LABELS),
         metric = factor(metric, levels = MET_R))
write.csv(both, file.path(DDIR, "int_impact_summary.csv"), row.names = FALSE)
seg  <- both %>% select(scenario, metric, model, mid) %>% pivot_wider(names_from = model, values_from = mid)
XCOL <- c(IBM = 1.08, blink = 1.22)                 # value columns to the right of the data
vals <- both %>% mutate(x = XCOL[model], lab = scales::percent(mid, accuracy = 1))
hdr  <- data.frame(x = XCOL, lab = names(XCOL))
xmin <- min(-0.04, floor(min(c(both$lo, both$mid), na.rm = TRUE) * 20) / 20 - 0.03)

g <- ggplot(both, aes(y = scenario)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.5) +
  geom_segment(data = seg, aes(x = IBM, xend = blink, yend = scenario), colour = GRID,
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
       subtitle = "Relative to the three pre-deployment years of the same run. Circle = IBM median with 10\u201390% replicate range; triangle = blink",
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
