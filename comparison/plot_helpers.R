# Shared plot builders for the IBM-vs-ODE comparison. Sourced by both
# run_comparison.R (live) and render_plots.R (re-render from saved CSVs), so the
# two never drift. Each builder takes a tidy/ wide data frame and returns a ggplot.
#
# Design (colourblind- and greyscale-safe): the IBM/ODE distinction is carried by
# THREE redundant channels — colour, line type, and point shape — so the series
# stay distinguishable in monochrome and under any colour-vision deficiency.

suppressMessages({library(ggplot2); library(tidyr)})

BURN_Y <- 30L          # IBM burn-in years (must match run_comparison.R)
POP    <- 10000L       # IBM population
AGE_EDGES <- c(0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 85)
band_lo <- utils::head(AGE_EDGES, -1) * 365
band_hi <- utils::tail(AGE_EDGES, -1) * 365
# SMC round timesteps (deterministic; mirrors scenario E in run_comparison.R)
SMC_ROUNDS <- as.vector(sapply(0:2, function(y) BURN_Y * 365 + y * 365 + c(0, 30, 60, 90) + 200))

# red vs dark-navy: distinct hue AND luminance (L ~0.14 vs ~0.06), so the pair
# survives greyscale/monochrome, not just colour-vision deficiency. Shape + line
# type add two further redundant channels.
cols   <- c(IBM = "#C0392B", ODE = "#12436D")   # red / dark navy
shapes <- c(IBM = 16, ODE = 17)                 # filled circle / triangle
ltys   <- c(IBM = 2, ODE = 1)                   # dashed / solid
LM_CAP <- "LM = light-microscopy-detectable prevalence (proportion 0-1)"

theme_cmp <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 plot.caption = ggplot2::element_text(hjust = 0, colour = "grey40"))

# colour + linetype scales (both series in every plot); point plots add the shape
# scale themselves so line-only plots don't warn about an unused shape scale.
.model_scales <- function()
  list(scale_colour_manual(values = cols), scale_linetype_manual(values = ltys),
       labs(colour = NULL, shape = NULL, linetype = NULL))
.shape_scale <- function() scale_shape_manual(values = shapes)

.longify <- function(d) pivot_longer(d, c(IBM, ODE), names_to = "model", values_to = "y")

# Incidence time series: clinical (per person-year) + severe (per 1000 person-yr)
# faceted vertically, IBM vs ODE. `d` is long: year, inc, model, metric. A single
# `vlines` value is drawn dark + labelled (a key event); many values are drawn
# light + unlabelled (e.g. repeated SMC rounds). The IBM line is drawn at reduced
# opacity so its Monte-Carlo severe-panel noise doesn't obscure the ODE mean.
plot_incidence <- function(d, title, subtitle, vlines = NULL, vlabel = NULL, xmin = NULL) {
  if (!is.null(xmin)) d <- subset(d, year >= xmin)
  d$metric <- factor(d$metric, levels = INC_LEVELS)   # defensive re-factor (read.csv drops it)
  g <- ggplot(d, aes(year, inc, colour = model))
  if (!is.null(vlines)) g <- g + geom_vline(xintercept = vlines, linetype = 3,
    colour = if (length(vlines) == 1) "grey40" else "grey75")
  if (!is.null(vlabel)) g <- g + annotate("text", x = vlines[1], y = Inf,
    label = vlabel, hjust = -0.05, vjust = 1.6, size = 3, colour = "grey30")
  g + geom_line(aes(linetype = model, alpha = model), linewidth = .9) +
    scale_alpha_manual(values = c(IBM = .6, ODE = 1), guide = "none") + .model_scales() +
    facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left") +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = title, subtitle = subtitle, x = "year", y = NULL,
         caption = paste("Severe events are rare; the IBM severe series reflects Monte-Carlo",
                         "sampling noise in the under-5 subpopulation.")) +
    theme_cmp + theme(strip.placement = "outside", strip.background = element_blank(),
                      strip.text.y.left = element_text(angle = 90))
}
# metric factor with unit-bearing labels, clinical first
INC_LEVELS <- c("Clinical incidence\n(per person-year)",
                "Severe incidence\n(per 1000 person-years)")
# scenario metadata shared by run_incidence.R (live) and render_plots.R (re-plot)
INC_SPECS <- list(
  nets = list(file = "inc_nets",
    title = "Bed-net campaign: under-5 incidence (80% coverage at year 30)",
    subtitle = "monthly, annualised; ages 0-5", vlines = BURN_Y, vlabel = "nets deployed",
    xmin = BURN_Y - 2),
  seasonal = list(file = "inc_seasonal",
    title = "Seasonal transmission: under-5 incidence cycle (EIR~20)",
    subtitle = "monthly, annualised; ages 0-5 (final 3 years)", vlines = NULL, vlabel = NULL,
    xmin = BURN_Y),
  smc = list(file = "inc_smc",
    title = "Seasonal SMC: under-5 incidence (4 rounds/yr x 3 yr from year 30)",
    subtitle = "monthly, annualised; ages 0-5 (SMC targets 0.25-5y)",
    vlines = SMC_ROUNDS / 365, vlabel = NULL, xmin = BURN_Y - 2))
render_incidence <- function(spec, d)
  plot_incidence(d, spec$title, spec$subtitle, spec$vlines, spec$vlabel, spec$xmin)

# A — equilibrium PfPR(2-10) vs EIR (log-x). Points at nominal input EIR.
plot_A <- function(d) {
  dl <- pivot_longer(d[, c("init_EIR", "IBM", "ODE")], c(IBM, ODE),
                     names_to = "model", values_to = "pfpr")
  ggplot(dl, aes(init_EIR, pfpr, colour = model)) +
    geom_line(aes(linetype = model), linewidth = 1) +
    geom_point(aes(shape = model), size = 2.6) +
    scale_x_log10() + .model_scales() + .shape_scale() + coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Equilibrium PfPR(2-10) vs EIR",
         subtitle = sprintf("IBM burned in %d y, pop %d; ODE at fixed point. x = nominal (input) EIR.",
                            BURN_Y, POP),
         x = "EIR (bites/adult/year, log scale)", y = "LM prevalence (ages 2-10)",
         caption = LM_CAP) + theme_cmp
}

# B — age-prevalence profile at EIR 20
plot_B <- function(d) {
  dl <- .longify(d); names(dl)[names(dl) == "y"] <- "pfpr"
  ggplot(dl, aes(age_mid, pfpr, colour = model)) +
    geom_line(aes(linetype = model), linewidth = 1) + geom_point(aes(shape = model), size = 2) +
    .model_scales() + .shape_scale() + coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Age-prevalence profile (EIR 20)",
         subtitle = sprintf("IBM burned in %d y, pop %d", BURN_Y, POP),
         x = "age (years)", y = "LM prevalence", caption = LM_CAP) + theme_cmp
}

# C — bed-net campaign deployed at year BURN_Y
plot_C <- function(d) {
  ggplot(subset(d, year >= BURN_Y - 2), aes(year, pfpr, colour = model)) +
    geom_vline(xintercept = BURN_Y, linetype = 3, colour = "grey40") +
    annotate("text", x = BURN_Y, y = Inf, label = "nets deployed", hjust = -0.05,
             vjust = 1.6, size = 3, colour = "grey30") +
    geom_line(aes(linetype = model), linewidth = .9) + .model_scales() +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = "Bed-net campaign (80% coverage at year 30)",
         subtitle = "monthly LM prevalence, ages 2-10", x = "year",
         y = "LM prevalence (ages 2-10)", caption = LM_CAP) + theme_cmp
}

# D — seasonal transmission, settled annual cycle
plot_D <- function(d) {
  ggplot(d, aes(doy, pfpr, colour = model)) +
    geom_line(aes(linetype = model), linewidth = .9) + .model_scales() +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = "Seasonal transmission: settled annual cycle (EIR~20)",
         subtitle = "final year, LM prevalence 2-10 (IBM 15-day smoothed)",
         x = "day of year", y = "LM prevalence (ages 2-10)", caption = LM_CAP) + theme_cmp
}

# E — seasonal SMC in under-5s (round markers shown)
plot_E <- function(d) {
  ggplot(subset(d, year >= BURN_Y - 2), aes(year, pfpr, colour = model)) +
    geom_vline(xintercept = SMC_ROUNDS / 365, linetype = 3, colour = "grey75") +
    geom_line(aes(linetype = model), linewidth = .9) + .model_scales() +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = "Seasonal SMC in under-5s (4 rounds/yr x 3 yr from year 30)",
         subtitle = "monthly LM prevalence; SMC targets 0.25-5y, nearest rendered band 0-5 shown",
         x = "year", y = "LM prevalence (ages 0-5)", caption = LM_CAP) + theme_cmp
}

# F1 — custom-demography equilibrium age structure, as DENSITY per year of age
# (fraction / band width) so unequal band widths don't distort the shape.
plot_F_structure <- function(ds) {
  wyr <- (band_hi - band_lo) / 365
  ds$IBM <- ds$IBM / wyr; ds$ODE <- ds$ODE / wyr
  dl <- .longify(ds); names(dl)[names(dl) == "y"] <- "dens"
  ggplot(dl, aes(age_mid, dens, colour = model)) +
    geom_line(aes(linetype = model), linewidth = 1) + geom_point(aes(shape = model), size = 2) +
    .model_scales() + .shape_scale() + coord_cartesian(ylim = c(0, NA)) +
    labs(title = "Custom demography: equilibrium age structure",
         subtitle = "high infant & elderly mortality (density normalised by band width)",
         x = "age (years)", y = "population density (fraction per year of age)") + theme_cmp
}

# F2 — custom-demography age-prevalence
plot_F_prev <- function(dp) {
  dl <- .longify(dp); names(dl)[names(dl) == "y"] <- "pfpr"
  ggplot(dl, aes(age_mid, pfpr, colour = model)) +
    geom_line(aes(linetype = model), linewidth = 1) + geom_point(aes(shape = model), size = 2) +
    .model_scales() + .shape_scale() + coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Custom demography: age-prevalence",
         subtitle = sprintf("IBM burned in %d y, pop %d", BURN_Y, POP),
         x = "age (years)", y = "LM prevalence", caption = LM_CAP) + theme_cmp
}
