# House style for the blink-vs-malariasimulation comparison figures.
#
# Sourced by run_replicates.R (live) and render_figures.R (re-render from saved
# CSVs), so the two never drift. Everything visual lives here: palette, theme,
# series scales, and the small helpers every figure shares.
#
# Design rules (see the dataviz method the figures were built against):
#   * Two series, IBM vs blink, encoded THREE ways -- colour, line type, point
#     shape -- so every panel reads in greyscale and under colour-vision
#     deficiency. Palette validated: indigo/coral pass CVD dE 27 (protan) and
#     normal-vision dE 38; both >= 3:1 on white.
#   * The IBM is stochastic. It is drawn as the median of N replicates with a
#     10-90% envelope at ~15% opacity, never as one noisy realisation.
#   * Thin marks, hairline SOLID gridlines one step off the surface, no panel
#     border, legend at the top, titles left-aligned. Text never wears a series
#     colour.
#   * One y-axis per panel. Two measures = two panels (patchwork), never a dual axis.

suppressMessages({library(ggplot2); library(patchwork)})

## ---- palette (site colours: _pkgdown.yml primary / danger) ------------------
COL <- c(IBM = "#E5533F", blink = "#4338CA")
LTY <- c(IBM = "22",      blink = "solid")        # dashed / solid
SHP <- c(IBM = 21,        blink = 24)             # filled circle / triangle
INK   <- "#111827"; INK2 <- "#52514E"; MUTED <- "#898781"
## reference lines (1:1 agreement). Amber, because it is the complement of the
## indigo hex ramp and is NOT a series colour, so it can never be read as a
## model; drawn over a white halo so it stays visible on both the white surface
## and the darkest hex fill (5.4:1 on white).
REF   <- "#B45309"
GRID  <- "#E5E7EB"; AXIS <- "#C9CCD1"; SURFACE <- "#FFFFFF"
ENV_ALPHA <- 0.16                                 # IBM envelope wash

FONT <- if ("Segoe UI" %in% systemfonts::system_fonts()$family) "Segoe UI" else "sans"

theme_cmp <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = FONT) +
    theme(
      text = element_text(colour = INK),
      plot.title = element_text(face = "bold", size = rel(1.15), hjust = 0,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(colour = INK2, size = rel(0.95), hjust = 0,
                                   margin = margin(b = 10)),
      plot.caption = element_text(colour = MUTED, size = rel(0.78), hjust = 0,
                                  margin = margin(t = 10)),
      plot.title.position = "plot", plot.caption.position = "plot",
      panel.grid.major = element_line(colour = GRID, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = AXIS, linewidth = 0.4),
      axis.ticks = element_blank(),
      axis.title = element_text(colour = INK2, size = rel(0.9)),
      axis.text = element_text(colour = INK2, size = rel(0.85)),
      legend.position = "top", legend.justification = "left",
      legend.title = element_blank(), legend.key.width = unit(1.6, "lines"),
      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(6, "pt"),
      strip.text = element_text(face = "bold", colour = INK, hjust = 0,
                                size = rel(0.95), margin = margin(b = 6)),
      strip.background = element_blank(), strip.placement = "outside",
      plot.background = element_rect(fill = SURFACE, colour = NA),
      panel.spacing = unit(1.4, "lines"),
      plot.margin = margin(10, 14, 8, 10)
    )
}

## series scales -- every panel that draws both models uses exactly these
scale_models <- function(shapes = TRUE, lines = TRUE) {
  s <- list(scale_colour_manual(values = COL, breaks = c("IBM", "blink")),
            scale_fill_manual(values = COL, breaks = c("IBM", "blink"), guide = "none"),
            labs(colour = NULL, linetype = NULL, shape = NULL))
  if (lines)  s <- c(s, list(scale_linetype_manual(values = LTY, breaks = c("IBM", "blink"))))
  if (shapes) s <- c(s, list(scale_shape_manual(values = SHP, breaks = c("IBM", "blink"))))
  s
}
## a legend that shows line + point together, in the same order everywhere
guide_models <- function() guides(
  colour = guide_legend(override.aes = list(linewidth = 0.9, size = 2.6)))
## captions do not wrap on their own; fold them at the width a 10-inch figure holds
cap <- function(..., width = 135) paste(strwrap(paste(...), width = width), collapse = "\n")

## ---- data helpers -----------------------------------------------------------
## IBM replicate summary: median and 10-90% band per x, for a value column.
envelope <- function(d, by, value = "y") {
  stopifnot(all(c(by, "rep", value) %in% names(d)))
  d <- d[is.finite(d[[value]]), ]
  agg <- function(f) aggregate(d[[value]], d[by], f)
  out <- agg(stats::median); names(out)[ncol(out)] <- "mid"
  out$lo <- agg(function(v) unname(stats::quantile(v, 0.10)))[[length(by) + 1]]
  out$hi <- agg(function(v) unname(stats::quantile(v, 0.90)))[[length(by) + 1]]
  out$model <- "IBM"
  out
}

## the two marks every time-series panel draws: IBM envelope + median, blink line
geom_ibm_envelope <- function(mapping_x, lwd = 0.8) list(
  geom_ribbon(aes(x = {{ mapping_x }}, ymin = lo, ymax = hi, fill = model),
              alpha = ENV_ALPHA, colour = NA),
  geom_line(aes(x = {{ mapping_x }}, y = mid, colour = model, linetype = model),
            linewidth = lwd, lineend = "round"))

## save to BOTH homes: man/figures (README, GitHub) and vignettes (pkgdown article)
save_fig <- function(g, name, width, height, dpi = 200) {
  root <- "C:/Users/pwinskil/Documents/dev/blink2/blink"
  for (dir in c("man/figures", "vignettes")) {
    f <- file.path(root, dir, paste0("cmp_", name, ".png"))
    ggsave(f, g, width = width, height = height, dpi = dpi, device = ragg::agg_png,
           bg = SURFACE)
  }
  invisible(g)
}

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
