# Generate the blink model flow diagram used at the top of vignette("model").
#
#   Rscript pkgdown/diagrams/make_model_flow.R
#
# Needs flodia (GitHub only, deliberately NOT a package dependency):
#   remotes::install_github("mrc-ide/flodia")
#
# Outputs vignettes/model_flow.png and man/figures/model_flow.png, both
# committed, so neither the vignette nor the README needs flodia at build time.
#
# ---------------------------------------------------------------------------
# DESIGN NOTES
#
# * SCOPE. This is an overview. Antimalarial resistance is deliberately out of
#   scope, so the treated state is drawn as a single T. The model itself splits
#   it into T and T_slow when resistance is switched on -- see the vignette
#   appendix, section B.4.
# * SIZE. The figure is displayed at ~520-700 px, not the 2400 px it is authored
#   at. Sized for the small end: no legend blocks (they live in the caption), no
#   compound algebra on the arrows, canvas kept narrow so the same cex renders
#   larger. ALWAYS check a 700 px render before committing.
# * SOURCES. Flows leave the set of compartments they actually leave. Where a
#   quantity is a sum over several compartments -- the at-risk pool, the
#   infectious pool, the adult mosquito pool -- the arrow departs from a
#   bracket, never from one representative box.
# * GRAMMAR. Solid black = a flow of individuals. Dashed colour = not a flow:
#   the two transmission couplings set a rate in the other panel, and the
#   chemoprevention pulse is an instantaneous state jump.
# * GEOMETRY. Every path is orthogonal, and every arrow terminates on something.
#   No line may pass through a box; no two flows may share a lane. S takes its
#   three returns on three different sides: P from above, P_c from below, U from
#   the left.
# * CLASS is encoded by border weight as well as fill, so it survives greyscale
#   and colour-vision deficiency.
# ---------------------------------------------------------------------------

library(flodia)

## ---- palette ---------------------------------------------------------------
col_sus   <- light_palette("bu")
col_clin  <- light_palette("rd")
col_asym  <- light_palette("oryl")
col_treat <- light_palette("gn")
col_proph <- light_palette("pu")
col_aq    <- light_palette("gybr")
col_mosq  <- light_palette("gnbu")

COUPLE <- "#464273"   # transmission coupling
PULSE  <- "#B03A22"   # discrete chemoprevention pulse
GREY   <- "#666666"   # all annotation text (5.7:1 on white)
POOL   <- "#7A6E7E"   # bracket / grouping rules
INK    <- "grey15"

RX <- 0.30; RY <- 0.20

nd <- function(x, y, label, col, heavy = FALSE) {
  node(x = x, y = y, rx = RX, ry = RY, label = label, node_col = col,
       border_col = INK, label_cex = 1.45, lwd = if (heavy) 2.6 else 1.0)
}
pt <- function(x, y) node(x = x, y = y, r = 0)

## orthogonal polyline; arrowhead on the final segment unless head = FALSE
path <- function(xs, ys, col = INK, lty = 1, lwd = 1.2, head = TRUE) {
  n <- length(xs)
  if (n > 2) segments(xs[1:(n - 2)], ys[1:(n - 2)], xs[2:(n - 1)],
                      ys[2:(n - 1)], col = col, lty = lty, lwd = lwd)
  if (head) arrows(xs[n - 1], ys[n - 1], xs[n], ys[n], length = 0.11,
                   angle = 20, col = col, lty = lty, lwd = lwd)
  else segments(xs[n - 1], ys[n - 1], xs[n], ys[n], col = col, lty = lty,
                lwd = lwd)
}

bracket <- function(x0, x1, y0, y1) {
  rect(x0, y0, x1, y1, border = POOL, lty = 3, lwd = 1.8,
       col = adjustcolor(POOL, alpha.f = 0.05))
}

death <- function(x, y0, len = 0.30) {
  arrows(x, y0, x, y0 - len, length = 0.09, angle = 20, col = INK, lwd = 1.1)
}

lab <- function(x, y, txt, col = GREY, cex = 0.95, adj = c(0.5, 0.5), font = 1) {
  text(x, y, txt, col = col, cex = cex, adj = adj, font = font)
}

model_flow <- function() {

  SPINE    <- 2.10   # x of the vertical infection spine
  COL      <- 3.60   # x of the infection-outcome column
  FOOT     <- 1.85   # y of the at-risk feeder lane and the spine foot
  LANE_U   <- 1.26   # y of the U -> S return lane
  LANE_CP  <- 1.00   # y of the chemoprevention lane
  LANE_EIR <- 0.22   # y of the EIR coupling lane
  MY       <- -0.35  # y of the mosquito row

  ## ======================= HUMAN =========================================
  bracket(3.24, 5.62, 2.04, 4.32)                    # infectious pool

  Ph <- nd(0.80, 4.15, expression(bold(P)),    col_proph)
  S  <- nd(0.80, 3.15, expression(bold(S)),    col_sus)
  Pc <- nd(0.80, 2.30, expression(bold(P[c])), col_proph)
  Tr <- nd(COL,  4.00, expression(bold(T)),    col_treat, heavy = TRUE)
  D  <- nd(COL,  3.15, expression(bold(D)),    col_clin,  heavy = TRUE)
  A  <- nd(COL,  2.30, expression(bold(A)),    col_asym,  heavy = TRUE)
  U  <- nd(5.25, 2.30, expression(bold(U)),    col_asym,  heavy = TRUE)

  lab(4.66, 3.66, "infectious to\nmosquitoes", col = POOL, cex = 0.92, font = 3)

  ## --- infection, out of the at-risk pool S + A + U ---
  flowx(from = S, to = pt(SPINE, S$y), label = expression(Lambda),
        label_pos = 0.55, label_cex = 1.2)
  segments(SPINE, FOOT, SPINE, Tr$y, col = INK, lwd = 1.2)     # the spine
  path(c(A$x, A$x, SPINE), c(A$y0, FOOT, FOOT), head = FALSE)  # A is at risk
  path(c(U$x, U$x, SPINE), c(U$y0, FOOT, FOOT), head = FALSE)  # so is U
  lab(SPINE + 0.20, FOOT - 0.26,
      expression(paste(Lambda, " acts on ", S + A + U)), adj = c(0, 0.5),
      cex = 0.92)

  flowx(from = pt(SPINE, Tr$y), to = Tr, label = "treated",      label_cex = 1.0)
  flowx(from = pt(SPINE, D$y),  to = D,  label = "untreated",    label_cex = 1.0)
  flowx(from = pt(SPINE, A$y),  to = A,  label = "asymptomatic", label_cex = 1.0)

  ## --- recovery ---
  flowy(from = D, to = A, label = expression(r[D]), label_cex = 1.2)
  flowx(from = A, to = U, label = expression(r[A]), label_cex = 1.2)

  ## --- returns to S, four separate lanes ---
  path(c(Tr$x, Tr$x, Ph$x, Ph$x), c(Tr$y1, 4.62, 4.62, Ph$y1))
  lab(2.36, 4.78, expression(r[T]), cex = 1.2)
  flowy(from = Ph, to = S, label = expression(r[P]),      label_cex = 1.2)
  flowy(from = Pc, to = S, label = expression(r[P[c]]),   label_cex = 1.2)
  path(c(U$x1, 6.15, 6.15, 0.30, 0.30, S$x0),
       c(U$y, U$y, LANE_U, LANE_U, S$y, S$y))
  lab(5.92, LANE_U + 0.19, expression(r[U]), cex = 1.2)

  ## --- chemoprevention: a discrete jump, out of S, U, A, D and T ---
  path(c(4.10, 4.10, Pc$x, Pc$x), c(2.04, LANE_CP, LANE_CP, Pc$y0),
       col = PULSE, lty = 2, lwd = 1.6)
  lab(3.45, 0.61,
      expression(atop(italic("MDA / SMC / PMC"),
                      italic("from S, U, A, D, T"))),
      col = PULSE, cex = 0.92)

  ## ======================= MOSQUITO ======================================
  E  <- nd(0.90, MY, expression(bold(E)),    col_aq)
  L  <- nd(2.20, MY, expression(bold(L)),    col_aq)
  Pl <- nd(3.50, MY, expression(bold(P[L])), col_aq)
  Sm <- nd(4.80, MY, expression(bold(S[M])), col_mosq)
  Em <- nd(6.10, MY, expression(bold(E[M])), col_mosq, heavy = TRUE)
  Im <- nd(7.75, MY, expression(bold(I[M])), col_mosq, heavy = TRUE)

  bracket(4.44, 8.12, MY - 0.58, MY + 0.30)          # adult pool

  flowx(from = E,  to = L,  label = expression(1/d[E]), label_cex = 1.2)
  flowx(from = L,  to = Pl, label = expression(1/d[L]), label_cex = 1.2)
  flowx(from = Pl, to = Sm, label = expression(paste(frac(1, 2), " ", 1/d[P])),
        label_cex = 1.1)
  flowx(from = Sm, to = Em, label = expression(Lambda[s]^M), label_cex = 1.2)
  flowx(from = Em, to = Im, label = "survives EIP", label_cex = 0.88)

  for (b in list(E, L, Pl, Sm, Em, Im)) death(b$x, b$y0)
  lab(0.36, MY - 0.78,
      expression(paste("death; larvae ", symbol("\265"), " ", n[L]/K[s](t),
                       "  (seasonality)")), adj = c(0, 0.5), cex = 0.92)
  lab(8.18, MY - 0.34, expression(mu[s](t)), adj = c(0, 0.5), cex = 1.1)

  ## oviposition: from ALL adults, entering E from the left
  path(c(5.45, 5.45, 0.35, 0.35, E$x0),
       c(MY - 0.58, MY - 1.18, MY - 1.18, MY, MY))
  lab(2.90, MY - 1.03, expression(beta[s]), cex = 1.2)

  ## ======================= COUPLING ======================================
  ## dashed: these set a rate in the other panel, they are not flows
  path(c(Im$x, Im$x, SPINE, SPINE), c(Im$y1, LANE_EIR, LANE_EIR, FOOT),
       col = COUPLE, lty = 2, lwd = 1.8)
  lab(6.00, LANE_EIR + 0.22, expression(paste("EIR, lagged ", tau[E])),
      col = COUPLE, cex = 1.0)

  path(c(Sm$x, Sm$x), c(2.04, Sm$y1), col = COUPLE, lty = 2, lwd = 1.8)
  lab(Sm$x - 0.16, 1.68, expression(Lambda[s]^M), col = COUPLE, cex = 1.1,
      adj = c(1, 0.5))

  ## ======================= ANNOTATION ====================================
  lab(0.28, 4.78, "HUMAN", adj = c(0, 0.5), font = 2, cex = 1.3)
  lab(0.28, MY + 0.62, "MOSQUITO", adj = c(0, 0.5), font = 2, cex = 1.3)

  list(x0 = 0.20, x1 = 8.62, y0 = -1.72, y1 = 5.00)
}

flodia_png(model_flow, filepath = "vignettes/model_flow.png",
           width = 2400, res = 250)
file.copy("vignettes/model_flow.png", "man/figures/model_flow.png",
          overwrite = TRUE)
message("wrote vignettes/model_flow.png and man/figures/model_flow.png")
