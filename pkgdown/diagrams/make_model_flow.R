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
# DESIGN NOTES (see git history for the review that produced them)
#
# * The figure is displayed at ~520-700 px, not the 2400 px it is authored at.
#   Everything is sized for the small end: no explanatory legend blocks (they
#   live in the caption), no compound algebra on the arrows, and a canvas kept
#   narrow so the same cex renders larger.
# * Flows leave the set of compartments they actually leave. Three quantities
#   are sums over several compartments -- the at-risk pool, the infectious pool
#   and the adult mosquito pool -- and each is drawn as a labelled bracket that
#   the arrow departs from, never as a single representative box.
# * Two visual grammars, kept apart: solid black = a flow of individuals;
#   dashed colour = something that is not a flow (the chemoprevention state
#   jump, and the two transmission couplings, which set a rate rather than
#   moving anyone).
# * Compartment class is encoded by border weight/style as well as fill, so it
#   survives greyscale and colour-vision deficiency -- the standard the README
#   already claims for the comparison figures.
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

COUPLE <- "#464273"   # transmission coupling (darkened for contrast)
PULSE  <- "#B03A22"   # discrete chemoprevention pulse (>= 4.5:1 on white)
GREY   <- "#666666"   # all annotation text (5.7:1); one tone, not five
POOL   <- "#7A6E7E"   # bracket / grouping rules

## Border encoding is redundant with fill, so class survives greyscale and CVD.
##   susceptible / prophylaxis : thin solid      infectious : heavy solid
##   treated                   : dashed          aquatic    : thin solid
BORDER <- list(plain = list(lty = 1, lwd = 1.0),
               inf   = list(lty = 1, lwd = 2.6),
               treat = list(lty = 5, lwd = 2.0))

nd <- function(x, y, label, col, cls = "plain", rx = 0.30, ry = 0.20) {
  b <- BORDER[[cls]]
  node(x = x, y = y, rx = rx, ry = ry, label = label, node_col = col,
       border_col = "grey15", label_cex = 1.45, lty = b$lty, lwd = b$lwd)
}
pt <- function(x, y) node(x = x, y = y, r = 0)

## a labelled bracket: the device for "this flow leaves a SET of compartments"
bracket <- function(x0, x1, y0, y1, label, lx, ly, cex = 0.95) {
  rect(x0, y0, x1, y1, border = POOL, lty = 3, lwd = 1.8,
       col = adjustcolor(POOL, alpha.f = 0.05))
  text(lx, ly, label, col = POOL, font = 3, cex = cex, adj = c(0.5, 0.5))
}

## a death outflow: short stub to nowhere, the anchor vector control acts on
death <- function(x, y0, len = 0.42, label = NULL, lx = NULL, ly = NULL) {
  arrows(x, y0, x, y0 - len, length = 0.10, angle = 20, col = "grey15", lwd = 1.1)
  if (!is.null(label)) text(lx, ly, label, col = GREY, cex = 0.92, adj = c(0, 0.5))
}

model_flow <- function() {

  SPINE <- 2.10
  COL   <- 3.55

  ## ======================= HUMAN =========================================
  ## infectious pool: exactly the states entering `inf` in the odin
  bracket(3.24, 5.94, 1.72, 5.02, "infectious to mosquitoes", 4.78, 3.02)

  S  <- nd(0.80, 3.55, expression(bold(S)),    col_sus,   "plain")
  Ts <- nd(COL,  4.70, expression(bold(T[s])), col_treat, "treat")
  Tr <- nd(COL,  3.90, expression(bold(T)),    col_treat, "treat")
  D  <- nd(COL,  3.10, expression(bold(D)),    col_clin,  "inf")
  A  <- nd(COL,  2.30, expression(bold(A)),    col_asym,  "inf")
  U  <- nd(5.15, 2.30, expression(bold(U)),    col_asym,  "inf")
  Ph <- nd(6.85, 4.30, expression(bold(P)),    col_proph, "plain")
  Pc <- nd(0.80, 1.75, expression(bold(P[c])), col_proph, "plain")

  ## --- infection: acts on the at-risk pool S + A + U ---
  FOOT <- 1.98
  flowx(from = S, to = pt(SPINE, S$y), label = expression(Lambda),
        label_pos = 0.62, label_cex = 1.15)
  flowy(from = pt(SPINE, FOOT), to = pt(SPINE, Ts$y), arr_width = 0)
  ## A and U are at risk too. One feeder lane along the foot of the spine,
  ## tapped from both, so the four splits are visibly fed by all three.
  segments(A$x, A$y0, A$x, FOOT, col = "grey15", lwd = 1)
  segments(U$x, U$y0, U$x, FOOT, col = "grey15", lwd = 1)
  arrows(U$x, FOOT, SPINE, FOOT, length = 0.10, angle = 20,
         col = "grey15", lwd = 1)
  text(SPINE + 0.16, FOOT - 0.46,
       expression(paste(Lambda, " acts on ", S + A + U)),
       col = GREY, cex = 0.95, adj = c(0, 0.5))

  flowx(from = pt(SPINE, Ts$y), to = Ts, label = "treated, slow",
        label_cex = 0.98)
  flowx(from = pt(SPINE, Tr$y), to = Tr, label = "treated", label_cex = 0.98)
  flowx(from = pt(SPINE, D$y),  to = D,  label = "untreated", label_cex = 0.98)
  flowx(from = pt(SPINE, A$y),  to = A,  label = "asymptomatic",
        label_cex = 0.98)

  ## --- recovery ---
  flowy(from = D, to = A, label = expression(r[D]), label_cex = 1.15)
  flowx(from = A, to = U, label = expression(r[A]), label_cex = 1.15)
  ## T_s into P's LEFT face, T into its bottom: one flow per side of P
  turnx(from = Ts, mid_x = 6.30, to = Ph, label = expression(r[T]*"(slow)"),
        label_x = 6.02, label_y = 4.86, label_cex = 1.15)
  bendx(from = Tr, to = Ph, label_to = expression(r[T]), label_cex = 1.15)

  ## --- returns to S: top, bottom, left. No two share a lane. ---
  turny(from = Ph, mid_y = Ts$y1 + 0.36, to = S, label = expression(r[P]),
        label_cex = 1.15)
  turny(from = U,  mid_y = 1.30, to = S, label = expression(r[U]),
        label_x = 4.30, label_y = 1.44, label_cex = 1.15)
  turnx(from = Pc, mid_x = 0.30, to = S, label = expression(r[P[c]]),
        label_x = 0.52, label_y = 2.62, label_cex = 1.15)

  ## --- chemoprevention: a discrete jump, out of S/U/A/D/T ---
  segments(1.10, 0.92, 4.10, 0.92, col = PULSE, lty = 2, lwd = 1.6)
  segments(4.10, 0.92, 4.10, 1.72, col = PULSE, lty = 2, lwd = 1.6)
  arrows(1.10, 0.92, Pc$x, Pc$y0, length = 0.12, angle = 20, col = PULSE,
         lty = 2, lwd = 1.6)
  text(2.55, 0.72, expression(italic("MDA / SMC / PMC: a fraction of S, U, A, D, T")),
       col = PULSE, cex = 0.95, adj = c(0.5, 0.5))

  ## ======================= MOSQUITO ======================================
  MY <- -0.30
  E  <- nd(0.80, MY, expression(bold(E)),    col_aq,   "plain")
  L  <- nd(2.10, MY, expression(bold(L)),    col_aq,   "plain")
  Pl <- nd(3.40, MY, expression(bold(P[L])), col_aq,   "plain")
  Sm <- nd(4.70, MY, expression(bold(S[M])), col_mosq, "plain")
  Em <- nd(6.00, MY, expression(bold(E[M])), col_mosq, "inf")
  Im <- nd(7.30, MY, expression(bold(I[M])), col_mosq, "inf")

  bracket(4.34, 7.66, MY - 0.30, MY + 0.34, "", 6.00, MY)

  flowx(from = E,  to = L,  label = expression(1/d[E]), label_cex = 1.15)
  flowx(from = L,  to = Pl, label = expression(1/d[L]), label_cex = 1.15)
  flowx(from = Pl, to = Sm, label = expression(paste(frac(1, 2), " ", 1/d[P])),
        label_cex = 1.05)
  flowx(from = Sm, to = Em, label = expression(Lambda[s]^M), label_cex = 1.15)
  flowx(from = Em, to = Im, label = "survives EIP", label_cex = 0.98)

  ## oviposition leaves ALL adults, not S_M
  turny(from = pt(6.00, MY - 0.30), mid_y = MY - 0.80, to = E,
        label = expression(beta[s]), label_cex = 1.15)

  ## deaths: the anchor bed nets, IRS and seasonality act on
  for (b in list(E, L, Pl)) death(b$x, b$y0, 0.30)
  text(0.28, MY - 1.04,
       expression(paste("death; larvae ", symbol("\265"), " ", n[L]/K[s](t),
                        "  (seasonality)")),
       col = GREY, cex = 0.92, adj = c(0, 0.5))
  for (b in list(Sm, Em, Im)) death(b$x, b$y0 - 0.30, 0.28)
  text(7.66, MY - 0.72, expression(mu[s](t)), col = GREY, cex = 1.05,
       adj = c(0, 0.5))

  ## ======================= COUPLING ======================================
  ## dashed: these set a rate in the other panel, they are not flows
  segments(Im$x, Im$y1 + 0.34, Im$x, 1.05, col = COUPLE, lty = 2, lwd = 1.8)
  segments(Im$x, 1.05, SPINE + 0.55, 1.05, col = COUPLE, lty = 2, lwd = 1.8)
  arrows(SPINE + 0.55, 1.05, SPINE + 0.10, 1.05, length = 0.11, angle = 20,
         col = COUPLE, lty = 2, lwd = 1.8)
  text(6.42, 1.28, expression(paste("EIR, lagged ", tau[E])),
       col = COUPLE, cex = 1.0, adj = c(0, 0.5))

  segments(4.70, 1.72, 4.70, 0.62, col = COUPLE, lty = 2, lwd = 1.8)
  arrows(4.70, 0.62, 4.70, MY + 0.40, length = 0.11, angle = 20,
         col = COUPLE, lty = 2, lwd = 1.8)

  ## ======================= ANNOTATION ====================================
  text(0.28, 5.30, "HUMAN",    adj = c(0, 0.5), font = 2, cex = 1.3, col = GREY)
  text(0.28, MY + 0.78, "MOSQUITO", adj = c(0, 0.5), font = 2, cex = 1.3,
       col = GREY)

  list(x0 = 0.20, x1 = 7.95, y0 = -1.55, y1 = 5.45)
}

flodia_png(model_flow, filepath = "vignettes/model_flow.png",
           width = 2400, res = 250)
file.copy("vignettes/model_flow.png", "man/figures/model_flow.png",
          overwrite = TRUE)
message("wrote vignettes/model_flow.png and man/figures/model_flow.png")
