# Generate the blink model flow diagram used at the top of vignette("model").
#
#   Rscript pkgdown/diagrams/make_model_flow.R
#
# Needs flodia (GitHub only, deliberately NOT a package dependency):
#   remotes::install_github("mrc-ide/flodia")
#
# Output: vignettes/model_flow.png, committed to the repo so that building the
# vignette needs nothing beyond base R.
#
# Layout notes: S receives its three returns on three different sides (r_P from
# above, r_U from below, r_Pc from the left) so no two return lanes overlap. The
# two transmission-coupling arrows run in opposite directions between the panels
# and necessarily cross once.

library(flodia)

## ---- palette ---------------------------------------------------------------
col_sus   <- light_palette("bu")
col_clin  <- light_palette("rd")
col_asym  <- light_palette("oryl")
col_treat <- light_palette("gn")
col_proph <- light_palette("pu")
col_aq    <- light_palette("gybr")
col_mosq  <- light_palette("gnbu")

COUPLE <- mid_palette("bupu")   # transmission coupling
PULSE  <- mid_palette("rdor")   # discrete chemoprevention pulse

nd <- function(x, y, label, col, rx = 0.26, ry = 0.17) {
  node(x = x, y = y, rx = rx, ry = ry, label = label, node_col = col,
       border_col = "grey20", label_cex = 1.3)
}
pt <- function(x, y) node(x = x, y = y, r = 0)   # invisible junction

model_flow <- function() {

  SPINE <- 2.45     # x of the vertical infection spine
  COL   <- 3.60     # x of the infection-outcome column

  ## ======================= HUMAN =========================================
  S  <- nd(0.85, 3.10, expression(bold(S)),    col_sus)
  Ts <- nd(COL,  4.40, expression(bold(T[s])), col_treat)
  Tr <- nd(COL,  3.55, expression(bold(T)),    col_treat)
  D  <- nd(COL,  2.70, expression(bold(D)),    col_clin)
  A  <- nd(COL,  1.85, expression(bold(A)),    col_asym)
  U  <- nd(5.30, 1.85, expression(bold(U)),    col_asym)
  Ph <- nd(5.95, 3.975, expression(bold(P)),   col_proph)
  Pc <- nd(0.85, 0.95, expression(bold(P[c])), col_proph)

  SPINE_FOOT <- 1.23

  ## --- infection: S -> spine -> the four outcomes ---
  flowx(from = S, to = pt(SPINE, S$y), arr_width = 0,
        label = expression(Lambda), label_pos = 0.80, label_cex = 1.05)
  flowy(from = pt(SPINE, SPINE_FOOT), to = pt(SPINE, Ts$y), arr_width = 0)

  flowx(from = pt(SPINE, Ts$y), to = Ts,
        label = expression(f[t]^eff ~ SPC %.% h^c))
  flowx(from = pt(SPINE, Tr$y), to = Tr,
        label = expression(f[t]^eff * (1 - SPC) %.% h^c))
  flowx(from = pt(SPINE, D$y), to = D,
        label = expression((1 - f[t]^eff) %.% h^c))
  flowx(from = pt(SPINE, A$y), to = A, label = expression(h^a))

  ## --- recovery ---
  flowy(from = D, to = A, label = expression(r[D]))
  flowx(from = A, to = U, label = expression(r[A]))
  bendx(from = Tr, to = Ph, label_to = expression(r[T]))
  bendx(from = Ts, to = Ph, label_to = expression(r[T]^slow))

  ## --- returns to S, one per side ---
  turny(from = Ph, mid_y = Ts$y1 + 0.34, to = S, label = expression(r[P]))
  turny(from = U,  mid_y = 1.52,         to = S, label = expression(r[U]),
        label_x = 4.55, label_y = 1.66)
  flowy(from = Pc, to = S, label = expression(r[P[c]]))

  ## --- chemoprevention: discrete pulse, dashed and coloured ---
  bendy(from = A, to = Pc, pos_from = 0.50, arr_lty = 2, arr_col = PULSE,
        label_to = expression(italic("MDA / SMC / PMC pulse")),
        label_col = PULSE, label_to_x = 1.55, label_to_y = 1.28)

  ## ======================= MOSQUITO ======================================
  MY <- -0.35
  E  <- nd(0.85, MY, expression(bold(E)),    col_aq)
  L  <- nd(2.10, MY, expression(bold(L)),    col_aq)
  Pl <- nd(3.35, MY, expression(bold(P[L])), col_aq)
  Sm <- nd(4.75, MY, expression(bold(S[M])), col_mosq)
  Em <- nd(6.15, MY, expression(bold(E[M])), col_mosq, rx = 0.30)
  Im <- nd(7.45, MY, expression(bold(I[M])), col_mosq)

  flowx(from = E,  to = L,  label = expression(1/d[E]))
  flowx(from = L,  to = Pl, label = expression(1/d[L]))
  flowx(from = Pl, to = Sm, label = expression(1/(2 * d[P])))
  flowx(from = Sm, to = Em, label = expression(Lambda^M))
  flowx(from = Em, to = Im, label = expression(e^{-mu[s] * tau[M]}))
  turny(from = Sm, mid_y = Sm$y0 - 0.48, to = E, label = expression(beta[s]))

  ## ======================= COUPLING ======================================
  ## humans -> mosquitoes: a stub out of the infectious pool, then down to S_M
  TAP <- U$x1 + 0.45
  flowx(from = U, to = pt(TAP, U$y), arr_width = 0, arr_col = COUPLE)
  turny(from = pt(TAP, U$y), mid_y = 0.66, to = pt(Sm$x, Sm$y1),
        arr_col = COUPLE, label_col = COUPLE, label_cex = 1.05,
        label = expression(Lambda^M == a[s] %.% bar(inf)),
        label_x = 6.75, label_y = 1.30)

  ## mosquitoes -> humans: EIR drives the infection spine
  turny(from = pt(Im$x, Im$y1), mid_y = 0.22, to = pt(SPINE, SPINE_FOOT),
        arr_col = COUPLE, label_col = COUPLE, label_cex = 1.05,
        label = expression(EIR == sum(a[s] * I[M], s, )),
        label_x = 3.60, label_y = 0.42)

  ## ======================= ANNOTATION ====================================
  text(0.28, Ts$y1 + 0.34, "HUMAN", adj = c(0, 0.5), font = 2, cex = 1.2,
       col = "grey35")
  text(0.28, Sm$y0 - 0.48, "MOSQUITO", adj = c(0, 0.5), font = 2, cex = 1.2,
       col = "grey35")
  text(8.20, Sm$y0 - 1.05,
       "every human compartment also ages between age groups, and loses people to death",
       adj = c(1, 0.5), font = 3, cex = 0.86, col = "grey45")

  list(x0 = 0.20, x1 = 8.20, y0 = -1.40, y1 = 5.05)
}

## vignettes/ for the article, man/figures/ for the README (GitHub renders from
## there); both are committed so neither needs flodia at build time.
flodia_png(model_flow, filepath = "vignettes/model_flow.png",
           width = 2600, res = 230)
file.copy("vignettes/model_flow.png", "man/figures/model_flow.png",
          overwrite = TRUE)
message("wrote vignettes/model_flow.png and man/figures/model_flow.png")
