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
# INVARIANTS. Break one of these and the figure states something false.
#
# SCOPE  Antimalarial resistance is out of scope, so the treated state is a
#        single T. The model splits it into T/T_slow when resistance is on;
#        vignette appendix B.4 documents that and section 2 discloses it.
#
# FIDELITY
#   A's only infection outflow is h_c (odin :255), so A's at-risk feeder joins
#   the spine ABOVE the asymptomatic branch. U loses the full FOI (:261), so U's
#   feeder joins below it. Getting this backwards draws a phantom A -> A flow.
#   Chemoprevention clears S, U, A, D, T -- S is outside the infectious bracket,
#   so the pulse must show an S source of its own.
#   Larval death is affine, me*(1 + nL/K), not proportional; pupal death is
#   density-independent. E is EARLY LARVAL, not eggs. E_M is exposed, not
#   infectious. Both couplings are lagged.
#
# GRAMMAR is three-part, not two:
#   solid black + filled triangle = a rate in the ODE, people move
#   dashed purple                 = a scalar coupling, nobody moves
#   dashed red                    = people move, instantaneously (a state jump)
#   dotted grey                   = a grouping, never a compartment; always labelled
#   Set membership is never drawn in flow ink, and there is ONE arrowhead glyph
#   for flows -- flodia's filled triangle, which hand-drawn paths must match.
#
# GEOMETRY  No line through a box. No two flows sharing a lane. Every arrowhead
#   lands on something. All paths orthogonal. Brackets are DERIVED from the boxes
#   they enclose and drawn BEFORE any node, so they never tint a fill.
#
# RENDERING  Authored at width 2400, res 250 (9.6 in). To check it at the size
#   readers actually get, render at width 700, res 73 (9.59 in) -- NOT
#   flodia_png's default res, which gives 3.5 in and inflates every cex ~2.7x.
#
# KNOWN, DEFERRED  The fill palette does not survive greyscale (six of seven
#   fills sit inside a 27-level luminance band) and S/P/P_c collapse under
#   deuteranopia. Fixing that means abandoning flodia's light_palette() for a
#   luminance ladder. Human ageing/death/births are disclosed in prose rather
#   than drawn, to protect legibility at 520 px.
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

COUPLE <- "#464273"   # scalar coupling      (9.2:1 on white)
PULSE  <- "#B03A22"   # discrete state jump  (6.0:1)
GREY   <- "#666666"   # annotation text      (5.7:1)
POOL   <- "#5F5563"   # groupings            (7.0:1, darker than before so the
                      # label sitting in the tint clears AA comfortably)
INK    <- "grey15"

RX <- 0.30; RY <- 0.20
CEX_NODE <- 1.45; CEX_RATE <- 1.20; CEX_NOTE <- 1.10

## ---- helpers ---------------------------------------------------------------
nd <- function(x, y, label, col, heavy = FALSE) {
  node(x = x, y = y, rx = RX, ry = RY, label = label, node_col = col,
       border_col = INK, label_cex = CEX_NODE, lwd = if (heavy) 2.6 else 1.0)
}

## a bare coordinate pair. NOT flodia::node(r = 0), which draws a degenerate
## rect in flodia's default grey. flowx/flowy read only these six fields.
pt <- function(x, y) list(x = x, y = y, x0 = x, x1 = x, y0 = y, y1 = y)

## the one arrowhead glyph, matching flodia's own flows
head_at <- function(x0, y0, x1, y1, col) {
  shape::Arrowhead(x1, y1, angle = atan2(y1 - y0, x1 - x0) * 180 / pi,
                   arr.type = "triangle", arr.length = 0.15, arr.width = 0.1,
                   lcol = col, arr.col = col)
}

## orthogonal polyline; filled-triangle head on the last segment unless head = FALSE
path <- function(xs, ys, col = INK, lty = 1, lwd = 1.2, head = TRUE) {
  stopifnot(length(xs) == length(ys), length(xs) >= 2)
  n <- length(xs)
  segments(xs[-n], ys[-n], xs[-1], ys[-1], col = col, lty = lty, lwd = lwd)
  if (head) head_at(xs[n - 1], ys[n - 1], xs[n], ys[n], col)
}

## a bracket DERIVED from the boxes it encloses -- never eyeballed
hull <- function(x, y, pad = 0.09) {
  list(x0 = min(x) - RX - pad, x1 = max(x) + RX + pad,
       y0 = min(y) - RY - pad, y1 = max(y) + RY + pad)
}
bracket <- function(h) {
  rect(h$x0, h$y0, h$x1, h$y1, border = POOL, lty = 3, lwd = 1.8, col = NA)
}

death <- function(x, y0, len = 0.34) path(c(x, x), c(y0, y0 - len))

lab <- function(x, y, txt, col = GREY, cex = CEX_NOTE, adj = c(0.5, 0.5),
                font = 1) {
  text(x, y, txt, col = col, cex = cex, adj = adj, font = font)
}

## flow-label defaults. flodia otherwise uses pure-black italic, which reads as
## a different class of symbol from every hand-drawn label.
fl <- function(...) {
  a <- list(...)
  d <- list(label_col = INK, label_font = 1, label_cex = CEX_RATE,
            arr_col = INK)
  d[names(a)] <- a      # explicit args win, no duplicate-formal error
  d
}

model_flow <- function() {

  ## ---- layout table: every position named exactly once --------------------
  SPINE  <- 2.10
  COL    <- 3.60                                # T / D / A column
  UX     <- 5.25                                # U
  LEFT   <- 0.80                                # P / S / P_c column
  Y_T    <- 4.00; Y_D <- 3.15; Y_A <- 2.30      # right column rows
  FOOT_U <- 1.86                                # U joins below the h_a branch
  FOOT_A <- 2.72                                # A joins ABOVE it
  LANE_U <- 1.30; LANE_CP <- 0.94; LANE_EIR <- 0.28
  MY     <- -0.55
  MX     <- c(E = 0.80, L = 2.15, Pl = 3.50, Sm = 4.85, Em = 6.20, Im = 7.90)

  INFECT <- hull(c(COL, UX), c(Y_T, Y_A))
  ADULT  <- hull(MX[c("Sm", "Im")], c(MY, MY))

  ## brackets first, so they never tint a node fill
  bracket(INFECT)
  bracket(ADULT)

  ## ======================= HUMAN =========================================
  Ph <- nd(LEFT, Y_T,  expression(bold(P)),    col_proph)
  S  <- nd(LEFT, Y_D,  expression(bold(S)),    col_sus)
  Pc <- nd(LEFT, Y_A,  expression(bold(P[c])), col_proph)
  Tr <- nd(COL,  Y_T,  expression(bold(T)),    col_treat, heavy = TRUE)
  D  <- nd(COL,  Y_D,  expression(bold(D)),    col_clin,  heavy = TRUE)
  A  <- nd(COL,  Y_A,  expression(bold(A)),    col_asym,  heavy = TRUE)
  U  <- nd(UX,   Y_A,  expression(bold(U)),    col_asym,  heavy = TRUE)

  lab(INFECT$x1 + 0.14, Y_D + 0.50, "infectious to\nmosquitoes", col = POOL,
      cex = 1.0, font = 3, adj = c(0, 0.5))

  ## --- infection ----------------------------------------------------------
  do.call(flowx, fl(from = S, to = pt(SPINE, Y_D), label = ""))
  segments(SPINE, FOOT_U, SPINE, Y_T, col = INK, lwd = 1.2)      # the spine
  lab(SPINE - 0.14, Y_T + 0.26, expression(Lambda), col = INK, cex = CEX_RATE,
      adj = c(1, 0.5))

  ## at-risk membership, in POOL dotted so it is never read as a flow
  path(c(A$x, A$x, SPINE), c(A$y1, FOOT_A, FOOT_A), col = POOL, lty = 3,
       head = FALSE)
  path(c(U$x, U$x, SPINE), c(U$y0, FOOT_U, FOOT_U), col = POOL, lty = 3,
       head = FALSE)
  lab(SPINE + 0.16, FOOT_U - 0.26, expression(paste("at risk: ", S + A + U)),
      col = POOL, cex = 1.0, adj = c(0, 0.5))

  do.call(flowx, fl(from = pt(SPINE, Y_T), to = Tr, label = "treated"))
  do.call(flowx, fl(from = pt(SPINE, Y_D), to = D,  label = "untreated"))
  do.call(flowx, fl(from = pt(SPINE, Y_A), to = A,  label = "asymptomatic"))

  ## --- recovery -----------------------------------------------------------
  do.call(flowy, fl(from = D, to = A, label = expression(r[D])))
  do.call(flowx, fl(from = A, to = U, label = expression(r[A])))

  ## --- returns to S: four lanes, none shared ------------------------------
  path(c(Tr$x, Tr$x, Ph$x, Ph$x), c(Tr$y1, Y_T + 0.66, Y_T + 0.66, Ph$y1))
  lab(SPINE + 0.34, Y_T + 0.82, expression(r[T]), col = INK, cex = CEX_RATE)
  do.call(flowy, fl(from = Ph, to = S, label = expression(r[P])))
  do.call(flowy, fl(from = Pc, to = S, label = expression(r[P[c]])))
  path(c(U$x1, 6.20, 6.20, 0.30, 0.30, S$x0),
       c(U$y, U$y, LANE_U, LANE_U, S$y, S$y))
  lab(7.05, LANE_U - 0.22, expression(r[U]), col = INK, cex = CEX_RATE)

  ## --- chemoprevention: a state jump, from S AND the infectious states -----
  path(c(4.10, 4.10, Pc$x, Pc$x), c(INFECT$y0, LANE_CP, LANE_CP, Pc$y0),
       col = PULSE, lty = 4, lwd = 1.7)
  path(c(1.45, 1.45), c(S$y0, LANE_CP), col = PULSE, lty = 4, lwd = 1.7,
       head = FALSE)                                  # S is a source too
  lab(3.00, LANE_CP - 0.32, expression(italic("MDA / SMC / PMC")),
      col = PULSE, cex = CEX_NOTE)

  ## ======================= MOSQUITO ======================================
  E  <- nd(MX[["E"]],  MY, expression(bold(E)),    col_aq)
  L  <- nd(MX[["L"]],  MY, expression(bold(L)),    col_aq)
  Pl <- nd(MX[["Pl"]], MY, expression(bold(P[L])), col_aq)
  Sm <- nd(MX[["Sm"]], MY, expression(bold(S[M])), col_mosq)
  Em <- nd(MX[["Em"]], MY, expression(bold(E[M])), col_mosq)
  Im <- nd(MX[["Im"]], MY, expression(bold(I[M])), col_mosq, heavy = TRUE)

  lab(mean(MX[c("Sm", "Im")]), ADULT$y1 + 0.20, "adult females: all lay eggs",
      col = POOL, cex = 1.0, font = 3)

  do.call(flowx, fl(from = E,  to = L,  label = expression(1/d[E])))
  do.call(flowx, fl(from = L,  to = Pl, label = expression(1/d[L])))
  do.call(flowx, fl(from = Pl, to = Sm, label = expression(1/(2 * d[P])), label_cex = 1.0))
  do.call(flowx, fl(from = Sm, to = Em, label = expression(Lambda[s]^M)))
  do.call(flowx, fl(from = Em, to = Im, label = "survives EIP", label_cex = 1.0))

  for (b in list(E, L, Pl, Sm, Em, Im)) death(b$x, b$y0)
  lab(0.52, MY - 1.06,
      expression(paste("death; larval death rises with ", n[L]/K[s](t),
                       " -- seasonality")), adj = c(0, 0.5), cex = 1.0)
  lab(MX[["Im"]] + 0.48, MY - 0.36, expression(mu[s](t)), adj = c(0, 0.5),
      cex = CEX_RATE)

  ## oviposition: from ALL adults, entering E from the left
  path(c(5.52, 5.52, 0.30, 0.30, E$x0),
       c(ADULT$y0, MY - 1.50, MY - 1.50, MY, MY))
  lab(3.00, MY - 1.34, expression(beta[s]), col = INK, cex = CEX_RATE)

  ## ======================= COUPLING ======================================
  path(c(Im$x, Im$x, SPINE, SPINE), c(Im$y1, LANE_EIR, LANE_EIR, FOOT_U),
       col = COUPLE, lty = 2, lwd = 1.8)
  lab(6.10, LANE_EIR + 0.26,
      expression(paste("EIR = ", a[s](t), " ", I[M], ", lagged ", tau[E])),
      col = COUPLE, cex = CEX_NOTE)

  path(c(Sm$x, Sm$x), c(INFECT$y0, Sm$y1), col = COUPLE, lty = 2, lwd = 1.8)
  lab(Sm$x + 0.34, FOOT_U - 0.26,
      expression(paste(Lambda[s]^M, " = ", a[s](t), " infectivity, lagged ",
                       tau[l])), col = COUPLE, cex = CEX_NOTE, adj = c(0, 0.5))

  ## ======================= ANNOTATION ====================================
  lab(0.28, Y_T + 0.98, "HUMAN", adj = c(0, 0.5), font = 2, cex = 1.3)
  lab(0.28, MY + 0.88, "MOSQUITO", adj = c(0, 0.5), font = 2, cex = 1.3)
  lab(0.30, MY - 1.92,
      expression(italic("human boxes also age, die, and are replenished by births")),
      adj = c(0, 0.5), cex = 0.98)

  list(x0 = 0.20, x1 = 8.60, y0 = -2.42, y1 = 5.14)
}

## ---- write ------------------------------------------------------------------
stopifnot(file.exists("DESCRIPTION"))
tmp <- tempfile(fileext = ".png")
flodia_png(model_flow, filepath = tmp, width = 2400, res = 250)
stopifnot(file.copy(tmp, "vignettes/model_flow.png", overwrite = TRUE),
          file.copy(tmp, "man/figures/model_flow.png", overwrite = TRUE))
message("wrote vignettes/model_flow.png and man/figures/model_flow.png")
