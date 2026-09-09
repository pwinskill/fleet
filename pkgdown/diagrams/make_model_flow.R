# Generate the blink model flow diagram used by the README and vignette("model").
#
#   Rscript pkgdown/diagrams/make_model_flow.R
#
# Needs flodia (GitHub only, deliberately NOT a package dependency):
#   remotes::install_github("mrc-ide/flodia")
#
# Outputs vignettes/model_flow.png and man/figures/model_flow.png, both
# committed, so neither the vignette nor the README needs flodia at build time.
# Set BLINK_DIAGRAM_OUT to write elsewhere, plus a 700-px proof (development).
#
# ---------------------------------------------------------------------------
# INVARIANTS. Break one of these and the figure states something false.
#
# SCOPE  Antimalarial resistance is out of scope, so the treated state is a
#        single T. The model splits it into T/T_slow when resistance is on;
#        vignette appendix B.4 documents that and section 2 discloses it.
#
# FIDELITY
#   The infection hazard acts on S + A + U. A's ONLY infection outflow is the
#   -h_c[i,j]*A[i,j] term in deriv(A) -- there is no A -> A arrow -- and U loses
#   the full FOI in deriv(U). (Cite the equation, not a line number: bare line
#   numbers into a file that is actively edited rot, and these already had.)
#   That membership is stated in WORDS beside the hazard, never in ink: a line
#   from A back to the spine would read as a phantom A -> A flow, and the dotted
#   rings this used to carry cost a whole grammar element for one fact.
#   Chemoprevention clears S, U, A, D, T, P and P_c -- S is outside the
#   infectious pool, so the pulse must show an S source of its own.
#   Larval death is affine, me*(1 + n_L/K), not proportional; pupal death is
#   density-independent. E is EARLY LARVAL, not eggs. E_M is exposed, not
#   infectious, and the EIP chain is loss-free with survival applied at exit.
#   ALL adult females lay eggs, so oviposition leaves the adult grouping.
#   Both couplings are lagged. P and P_c are Erlang CHAINS, not single
#   compartments (NEWS 2026-09: k_P = 1 for AL, 16 SP-AQ; k_Pc = 14 SP-AQ).
#
# GRAMMAR is four-part, and the key states all of it:
#   solid black + filled triangle = a rate in the ODE, people move
#   dashed indigo                 = a scalar coupling, nobody moves
#   dashed red                    = people move, instantaneously (a state jump)
#   dotted grey                   = a grouping, never a compartment; always
#                                   labelled
#   offset layers                 = the array dimensions, and ONLY those
#   Set membership is never drawn in flow ink, and there is ONE arrowhead glyph
#   for flows -- flodia's filled triangle, which hand-drawn paths must match.
#
# GEOMETRY  No line through a box. No two flows sharing a lane. Every arrowhead
#   lands on something. All paths orthogonal. Flow labels sit BESIDE their line
#   (LGAP clear of it), never on it. Panels, decks and regions are DERIVED from
#   the boxes they enclose and drawn BEFORE any node, so they never tint a
#   fill. The transmission cycle runs counter-clockwise -- humans left to
#   right, mosquitoes right to left -- so the two couplings never cross each
#   other. The only crossings in the figure are where the EIR coupling meets
#   the two return lanes; both are drawn as hops (a gap in the EIR line).
#
# DIMENSIONS  The layered idiom appears TWICE and means one thing both times:
#   the array dimensions. A deck behind the human panel, its two visible edges
#   labelled one dimension each (age i, biting heterogeneity j), and a deck
#   behind the mosquito panel (species s). It used to also sit behind P, P_c and
#   the EIP to mean "Erlang chain" -- one idiom carrying two unrelated meanings,
#   which is what made the figure hard to read. The chains are now disclosed in
#   the caption and vignette("model") instead, where the stage counts belong
#   anyway (they vary by drug). Immunity likewise: it is algebraic, not a state,
#   so it is named on the hazard rather than boxed in the busiest part of the
#   panel.
#
# RENDERING  Authored at width 2400, res 250 (9.6 in). Check it at the size
#   readers actually get by rendering at width 700, res 73 (9.59 in) -- NOT
#   flodia_png's default res, which gives 3.5 in and inflates every cex ~2.7x.
#   The fills are a luminance ladder spanning 36 levels (flodia's light_palette
#   packs seven fills into 27): every SPATIALLY ADJACENT pair differs by at least
#   13 levels, so the figure survives greyscale, and every state is named too --
#   colour never carries meaning alone.
# ---------------------------------------------------------------------------

library(flodia)

## ---- palette: a luminance ladder -------------------------------------------
SUS   <- "#EAF0F9"   # susceptible          luminance 94
SUBP  <- "#F2DFAF"   # sub-patent                      88
TREAT <- "#BCDCBC"   # treated                         83
PROPH <- "#CBBCE2"   # prophylaxis                     76
AQ    <- "#E0D9CB"   # aquatic mosquito                85
ASYM  <- "#E3BC6A"   # asymptomatic                    75
MOSQ  <- "#8FBAC5"   # adult mosquito                  70
CLIN  <- "#CE8484"   # clinical disease                58

PANEL  <- "#FFFFFF"
DECK1  <- "#F1F1EE"
DECK2  <- "#E4E4E0"
POOLT  <- "#F7F4EF"

COUPLE <- "#463C8A"   # scalar coupling      (8.6:1 on white)
PULSE  <- "#B03A22"   # discrete state jump  (6.0:1)
GREY   <- "#5B5B5B"   # annotation text      (7.0:1)
POOL   <- "#6B6B6B"   # groupings            (5.6:1)
INK    <- "grey15"

RX <- 0.34; RY <- 0.26                        # box half-widths
CEX_NODE <- 1.5; CEX_RATE <- 1.15; CEX_NOTE <- 1.0; CEX_SMALL <- 0.85
LGAP <- 0.19          # every flow label sits this far off its line

## ---- helpers ---------------------------------------------------------------
nd <- function(x, y, label, col, heavy = FALSE) {
  node(x = x, y = y, rx = RX, ry = RY, label = label, node_col = col,
       border_col = INK, label_cex = CEX_NODE, lwd = if (heavy) 2.4 else 1.0)
}

## the one arrowhead glyph, matching flodia's own flows
head_at <- function(x0, y0, x1, y1, col) {
  shape::Arrowhead(x1, y1, angle = atan2(y1 - y0, x1 - x0) * 180 / pi,
                   arr.type = "triangle", arr.length = 0.15, arr.width = 0.1,
                   lcol = col, arr.col = col)
}

## orthogonal polyline; filled-triangle head on the last segment unless head = FALSE
path <- function(xs, ys, col = INK, lty = 1, lwd = 1.3, head = TRUE) {
  stopifnot(length(xs) == length(ys), length(xs) >= 2)
  n <- length(xs)
  segments(xs[-n], ys[-n], xs[-1], ys[-1], col = col, lty = lty, lwd = lwd)
  if (head) head_at(xs[n - 1], ys[n - 1], xs[n], ys[n], col)
}

## a vertical that HOPS over the lanes it must cross (circuit-diagram convention)
vhop <- function(x, ya, yb, hops, col, lty = 1, lwd = 1.3, gap = 0.10) {
  lo <- min(ya, yb); hi <- max(ya, yb)
  cuts <- sort(unlist(lapply(hops, function(h) c(h - gap, h + gap))))
  cuts <- cuts[cuts > lo & cuts < hi]
  pts <- c(lo, cuts, hi)
  for (k in seq(1, length(pts) - 1, by = 2))
    segments(x, pts[k], x, pts[k + 1], col = col, lty = lty, lwd = lwd)
}

## a region DERIVED from the boxes it encloses -- never eyeballed
hull <- function(x, y, pad = 0.12) {
  list(x0 = min(x) - RX - pad, x1 = max(x) + RX + pad,
       y0 = min(y) - RY - pad, y1 = max(y) + RY + pad)
}
region <- function(h, fill = POOLT) {
  rect(h$x0, h$y0, h$x1, h$y1, border = POOL, lty = 3, lwd = 1.6, col = fill)
}

## a deck of offset panels behind a panel: the ONLY layered idiom in the figure,
## and it means one thing -- these boxes are arrays over the grid named on it
deck <- function(p, n = 2, dx = 0.20, dy = 0.20) {
  cols <- c(DECK2, DECK1)
  for (k in seq(n, 1))
    rect(p$x0 + k * dx, p$y0 + k * dy, p$x1 + k * dx, p$y1 + k * dy,
         col = cols[n - k + 1], border = POOL, lwd = 1.1)
  rect(p$x0, p$y0, p$x1, p$y1, col = PANEL, border = POOL, lwd = 1.4)
}

death <- function(x, y0, len = 0.44) path(c(x, x), c(y0, y0 - len))

lab <- function(x, y, txt, col = GREY, cex = CEX_NOTE, adj = c(0.5, 0.5),
                font = 1) {
  text(x, y, txt, col = col, cex = cex, adj = adj, font = font)
}

model_flow <- function() {

  ## ---- layout table: every position named exactly once --------------------
  ## humans, left to right
  LEFT  <- 1.25                                  # P / S / P_c column
  SPINE <- 2.95                                  # the infection spine
  MID   <- 4.85                                  # T / D / A column
  UX    <- 6.75                                  # U
  Y_T <- 4.70; Y_S <- 3.35; Y_A <- 2.00          # rows (P on Y_T, P_c on Y_A)
  LANE_RT <- 5.36                                # T -> P, above the boxes
  LANE_CP <- 1.28                                # chemoprevention pulse
  LANE_U  <- 0.92                                # U -> S
  RET_X   <- 0.42                                # the U -> S riser
  CP_X    <- 5.20                                # the pulse drop off the pool
  PULSE_S <- 1.45                                # the pulse's own S source
  PC_X    <- 1.10                                # P_c -> S, clear of the pulse
  HP <- list(x0 = 0.30, x1 = 7.42, y0 = 0.55, y1 = 5.90)

  ## mosquitoes, right to left, so the transmission cycle closes without the
  ## two couplings crossing each other
  MY <- -1.75
  MX <- c(Im = 0.95, Em = 2.75, Sm = 4.35, Pl = 6.00, L = 7.30, E = 8.60)
  LANE_OVI <- -3.05; OVI_X <- 9.15
  MP <- list(x0 = 0.24, x1 = 9.32, y0 = -3.30, y1 = -0.98)

  ## couplings
  LANE_EIR <- -0.90; LANE_FOIM <- -0.55; FOIM_X <- 7.85

  POOLH  <- hull(c(MID, UX), c(Y_A, Y_T))                  # infectious to mosquitoes
  ADULTH <- hull(MX[c("Im", "Sm")], c(MY, MY), pad = 0.20) # adult females

  ## ---- backdrop: decks first, so nothing tints a fill ---------------------
  deck(HP)                                   # age i x heterogeneity j
  deck(MP, dy = -0.20)                       # species s
  region(POOLH)
  region(ADULTH)

  ## ======================= HUMANS ==========================================
  Ph <- nd(LEFT, Y_T, expression(bold(P)),    PROPH)
  S  <- nd(LEFT, Y_S, expression(bold(S)),    SUS)
  Pc <- nd(LEFT, Y_A, expression(bold(P[c])), PROPH)
  Tr <- nd(MID,  Y_T, expression(bold(T)),    TREAT, heavy = TRUE)
  D  <- nd(MID,  Y_S, expression(bold(D)),    CLIN,  heavy = TRUE)
  A  <- nd(MID,  Y_A, expression(bold(A)),    ASYM,  heavy = TRUE)
  U  <- nd(UX,   Y_A, expression(bold(U)),    SUBP,  heavy = TRUE)

  lab(5.10, POOLH$y1 + 0.24, "infectious to mosquitoes",
      col = POOL, cex = CEX_NOTE, font = 3, adj = c(0, 0.5))

  ## --- infection: one spine, three branches -------------------------------
  segments(S$x1, Y_S, SPINE, Y_S, col = INK, lwd = 1.3)  # the at-risk pool feeds it
  segments(SPINE, Y_A, SPINE, Y_T, col = INK, lwd = 1.3)       # the spine
  lab(SPINE + 0.16, 4.28, expression(Lambda[ij]), col = INK,
      cex = CEX_RATE + 0.2, adj = c(0, 0.5))
  lab(SPINE + 0.16, 4.02, "infection hazard", col = GREY,
      cex = CEX_SMALL, adj = c(0, 0.5))
  lab(SPINE + 0.16, 3.80, "on S, A, U", col = GREY,
      cex = CEX_SMALL, adj = c(0, 0.5))

  br <- function(y, txt, to_x) {
    path(c(SPINE, to_x), c(y, y))
    lab((SPINE + to_x) / 2, y + LGAP, txt, col = INK, cex = CEX_RATE)
  }
  br(Y_T, "treated",     MID - RX)
  br(Y_S, "untreated",   MID - RX)
  br(Y_A, "no symptoms", MID - RX)

  ## --- recovery -----------------------------------------------------------
  path(c(MID, MID), c(D$y0, A$y1))
  lab(MID + LGAP, (D$y0 + A$y1) / 2, expression(r[D]), col = INK,
      cex = CEX_RATE, adj = c(0, 0.5))
  path(c(A$x1, U$x0), c(Y_A, Y_A))
  lab((A$x1 + U$x0) / 2, Y_A + LGAP, expression(r[A]), col = INK, cex = CEX_RATE)

  ## --- returns to S: four lanes, none shared ------------------------------
  path(c(MID, MID, LEFT, LEFT), c(Tr$y1, LANE_RT, LANE_RT, Ph$y1))
  lab(3.05, LANE_RT + 0.16, expression(r[T]), col = INK, cex = CEX_RATE)
  path(c(LEFT, LEFT), c(Ph$y0, S$y1))
  lab(LEFT + LGAP, (Ph$y0 + S$y1) / 2, expression(r[P]), col = INK,
      cex = CEX_RATE, adj = c(0, 0.5))
  path(c(PC_X, PC_X), c(Pc$y1, S$y0))
  lab(PC_X - LGAP, (Pc$y1 + S$y0) / 2, expression(r[P[c]]), col = INK,
      cex = CEX_RATE, adj = c(1, 0.5))
  path(c(UX, UX, RET_X, RET_X, S$x0),
       c(U$y0, LANE_U, LANE_U, Y_S, Y_S))
  lab(UX - 0.14, LANE_U + 0.16, expression(r[U]), col = INK, cex = CEX_RATE,
      adj = c(1, 0.5))

  ## --- chemoprevention: a state jump, from S AND the infectious states -----
  path(c(CP_X, CP_X, LEFT, LEFT), c(POOLH$y0, LANE_CP, LANE_CP, Pc$y0),
       col = PULSE, lty = 4, lwd = 1.7)
  path(c(PULSE_S, PULSE_S), c(S$y0, Pc$y1), col = PULSE, lty = 4, lwd = 1.7)
  lab(3.30, LANE_CP + 0.18, expression(italic("MDA / SMC / PMC")),
      col = PULSE, cex = CEX_NOTE)

  ## --- what the layers mean, and what happens between them ----------------
  lab(HP$x0 + 0.14, HP$y1 - 0.16, "HUMANS", col = INK, cex = 1.25, font = 2,
      adj = c(0, 0.5))
  lab(1.60, HP$y1 - 0.16, "every box is an array over the layers behind it",
      col = GREY, cex = CEX_NOTE, font = 3, adj = c(0, 0.5))
  ## the deck is labelled on its own two edges, one dimension each, so the
  ## layering reads as the [age x heterogeneity] grid and nothing else
  lab(HP$x1 + 0.32, HP$y1 + 0.30, expression(paste("age  ", italic(i), " = 1 ... 52")),
      col = GREY, cex = CEX_SMALL, adj = c(1, 0.5))
  lab(HP$x1 + 0.12, HP$y1 + 0.10,
      expression(paste("biting heterogeneity  ", italic(j), " = 1 ... 5")),
      col = GREY, cex = CEX_SMALL, adj = c(1, 0.5))
  lab(HP$x0 + 0.30, HP$y0 + 0.16,
      expression(italic(paste("people age along ", italic(i), " (", r[i],
                              "), die from every box (", mu[i](t),
                              "), and are born into S"))),
      col = GREY, cex = CEX_SMALL, adj = c(0, 0.5))

  ## ======================= MOSQUITOES ======================================
  E  <- nd(MX[["E"]],  MY, expression(bold(E)),    AQ)
  L  <- nd(MX[["L"]],  MY, expression(bold(L)),    AQ)
  Pl <- nd(MX[["Pl"]], MY, expression(bold(P[L])), AQ)
  Sm <- nd(MX[["Sm"]], MY, expression(bold(S[M])), MOSQ)
  Em <- nd(MX[["Em"]], MY, expression(bold(E[M])), MOSQ)
  Im <- nd(MX[["Im"]], MY, expression(bold(I[M])), MOSQ, heavy = TRUE)

  fl <- function(from, to, txt, cex = CEX_RATE) {
    path(c(from$x0, to$x1), c(MY, MY))
    lab((from$x0 + to$x1) / 2, MY + LGAP, txt, col = INK, cex = cex)
  }
  fl(E,  L,  expression(1/d[E]))
  fl(L,  Pl, expression(1/d[L]))
  fl(Pl, Sm, expression(1/(2 * d[P])), cex = CEX_NOTE)
  fl(Sm, Em, expression(Lambda[s]^M))
  fl(Em, Im, "survives EIP", cex = 0.78)

  for (b in list(E, L, Pl, Sm, Em, Im)) death(b$x, b$y0)

  lab(MP$x0 + 0.14, MP$y1 - 0.14, "MOSQUITOES", col = INK, cex = 1.25,
      font = 2, adj = c(0, 0.5))
  lab(ADULTH$x1 + 0.14, MP$y1 - 0.14, "adult females -- all lay eggs",
      col = POOL, cex = CEX_NOTE, font = 3, adj = c(0, 0.5))
  lab(MP$x0 + 0.30, MP$y0 + 0.16,
      expression(italic(paste("death: adults ", mu[s](t), "; larvae rise with ",
                              n[L]/K[s](t), " -- seasonality enters here"))),
      col = GREY, cex = CEX_SMALL, adj = c(0, 0.5))
  lab(MP$x1 + 0.46, MP$y0 - 0.30,
      expression(paste("s = 1 ... ", n[v], " species")), col = GREY,
      cex = CEX_SMALL, adj = c(1, 0.5))

  ## oviposition: from ALL adult females, entering E from the right
  path(c(3.50, 3.50, OVI_X, OVI_X, E$x1),
       c(ADULTH$y0, LANE_OVI, LANE_OVI, MY, MY))
  lab(7.70, LANE_OVI + 0.18, expression(beta[s]), col = INK, cex = CEX_RATE)

  ## ======================= COUPLINGS =======================================
  ## EIR: infectious mosquitoes -> the human hazard. Hops the two return lanes.
  segments(Im$x, Im$y1, Im$x, LANE_EIR, col = COUPLE, lty = 2, lwd = 1.7)
  segments(Im$x, LANE_EIR, SPINE, LANE_EIR, col = COUPLE, lty = 2, lwd = 1.7)
  vhop(SPINE, LANE_EIR, Y_A, hops = c(LANE_U, LANE_CP), col = COUPLE,
       lty = 2, lwd = 1.7)
  head_at(SPINE, Y_A - 0.20, SPINE, Y_A, COUPLE)
  lab(1.20, LANE_EIR + 0.20,
      expression(paste("EIR = ", sum(a[s](t) * I[M], s), ",  lagged ", tau[E])),
      col = COUPLE, cex = CEX_NOTE, adj = c(0, 0.5))

  ## human infectivity -> the mosquito FOI
  path(c(POOLH$x1, FOIM_X, FOIM_X, MX[["Sm"]], MX[["Sm"]]),
       c(Y_S, Y_S, LANE_FOIM, LANE_FOIM, Sm$y1), col = COUPLE, lty = 2, lwd = 1.7)
  lab(6.10, LANE_FOIM + 0.20,
      expression(paste(Lambda[s]^M, " = ", a[s](t), " x infectivity,  lagged ", tau[l])),
      col = COUPLE, cex = CEX_NOTE)

  ## ======================= KEY =============================================
  KX <- 8.18; KT <- 8.80; ky <- 5.52
  lab(KX, ky + 0.44, "Key", col = INK, cex = 1.15, font = 2, adj = c(0, 0.5))
  krow <- function(y, txt, draw) {
    draw(y); lab(KT, y, txt, col = GREY, cex = CEX_SMALL, adj = c(0, 0.5))
  }
  krow(ky, "flow: a rate in the ODE",
       function(y) path(c(KX, KX + 0.48), c(y, y)))
  krow(ky - 0.46, "coupling: nobody moves",
       function(y) path(c(KX, KX + 0.48), c(y, y), col = COUPLE, lty = 2, lwd = 1.7))
  krow(ky - 0.92, "pulse: a state jump",
       function(y) path(c(KX, KX + 0.48), c(y, y), col = PULSE, lty = 4, lwd = 1.7))
  krow(ky - 1.38, "grouping, not a state",
       function(y) rect(KX, y - 0.13, KX + 0.48, y + 0.13, border = POOL,
                        lty = 3, lwd = 1.6, col = POOLT))
  ## the layered idiom appears once in the figure and means one thing
  krow(ky - 1.90, "layers: the array dimensions",
       function(y) {
         rect(KX + 0.16, y - 0.05, KX + 0.50, y + 0.15, col = DECK2, border = POOL)
         rect(KX + 0.08, y - 0.10, KX + 0.42, y + 0.10, col = DECK1, border = POOL)
         rect(KX, y - 0.15, KX + 0.34, y + 0.05, col = PANEL, border = INK)
       })

  list(x0 = 0.10, x1 = 10.85, y0 = -3.90, y1 = 6.45)
}

## ---- write ------------------------------------------------------------------
out <- Sys.getenv("BLINK_DIAGRAM_OUT", "")
tmp <- tempfile(fileext = ".png")
flodia_png(model_flow, filepath = tmp, width = 2400, res = 250)
if (nzchar(out)) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  stopifnot(file.copy(tmp, file.path(out, "model_flow.png"), overwrite = TRUE))
  flodia_png(model_flow, filepath = file.path(out, "model_flow_700.png"),
             width = 700, res = 73)
  message("wrote ", out)
} else {
  stopifnot(file.exists("DESCRIPTION"))
  stopifnot(file.copy(tmp, "vignettes/model_flow.png", overwrite = TRUE),
            file.copy(tmp, "man/figures/model_flow.png", overwrite = TRUE))
  message("wrote vignettes/model_flow.png and man/figures/model_flow.png")
}
