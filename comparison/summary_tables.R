# The numbers quoted in vignettes/comparison.Rmd, computed from the saved CSVs.
#
#   Rscript comparison/summary_tables.R          # writes comparison/data/tables.md
#   CMP_SMOKE=1 Rscript comparison/summary_tables.R
#
# Markdown tables + one-line statistics, so the article's figures and its prose
# come from the same data. Paste from tables.md; do not hand-edit numbers.

.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")
suppressMessages({library(dplyr); library(tidyr)})
ROOT <- "C:/Users/pwinskil/Documents/dev/blink2/blink"
source(file.path(ROOT, "comparison", "theme.R"))
SMOKE <- nzchar(Sys.getenv("CMP_SMOKE"))
DDIR  <- file.path(ROOT, "comparison", "data"); if (SMOKE) { DDIR <- file.path(DDIR, "smoke"); BURN_Y <- 1L }
rd <- function(part) read.csv(file.path(DDIR, paste0("rep_", part, ".csv")), stringsAsFactors = FALSE)
eq <- rd("eq"); age <- rd("age"); monthly <- rd("monthly"); doy <- rd("doy"); timing <- rd("timing")

out <- character()
say <- function(fmt, ...) out <<- c(out, if (...length()) sprintf(fmt, ...) else fmt)  # literal when no args
md_table <- function(df) {
  df[] <- lapply(df, as.character)
  say(paste0("| ", paste(names(df), collapse = " | "), " |"))
  say(paste0("|", paste(rep(" --- ", ncol(df)), collapse = "|"), "|"))
  for (i in seq_len(nrow(df))) say(paste0("| ", paste(df[i, ], collapse = " | "), " |"))
  say("")
}
q <- function(v, p) unname(quantile(v, p))
med_rng <- function(v, fmt = "%.3f") sprintf(paste0(fmt, " (", fmt, "\u2013", fmt, ")"), median(v), q(v, .1), q(v, .9))
pct <- function(x, d = 0) sprintf(paste0("%.", d, "f%%"), 100 * x)

## ---- 1. EIR grid --------------------------------------------------------------
say("## Equilibrium vs EIR (final 3 years)\n")
e <- eq %>% filter(grepl("^eir_", scenario)) %>% mutate(EIR = as.numeric(sub("eir_", "", scenario)))
ei <- e %>% filter(model == "IBM") %>% group_by(EIR) %>%
  summarise(`IBM realised EIR` = sprintf("%.1f", median(eir_realised)),
            `IBM PfPR (2-10)` = med_rng(pfpr_2_10),
            `IBM clinical (0-5, per child-year)` = med_rng(clin_0_5, "%.2f"),
            pf_i = median(pfpr_2_10), cl_i = median(clin_0_5), .groups = "drop")
eo <- e %>% filter(model == "blink") %>% transmute(EIR, `blink realised EIR` = sprintf("%.1f", eir_realised),
  `blink PfPR (2-10)` = sprintf("%.3f", pfpr_2_10), `blink clinical` = sprintf("%.2f", clin_0_5),
  pf_o = pfpr_2_10, cl_o = clin_0_5)
et <- left_join(ei, eo, by = "EIR") %>% arrange(EIR)
md_table(et %>% transmute(`init EIR` = EIR, `IBM realised EIR`, `blink realised EIR`,
                          `IBM PfPR (2-10)`, `blink PfPR (2-10)`,
                          `IBM clinical (0-5, per child-year)`, `blink clinical`))
say("max |blink - IBM median| PfPR(2-10): %.3f (at EIR %s); clinical incidence relative difference range: %s to %s\n",
    max(abs(et$pf_o - et$pf_i)), et$EIR[which.max(abs(et$pf_o - et$pf_i))],
    pct(min(et$cl_o / et$cl_i - 1), 1), pct(max(et$cl_o / et$cl_i - 1), 1))
in_band <- e %>% filter(model == "IBM") %>% group_by(EIR) %>%
  summarise(lo = q(pfpr_2_10, .1), hi = q(pfpr_2_10, .9), lo_c = q(clin_0_5, .1), hi_c = q(clin_0_5, .9), .groups = "drop") %>%
  left_join(eo, by = "EIR") %>% mutate(in_p = pf_o >= lo & pf_o <= hi, in_c = cl_o >= lo_c & cl_o <= hi_c)
say("blink inside the IBM 10-90%% band: PfPR at %d of %d EIRs; clinical at %d of %d\n",
    sum(in_band$in_p), nrow(in_band), sum(in_band$in_c), nrow(in_band))

## ---- 2. age profile -----------------------------------------------------------
say("## Age profile at EIR %s\n", EIR_REF)
a <- age %>% filter(scenario == paste0("eir_", EIR_REF))
ai <- a %>% filter(model == "IBM") %>% group_by(age_lo, age_hi, age_mid) %>%
  summarise(prev_i = median(prev), clin_i = median(clin), sev_i = median(sev),
            prev_lo = q(prev, .1), prev_hi = q(prev, .9), sev_lo = q(sev, .1), sev_hi = q(sev, .9), .groups = "drop")
ao <- a %>% filter(model == "blink") %>% select(age_mid, prev_o = prev, clin_o = clin, sev_o = sev)
at <- left_join(ai, ao, by = "age_mid") %>% arrange(age_lo)
md_table(at %>% transmute(`age band (y)` = sprintf("%g\u2013%g", age_lo, age_hi),
                          `IBM prevalence` = sprintf("%.3f", prev_i), `blink prevalence` = sprintf("%.3f", prev_o),
                          `IBM clinical` = sprintf("%.2f", clin_i), `blink clinical` = sprintf("%.2f", clin_o),
                          `IBM severe /1000` = sprintf("%.1f", sev_i), `blink severe /1000` = sprintf("%.1f", sev_o)))
yk <- at$age_hi <= 20                      # relative differences only where the rates are not tiny
say("max |blink - IBM| prevalence across bands: %.3f; under-20 bands: clinical relative diff range %s to %s, severe relative diff range %s to %s; blink severe inside IBM band in %d of %d bands\n",
    max(abs(at$prev_o - at$prev_i)), pct(min(at$clin_o[yk] / at$clin_i[yk] - 1), 1), pct(max(at$clin_o[yk] / at$clin_i[yk] - 1), 1),
    pct(min(at$sev_o[yk] / at$sev_i[yk] - 1), 0), pct(max(at$sev_o[yk] / at$sev_i[yk] - 1), 0),
    sum(at$sev_o >= at$sev_lo & at$sev_o <= at$sev_hi), nrow(at))

## ---- 2b. custom demography ---------------------------------------------------
if ("demography" %in% age$scenario) {
  say("## Custom demography at EIR %s\n", EIR_REF)
  dm <- age %>% filter(scenario == "demography")
  di <- dm %>% filter(model == "IBM") %>% group_by(age_lo, age_hi, age_mid) %>%
    summarise(pf_i = median(pop_frac), pr_i = median(prev), .groups = "drop")
  do <- dm %>% filter(model == "blink") %>% select(age_mid, pf_o = pop_frac, pr_o = prev)
  dt <- left_join(di, do, by = "age_mid") %>% arrange(age_lo)
  md_table(dt %>% transmute(`age band (y)` = sprintf("%g-%g", age_lo, age_hi),
                            `IBM population share` = pct(pf_i, 1), `blink population share` = pct(pf_o, 1),
                            `IBM prevalence` = sprintf("%.3f", pr_i), `blink prevalence` = sprintf("%.3f", pr_o)))
  u5 <- dt$age_hi <= 5
  say("under-5 share: IBM %s, blink %s; max |share diff| %.2f pp; max |prevalence diff| %.3f\n",
      pct(sum(dt$pf_i[u5]), 1), pct(sum(dt$pf_o[u5]), 1), 100 * max(abs(dt$pf_o - dt$pf_i)),
      max(abs(dt$pr_o - dt$pr_i)))
}

## ---- 3. seasonal cycle --------------------------------------------------------
say("## Seasonal cycle (final year, weekly bins)\n")
s <- doy %>% filter(scenario == "seasonal")
si <- s %>% filter(model == "IBM") %>% group_by(doy) %>% summarise(p = median(pfpr_2_10), c = median(clin_0_5), .groups = "drop")
so <- s %>% filter(model == "blink") %>% select(doy, p = pfpr_2_10, c = clin_0_5)
say("prevalence peak: IBM %.3f (day %d), blink %.3f (day %d); trough: IBM %.3f (day %d), blink %.3f (day %d)",
    max(si$p), round(si$doy[which.max(si$p)]), max(so$p), round(so$doy[which.max(so$p)]),
    min(si$p), round(si$doy[which.min(si$p)]), min(so$p), round(so$doy[which.min(so$p)]))
say("clinical peak (per child-year): IBM %.2f (day %d), blink %.2f (day %d); annual mean clinical IBM %.2f blink %.2f\n",
    max(si$c), round(si$doy[which.max(si$c)]), max(so$c), round(so$doy[which.max(so$c)]), mean(si$c), mean(so$c))
sea <- eq %>% filter(scenario == "seasonal")
say("seasonal realised EIR: IBM %.1f, blink %.1f (target %s); annual PfPR IBM %s, blink %.3f\n",
    median(sea$eir_realised[sea$model == "IBM"]), sea$eir_realised[sea$model == "blink"], EIR_REF,
    med_rng(sea$pfpr_2_10[sea$model == "IBM"]), sea$pfpr_2_10[sea$model == "blink"])

## ---- 4. country site files ----------------------------------------------------
vdir <- "C:/Users/pwinskil/Documents/dev/blink2/blink2_validate"
fs <- list.files(file.path(vdir, "results"), pattern = "_compare.rds$", full.names = TRUE)
if (length(fs)) {
  v <- bind_rows(lapply(fs, readRDS))
  ag <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    sprintf("n = %s, r = %.3f, slope = %.3f, relative bias = %s", format(length(x), big.mark = ","),
            cor(x, y), unname(coef(lm(y ~ x))[2]), pct(mean(y - x) / mean(x), 1)) }
  say("## Country site files\n")
  say("countries: %d; sub-sites: %d; years %d\u2013%d", length(unique(v$iso3c)),
      nrow(distinct(v, iso3c, name_1, urban_rural)), min(v$year), max(v$year))
  say("clinical: %s", ag(v$ms_clinical, v$mo_clinical))
  say("severe:   %s\n", ag(v$ms_severe, v$mo_severe))
}

## ---- 5. intervention impact ---------------------------------------------------
say("## Intervention impact: reduction over post-deployment years 0-3 vs pre-deployment years -3-0\n")
INT <- names(INT_LABELS)
red <- monthly %>% filter(scenario %in% INT) %>%
  mutate(phase = case_when(year >= BURN_Y - 3 & year < BURN_Y ~ "pre",
                           year >= BURN_Y & year < BURN_Y + 3 ~ "post", TRUE ~ NA_character_)) %>%
  filter(!is.na(phase)) %>% group_by(scenario, model, rep, phase) %>%
  summarise(pfpr = mean(pfpr_2_10), clin = mean(clin_0_5), .groups = "drop") %>%
  pivot_wider(names_from = phase, values_from = c(pfpr, clin)) %>%
  mutate(r_pfpr = 1 - pfpr_post / pfpr_pre, r_clin = 1 - clin_post / clin_pre)
ri <- red %>% filter(model == "IBM") %>% group_by(scenario) %>%
  summarise(`IBM prevalence reduction` = sprintf("%s (%s\u2013%s)", pct(median(r_pfpr)), pct(q(r_pfpr, .1)), pct(q(r_pfpr, .9))),
            `IBM clinical reduction` = sprintf("%s (%s\u2013%s)", pct(median(r_clin)), pct(q(r_clin, .1)), pct(q(r_clin, .9))),
            mp = median(r_pfpr), mc = median(r_clin), lp = q(r_pfpr, .1), hp = q(r_pfpr, .9), lc = q(r_clin, .1), hc = q(r_clin, .9),
            .groups = "drop")
ro <- red %>% filter(model == "blink") %>% transmute(scenario, `blink prevalence reduction` = pct(r_pfpr),
                                                   `blink clinical reduction` = pct(r_clin), op = r_pfpr, oc = r_clin)
rt <- left_join(ri, ro, by = "scenario") %>% mutate(scenario = factor(scenario, levels = INT)) %>% arrange(scenario)
md_table(rt %>% transmute(Scenario = sub("\n.*", "", INT_LABELS[as.character(scenario)]),
                          `IBM prevalence reduction`, `blink prevalence reduction`,
                          `IBM clinical reduction`, `blink clinical reduction`))
say("largest |blink - IBM median| gap: prevalence %.1f pp (%s), clinical %.1f pp (%s); blink inside the IBM band: prevalence %d/%d, clinical %d/%d\n",
    100 * max(abs(rt$op - rt$mp)), rt$scenario[which.max(abs(rt$op - rt$mp))],
    100 * max(abs(rt$oc - rt$mc)), rt$scenario[which.max(abs(rt$oc - rt$mc))],
    sum(rt$op >= rt$lp & rt$op <= rt$hp), nrow(rt), sum(rt$oc >= rt$lc & rt$oc <= rt$hc), nrow(rt))
## per-scenario post-deployment trajectory gap, years 0-6, as % of the IBM median (monthly, prevalence)
gap <- monthly %>% filter(scenario %in% INT, year >= BURN_Y, year < BURN_Y + 6) %>%
  group_by(scenario, model, year) %>% summarise(p = median(pfpr_2_10), c = median(clin_0_5), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = c(p, c)) %>% group_by(scenario) %>%
  summarise(`mean prevalence gap (blink - IBM, pp)` = sprintf("%+.1f", 100 * mean(p_blink - p_IBM)),
            `mean clinical gap (% of IBM)` = pct(mean(c_blink - c_IBM) / mean(c_IBM), 1), .groups = "drop")
md_table(gap)

## ---- 6. timing ----------------------------------------------------------------
say("## Run time\n")
tm <- timing %>% group_by(model) %>%
  summarise(runs = n(), `s per simulated year` = sprintf("%.2f", sum(elapsed_s) / sum(years)),
            `mean run (s)` = sprintf("%.1f", mean(elapsed_s)), `mean horizon (y)` = sprintf("%.0f", mean(years)), .groups = "drop")
md_table(tm)
say("IBM total CPU: %.1f h across %d runs; blink total: %.0f s across %d runs (IBM/blink per-year ratio %.0fx)",
    sum(timing$elapsed_s[timing$model == "IBM"]) / 3600, sum(timing$model == "IBM"),
    sum(timing$elapsed_s[timing$model == "blink"]), sum(timing$model == "blink"),
    (sum(timing$elapsed_s[timing$model == "IBM"]) / sum(timing$years[timing$model == "IBM"])) /
      (sum(timing$elapsed_s[timing$model == "blink"]) / sum(timing$years[timing$model == "blink"])))

writeLines(out, file.path(DDIR, "tables.md"))
cat(out, sep = "\n")
