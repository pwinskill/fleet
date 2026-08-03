# Regenerate the comparison PNGs from the saved CSVs (comparison/data/) using the
# shared builders in plot_helpers.R. No model runs — pure re-plot, so it is cheap
# to iterate on figure design without repeating the 30-year IBM burn-in sweep.
.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")
ROOT <- "C:/Users/pwinskil/Documents/dev/blink2/malariaode/comparison"
source(file.path(ROOT, "plot_helpers.R"))
DDIR <- file.path(ROOT, "data"); PDIR <- file.path(ROOT, "plots")
rd <- function(f) read.csv(file.path(DDIR, f))
sv <- function(f, g, w = 7, h = 5) ggplot2::ggsave(file.path(PDIR, f), g, width = w, height = h, dpi = 130)

sv("A_pfpr_eir.png",        plot_A(rd("A_pfpr_eir.csv")))
sv("B_age_profile.png",     plot_B(rd("B_age_profile.csv")))
sv("C_bednets.png",         plot_C(rd("C_bednets.csv")), w = 8)
sv("D_seasonal.png",        plot_D(rd("D_seasonal.csv")), w = 8)
sv("E_smc.png",             plot_E(rd("E_smc.csv")), w = 8)
sv("F_demog_structure.png", plot_F_structure(rd("F_demog_structure.csv")))
sv("F_demog_prevalence.png",plot_F_prev(rd("F_demog_prevalence.csv")))
# incidence figures (clinical + severe), from the shared INC_SPECS
for (spec in INC_SPECS)
  sv(paste0(spec$file, ".png"), render_incidence(spec, rd(paste0(spec$file, ".csv"))), w = 8, h = 6)
cat("re-rendered 10 PNGs from CSVs\n")
