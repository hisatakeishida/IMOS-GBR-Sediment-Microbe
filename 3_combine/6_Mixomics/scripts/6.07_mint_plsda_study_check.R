# ==============================================================================
# 6.07_mint_plsda_study_check.R
#
# Supporting analysis (not reported): MINT PLS-DA, and a check of the choice
# of study variable for MINT sPLS-DA (GBR sector versus sampling trip).
#   1. Setup and input
#   2. MINT PLS-DA (all features): number of components, final model, AUC
#   3. MINT sPLS-DA error rates with sector or trip as the study variable
#   4. Overlap of selected features with sector or trip as the study variable
#   5. Per-sector error rates of the 10-component MINT sPLS-DA
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.07_mint_plsda_study_check.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01)
# Output: figures (png) and tables (xlsx)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(ggVennDiagram)
  library(openxlsx)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
out_dir <- file.path(in_dir, "7_MINT_PLSDA_study_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col   <- "FishingStatus"
sector_col <- "GBRSector"
trip_col   <- "TripDescription"

inp    <- readRDS(file.path(in_dir, paste0("zoning_input_", tag, ".rds")))
X      <- inp$clr
Y      <- inp$meta[[zone_col]]
sector <- factor(inp$meta[[sector_col]])
trip   <- factor(inp$meta[[trip_col]])

save_png <- function(file, expr, width = 7.5, height = 7.5) {
  png(file.path(out_dir, file), width = width, height = height, units = "in", res = 300)
  on.exit(dev.off())
  force(expr)
}


# ==============================================================================
# 2. MINT PLS-DA
# ==============================================================================

mint_pls10      <- mint.plsda(X = X, Y = Y, study = trip, ncomp = 10)
perf_mint_pls10 <- perf(mint_pls10, auc = TRUE, progressBar = FALSE)
save_png(paste0("MINT_PLSDA_", tag, "_error_rate_trip.png"), plot(perf_mint_pls10))
print(perf_mint_pls10$global.error)

mint_pls      <- mint.plsda(X = X, Y = Y, study = sector, ncomp = 2)
perf_mint_pls <- perf(mint_pls, auc = TRUE, progressBar = FALSE)
auc_cv <- perf_mint_pls$auc$comp1[[1]]["Closed vs Open", "AUC"]
cat(sprintf("MINT PLS-DA (sector): BER comp 2 = %.3f; cross-validated AUC comp 1 = %.3f\n",
            perf_mint_pls$global.error$BER[2, "centroids.dist"], auc_cv))

roc_df <- perf_mint_pls$auc$comp1$df
p_roc <- ggplot(roc_df, aes(Specificity, Sensitivity)) +
  geom_line(colour = "red", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey") +
  labs(title = sprintf("MINT PLS-DA, cross-validated ROC (LOGOCV), AUC = %.3f", auc_cv)) +
  theme_bw()
ggsave(file.path(out_dir, paste0("MINT_PLSDA_", tag, "_LOGOCV_ROC.png")), p_roc, width = 7.5, height = 7.5)


# ==============================================================================
# 3. MINT sPLS-DA ERROR RATES: SECTOR OR TRIP AS STUDY
# ==============================================================================

mint_sector10 <- mint.splsda(X = X, Y = Y, study = sector, ncomp = 10)
mint_trip10   <- mint.splsda(X = X, Y = Y, study = trip,   ncomp = 10)
perf_sector10 <- perf(mint_sector10)
perf_trip10   <- perf(mint_trip10)
save_png(paste0("MINT_sPLSDA_", tag, "_error_rate_sector.png"), plot(perf_sector10))
save_png(paste0("MINT_sPLSDA_", tag, "_error_rate_trip.png"),   plot(perf_trip10))


# ==============================================================================
# 4. OVERLAP OF SELECTED FEATURES
# ==============================================================================
# 100 features on component 1 with sector or trip as the study variable

vars_sector <- selectVar(mint.splsda(X = X, Y = Y, study = sector, ncomp = 2, keepX = c(100, 10)), comp = 1)$name
vars_trip   <- selectVar(mint.splsda(X = X, Y = Y, study = trip,   ncomp = 2, keepX = c(100, 10)), comp = 1)$name
cat(sprintf("Shared features (component 1, keepX = 100): %d\n", length(intersect(vars_sector, vars_trip))))
ggsave(file.path(out_dir, paste0("MINT_sPLSDA_", tag, "_sector_vs_trip_venn.png")),
       ggVennDiagram(list(Sector = vars_sector, Trip = vars_trip)), width = 6, height = 5)


# ==============================================================================
# 5. PER-SECTOR ERROR RATES
# ==============================================================================

sector_overall <- imap_dfr(perf_sector10$study.specific.error, function(e, s) {
  as.data.frame(e$overall) %>% rownames_to_column("Component") %>% mutate(Sector = s, .before = 1)
})
sector_class <- imap_dfr(perf_sector10$study.specific.error, function(e, s) {
  as.data.frame(e$error.rate.class$centroids.dist) %>% rownames_to_column("Class") %>%
    mutate(Sector = s, .before = 1)
})

write.xlsx(list(MINT_PLSDA_global = as.data.frame(perf_mint_pls$global.error$BER),
                sector_overall = sector_overall, sector_by_class = sector_class,
                shared_features = tibble(Feature = intersect(vars_sector, vars_trip))),
           file.path(out_dir, paste0("MINT_study_check_", tag, ".xlsx")), rowNames = TRUE)
