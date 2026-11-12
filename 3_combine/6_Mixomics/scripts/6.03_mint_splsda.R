# ==============================================================================
# 6.03_mint_splsda.R
#
# Prediction of reef zoning status (NTMR vs fished) from community composition
# with MINT sPLS-DA, GBR sector as the study variable (Fig. S3b).
#   1. Setup and input
#   2. Tuning: number of features per component (LOGOCV)
#   3. Final two-component model and cross-validated performance
#   4. Sample plot (Fig. S3b)
#   5. Indicator features: direction of enrichment and consistency across sectors
#
# Tuning and evaluation use leave-one-group-out cross-validation (LOGOCV, one
# GBR sector left out at a time), balanced error rate (BER), and centroids
# distance. Balanced accuracy = 1 - BER. Direction of each selected feature is
# taken from the difference in mean CLR between no-take and fished reefs, as
# the sign of (s)PLS loadings is arbitrary.
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.03_mint_splsda.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01)
# Output: tuning and final model (rds), performance and indicators (xlsx),
#         sample plot (svg)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)   # v6.26.0
  library(BiocParallel)
  library(openxlsx)
  library(tidyverse)
})
select <- dplyr::select

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
profile_dir <- file.path(in_dir, "0_profiles")   # output of 6.01
out_dir <- file.path(in_dir, "2_MINT_sPLSDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col    <- "FishingStatus"
zone_closed <- "Closed"   # NTMR
zone_open   <- "Open"     # fished
sector_col  <- "GBRSector"
fish_cols   <- c(Closed = "#41ce80", Open = "#37648b")

ncomp_tune  <- 5
test_keepX  <- seq(10, 500, 10)
ncomp_final <- 2
dist_used   <- "centroids.dist"
param       <- MulticoreParam(workers = 8)

inp   <- readRDS(file.path(profile_dir, paste0("zoning_input_", tag, ".rds")))
X     <- inp$clr
Y     <- inp$meta[[zone_col]]
study <- factor(inp$meta[[sector_col]])
tax   <- as.data.frame(as(tax_table(inp$ps), "matrix"))[colnames(X), , drop = FALSE]

print(table(Y, study))


# ==============================================================================
# 2. TUNING
# ==============================================================================

tune_rds <- file.path(out_dir, paste0("tune_mint_splsda_", tag, ".rds"))
if (file.exists(tune_rds)) {
  tune_mint <- readRDS(tune_rds)
} else {
  tune_mint <- tune(X = X, Y = Y, study = study,
                    ncomp = ncomp_tune, test.keepX = test_keepX,
                    method = "mint.splsda", measure = "BER", dist = dist_used)
  saveRDS(tune_mint, tune_rds)
}

print(tune_mint$choice.ncomp)
print(tune_mint$choice.keepX)

svg(file.path(out_dir, paste0("tune_mint_splsda_", tag, ".svg")), width = 7.5, height = 7.5)
plot(tune_mint, sd = FALSE)
dev.off()


# ==============================================================================
# 3. FINAL MODEL AND PERFORMANCE
# ==============================================================================

keepX <- tune_mint$choice.keepX[seq_len(ncomp_final)]

final_mint <- mint.splsda(X = X, Y = Y, study = study,
                          ncomp = ncomp_final, keepX = keepX)
saveRDS(final_mint, file.path(out_dir, paste0("mint_splsda_final_", tag, ".rds")))

perf_mint <- perf(final_mint, dist = dist_used, auc = TRUE,
                  BPPARAM = param, progressBar = FALSE)
saveRDS(perf_mint, file.path(out_dir, paste0("perf_mint_splsda_final_", tag, ".rds")))

ber          <- perf_mint$global.error$BER[ncomp_final, dist_used]
balanced_acc <- 1 - ber
auc_cv       <- perf_mint$auc$comp1[[1]]["Closed vs Open", "AUC"]

cat(sprintf("%s: keepX = %s; balanced accuracy (LOGOCV) = %.1f%%; AUC (comp 1) = %.3f\n",
            tag, paste(keepX, collapse = ", "), 100 * balanced_acc, auc_cv))

# Per-sector error rates
sector_error <- imap_dfr(perf_mint$study.specific.error, function(e, s) {
  as.data.frame(e$error.rate.class[[dist_used]]) %>%
    rownames_to_column("Class") %>%
    mutate(Sector = s, .before = 1)
})

write.xlsx(list(summary = tibble(habitat = hab, profile = method,
                                 keepX = paste(keepX, collapse = ", "),
                                 BER = ber, balanced_accuracy = balanced_acc,
                                 AUC_comp1 = auc_cv),
                global_BER = as.data.frame(perf_mint$global.error$BER),
                sector_error = sector_error),
           file.path(out_dir, paste0("performance_mint_splsda_", tag, ".xlsx")))


# ==============================================================================
# 4. SAMPLE PLOT (Fig. S3b)
# ==============================================================================
# Global components; ellipses at 95% confidence; point shape by GBR sector

p_global <- plotIndiv(final_mint, study = "global", comp = c(1, 2),
                      legend = TRUE, ellipse = TRUE,
                      col = unname(fish_cols[levels(Y)]),
                      subtitle = paste("MINT sPLS-DA,", tag))
ggsave(file.path(out_dir, paste0("mint_splsda_", tag, "_global.svg")),
       p_global$graph, width = 7.5, height = 7.5)


# ==============================================================================
# 5. INDICATOR FEATURES
# ==============================================================================

direction_of <- function(d) if_else(d > 0, "NTMR-enriched", "Fished-enriched")

zone_clr_diff <- colMeans(X[Y == zone_closed, , drop = FALSE]) -
  colMeans(X[Y == zone_open, , drop = FALSE])

sel <- selectVar(final_mint, comp = 1)
indicators <- tibble(Feature = sel$name, loading = sel$value$value.var) %>%
  mutate(clr_diff = zone_clr_diff[Feature], direction = direction_of(clr_diff)) %>%
  left_join(rownames_to_column(tax, "Feature"), by = "Feature")

print(count(indicators, direction))

# Direction of each indicator within each sector
sector_diff <- sapply(levels(study), function(s) {
  i <- study == s
  colMeans(X[i & Y == zone_closed, indicators$Feature, drop = FALSE]) -
    colMeans(X[i & Y == zone_open,   indicators$Feature, drop = FALSE])
})

write.xlsx(list(indicators_comp1 = indicators,
                clr_diff_by_sector = rownames_to_column(as.data.frame(sector_diff), "Feature")),
           file.path(out_dir, paste0("indicators_mint_splsda_", tag, ".xlsx")))

p_sector <- sector_diff %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Sector", values_to = "clr_diff") %>%
  mutate(Direction = direction_of(clr_diff),
         Feature = factor(Feature, levels = indicators$Feature[order(indicators$clr_diff)])) %>%
  ggplot(aes(Sector, Feature, fill = Direction)) +
  geom_tile() +
  scale_fill_manual(values = c(`NTMR-enriched` = fish_cols[["Closed"]],
                               `Fished-enriched` = fish_cols[["Open"]])) +
  labs(y = "MINT sPLS-DA indicators (component 1)", title = paste("Indicator direction per sector,", tag)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(out_dir, paste0("indicators_by_sector_", tag, ".svg")),
       p_sector, width = 7.5, height = 7.5)
