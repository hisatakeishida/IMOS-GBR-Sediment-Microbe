# ==============================================================================
# 6.09_random_forest.R
#
# Supporting analysis (not reported): random forest classification of zoning
# status with leave-one-sector-out cross-validation, a permutation test, and
# overlap of important features with the MINT sPLS-DA indicators (6.03).
#   1. Setup and input
#   2. Leave-one-sector-out cross-validation
#   3. Permutation test
#   4. Feature importance and overlap with MINT sPLS-DA indicators
#
# Forests are grown on CLR values with class-balanced subsampling (equal
# numbers of NTMR and fished samples per tree). Accuracy is the mean over the
# held-out sectors of the proportion of samples predicted correctly. Null
# forests use the same settings as the observed forests, with zoning status
# shuffled among reefs within each sector (as in 6.04).
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.09_random_forest.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01); mint_splsda_final_<HAB>_<METHOD>.rds (6.03)
# Output: accuracies and importance (xlsx, rds), permutation histogram (png)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(randomForest)
  library(mixOmics)
  library(parallel)
  library(openxlsx)
  library(tidyverse)
})
select  <- dplyr::select

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
out_dir <- file.path(in_dir, "9_random_forest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col   <- "FishingStatus"
sector_col <- "GBRSector"
reef_col   <- "ReefName"
ntree      <- 1000
n_perm     <- 999
n_top      <- 350      # top features by importance compared with MINT indicators
n_cores    <- max(1, detectCores() - 4)

inp     <- readRDS(file.path(in_dir, paste0("zoning_input_", tag, ".rds")))
X       <- as.data.frame(inp$clr)
Y       <- inp$meta[[zone_col]]
sector  <- as.character(inp$meta[[sector_col]])
reef    <- inp$meta[[reef_col]]
sectors <- unique(sector)
tax     <- as.data.frame(as(tax_table(inp$ps), "matrix"))

reef_tbl <- inp$meta[!duplicated(reef), c(reef_col, zone_col, sector_col)]
shuffle_zone_by_reef <- function() {
  z <- ave(as.character(reef_tbl[[zone_col]]), reef_tbl[[sector_col]], FUN = sample)
  factor(z[match(reef, reef_tbl[[reef_col]])], levels = levels(Y))
}

grow_forest <- function(x, y, importance = FALSE) {
  randomForest(x = x, y = y, ntree = ntree, mtry = floor(sqrt(ncol(x))),
               importance = importance, strata = y, sampsize = rep(min(table(y)), nlevels(y)))
}

leave_one_sector_out <- function(y, importance = FALSE) {
  lapply(sectors, function(s) {
    test <- sector == s
    rf   <- grow_forest(X[!test, , drop = FALSE], y[!test], importance)
    list(sector = s, accuracy = mean(predict(rf, X[test, , drop = FALSE]) == y[test]),
         importance = if (importance) importance(rf) else NULL)
  })
}


# ==============================================================================
# 2. LEAVE-ONE-SECTOR-OUT CROSS-VALIDATION
# ==============================================================================

set.seed(123)
cv <- leave_one_sector_out(Y, importance = TRUE)
cv_accuracy <- tibble(Sector = map_chr(cv, "sector"), Accuracy = map_dbl(cv, "accuracy"))
observed_accuracy <- mean(cv_accuracy$Accuracy)
print(cv_accuracy)
cat(sprintf("Random forest mean accuracy: %.1f%% (SD %.1f%%)\n",
            100 * observed_accuracy, 100 * sd(cv_accuracy$Accuracy)))


# ==============================================================================
# 3. PERMUTATION TEST
# ==============================================================================

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
null_accuracy <- unlist(mclapply(seq_len(n_perm), function(i) {
  mean(map_dbl(leave_one_sector_out(shuffle_zone_by_reef()), "accuracy"))
}, mc.cores = n_cores, mc.set.seed = TRUE))

p_value <- (sum(null_accuracy >= observed_accuracy) + 1) / (n_perm + 1)
null_ci <- quantile(null_accuracy, c(0.025, 0.975))
cat(sprintf("Null mean %.1f%%; p = %.4f\n", 100 * mean(null_accuracy), p_value))

p_hist <- ggplot(tibble(accuracy = null_accuracy), aes(accuracy)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", alpha = 0.7) +
  geom_vline(xintercept = observed_accuracy, colour = "red3", linewidth = 1.5) +
  geom_vline(xintercept = null_ci, colour = "blue", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Accuracy (leave one sector out)", y = "Frequency",
       title = paste("Random forest permutation test,", tag),
       subtitle = sprintf("Observed %.1f%%, p = %.3f (%d permutations)",
                          100 * observed_accuracy, p_value, n_perm)) +
  theme_minimal()
ggsave(file.path(out_dir, paste0("RF_permutation_", tag, ".png")), p_hist, width = 7.5, height = 5)


# ==============================================================================
# 4. FEATURE IMPORTANCE AND OVERLAP WITH MINT sPLS-DA INDICATORS
# ==============================================================================

avg_importance <- Reduce(`+`, lapply(cv, `[[`, "importance")) / length(cv)
importance_tbl <- tibble(Feature = rownames(avg_importance),
                         MeanDecreaseAccuracy = avg_importance[, "MeanDecreaseAccuracy"],
                         MeanDecreaseGini     = avg_importance[, "MeanDecreaseGini"]) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  left_join(rownames_to_column(tax, "Feature"), by = "Feature")

final_mint      <- readRDS(file.path(in_dir, "2_MINT_sPLSDA", paste0("mint_splsda_final_", tag, ".rds")))
mint_indicators <- selectVar(final_mint, comp = 1)$name
top_rf          <- head(importance_tbl$Feature, n_top)
overlap         <- intersect(top_rf, mint_indicators)
cat(sprintf("MINT indicators (comp 1): %d; in RF top %d: %d (%.1f%% of indicators)\n",
            length(mint_indicators), n_top, length(overlap),
            100 * length(overlap) / length(mint_indicators)))

saveRDS(list(cv_accuracy = cv_accuracy, null_accuracy = null_accuracy, p_value = p_value,
             importance = importance_tbl),
        file.path(out_dir, paste0("RF_", tag, ".rds")))
write.xlsx(list(summary = tibble(observed_accuracy, null_mean = mean(null_accuracy),
                                 null_sd = sd(null_accuracy), null_q2.5 = null_ci[[1]],
                                 null_q97.5 = null_ci[[2]], p_value,
                                 mint_indicators = length(mint_indicators),
                                 overlap_top = length(overlap)),
                cv_accuracy = cv_accuracy, importance = importance_tbl,
                overlap = tibble(Feature = overlap)),
           file.path(out_dir, paste0("RF_", tag, ".xlsx")))
