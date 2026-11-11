# ==============================================================================
# 6.04_mint_permutation.R
#
# Permutation test of zoning status classification by MINT sPLS-DA
# (Fig. S3c for MAG-based and Fig. S3d for mOTU-based profiles).
#   1. Setup and input
#   2. Observed balanced accuracy (final model from 6.03)
#   3. Null distribution: zoning status permuted among reefs within each sector
#   4. p-value, summary, and histogram
#
# Zoning status is a reef-level attribute, so labels are shuffled among reefs
# within each GBR sector (all samples of a reef keep the same label). Each
# permutation refits MINT sPLS-DA with the tuned parameters (ncomp, keepX) and
# evaluates it by LOGOCV (leave one sector out), as for the observed model.
# p = (number of null accuracies >= observed + 1) / (number of permutations + 1)
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.04_mint_permutation.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01); final model and performance (6.03)
# Output: null accuracies (xlsx, rds) and histogram (svg)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(openxlsx)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir   <- "<path/to/6_Mixomics>"
mint_dir <- file.path(in_dir, "2_MINT_sPLSDA")
out_dir  <- file.path(in_dir, "3_permutation_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col   <- "FishingStatus"
sector_col <- "GBRSector"
reef_col   <- "ReefName"
dist_used  <- "centroids.dist"
n_perm     <- 999

inp   <- readRDS(file.path(in_dir, paste0("zoning_input_", tag, ".rds")))
X     <- inp$clr
Y     <- inp$meta[[zone_col]]
study <- factor(inp$meta[[sector_col]])
reef  <- inp$meta[[reef_col]]

final_mint <- readRDS(file.path(mint_dir, paste0("mint_splsda_final_", tag, ".rds")))
perf_mint  <- readRDS(file.path(mint_dir, paste0("perf_mint_splsda_final_", tag, ".rds")))
ncomp <- final_mint$ncomp
keepX <- final_mint$keepX


# ==============================================================================
# 2. OBSERVED BALANCED ACCURACY
# ==============================================================================

observed_accuracy <- 1 - perf_mint$global.error$BER[ncomp, dist_used]
cat(sprintf("%s observed balanced accuracy: %.1f%%\n", tag, 100 * observed_accuracy))


# ==============================================================================
# 3. NULL DISTRIBUTION
# ==============================================================================

reef_tbl <- inp$meta[!duplicated(reef), c(reef_col, zone_col, sector_col)]

shuffle_zone_by_reef <- function() {
  z <- ave(as.character(reef_tbl[[zone_col]]), reef_tbl[[sector_col]], FUN = sample)
  factor(z[match(reef, reef_tbl[[reef_col]])], levels = levels(Y))
}

set.seed(123)
null_accuracy <- rep(NA_real_, n_perm)
start_time    <- Sys.time()

for (i in seq_len(n_perm)) {
  if (i == 1 || i %% 50 == 0) {
    cat(sprintf("  permutation %d/%d (%.1f min)\n", i, n_perm,
                as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
  }
  Y_perm <- shuffle_zone_by_reef()
  null_accuracy[i] <- tryCatch({
    fit <- mint.splsda(X = X, Y = Y_perm, study = study, ncomp = ncomp, keepX = keepX)
    1 - perf(fit, dist = dist_used, progressBar = FALSE)$global.error$BER[ncomp, dist_used]
  }, error = function(e) {
    message("permutation ", i, " failed: ", conditionMessage(e))
    NA_real_
  })
}

n_valid       <- sum(!is.na(null_accuracy))
null_accuracy <- null_accuracy[!is.na(null_accuracy)]
cat(sprintf("Successful permutations: %d/%d\n", n_valid, n_perm))


# ==============================================================================
# 4. p-VALUE, SUMMARY, AND HISTOGRAM
# ==============================================================================

p_value <- (sum(null_accuracy >= observed_accuracy) + 1) / (length(null_accuracy) + 1)
null_ci <- quantile(null_accuracy, c(0.025, 0.975))

summary_tbl <- tibble(
  habitat = hab, profile = method,
  observed_accuracy = observed_accuracy,
  null_mean = mean(null_accuracy), null_sd = sd(null_accuracy),
  null_q2.5 = null_ci[[1]], null_q97.5 = null_ci[[2]],
  n_permutations = length(null_accuracy), p_value = p_value
)
print(summary_tbl)

saveRDS(list(summary = summary_tbl, null_accuracy = null_accuracy),
        file.path(out_dir, paste0("permutation_mint_splsda_", tag, ".rds")))
write.xlsx(list(summary = summary_tbl,
                null_accuracy = tibble(Permutation = seq_along(null_accuracy),
                                       Accuracy = null_accuracy)),
           file.path(out_dir, paste0("permutation_mint_splsda_", tag, ".xlsx")))

# Red solid line: observed accuracy; blue dashed lines: 2.5th and 97.5th
# percentiles of the null distribution
p_hist <- ggplot(tibble(accuracy = null_accuracy), aes(accuracy)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", alpha = 0.7) +
  geom_vline(xintercept = observed_accuracy, colour = "red3", linewidth = 1.5) +
  geom_vline(xintercept = null_ci, colour = "blue", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Balanced classification accuracy (LOGOCV)", y = "Frequency",
       title = paste("MINT sPLS-DA permutation test,", tag),
       subtitle = sprintf("Observed %.1f%%, p = %.3f (%d permutations)",
                          100 * observed_accuracy, p_value, length(null_accuracy))) +
  theme_minimal()
ggsave(file.path(out_dir, paste0("permutation_mint_splsda_", tag, ".svg")),
       p_hist, width = 7.5, height = 5)
