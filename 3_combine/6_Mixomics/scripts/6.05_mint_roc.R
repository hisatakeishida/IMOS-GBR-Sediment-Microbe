# ==============================================================================
# 6.05_mint_roc.R
#
# Cross-validated ROC for MINT sPLS-DA zoning classification (LOGOCV).
# Each GBR sector is held out in turn; the model is refitted on the remaining
# sectors with the tuned ncomp and keepX from 6.03, and the held-out sector is
# predicted as a new study. Score = predicted value for "Closed" on comp 1.
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.05_mint_roc.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01); final model (6.03)
# Output: held-out scores (xlsx) and ROC curve (svg)
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(pROC)
  library(openxlsx)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir      <- "<path/to/6_Mixomics>"
profile_dir <- file.path(in_dir, "0_profiles")
mint_dir    <- file.path(in_dir, "2_MINT_sPLSDA")

zone_col    <- "FishingStatus"
zone_closed <- "Closed"
zone_open   <- "Open"
sector_col  <- "GBRSector"
dist_used   <- "centroids.dist"
fish_cols   <- c(Closed = "#41ce80", Open = "#37648b")

inp   <- readRDS(file.path(profile_dir, paste0("zoning_input_", tag, ".rds")))
X     <- inp$clr
Y     <- factor(inp$meta[[zone_col]])
study <- factor(inp$meta[[sector_col]])

final_mint <- readRDS(file.path(mint_dir, paste0("mint_splsda_final_", tag, ".rds")))
ncomp <- final_mint$ncomp
keepX <- final_mint$keepX

# Held-out predictions, one sector at a time
cv_scores <- map_dfr(levels(study), function(s) {
  test <- study == s
  fit  <- mint.splsda(X = X[!test, , drop = FALSE], Y = droplevels(Y[!test]),
                      study = droplevels(study[!test]), ncomp = ncomp, keepX = keepX)
  pred <- predict(fit, newdata = X[test, , drop = FALSE],
                  study.test = factor(study[test]), dist = dist_used)
  tibble(sample = rownames(X)[test], Sector = s, truth = Y[test],
         score  = pred$predict[, zone_closed, 1])
})

roc_cv     <- roc(cv_scores$truth, cv_scores$score,
                  levels = c(zone_open, zone_closed), direction = "<", quiet = TRUE)
auc_pooled <- as.numeric(auc(roc_cv))
cat(sprintf("%s: pooled LOGOCV AUC (comp 1) = %.3f\n", tag, auc_pooled))

roc_df <- tibble(fpr = 1 - roc_cv$specificities, tpr = roc_cv$sensitivities) %>%
  arrange(fpr, tpr)

p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_path(colour = fish_cols[["Closed"]], linewidth = 1.2) +
  coord_equal() +
  labs(title = sprintf("MINT sPLS-DA, %s: cross-validated ROC (LOGOCV)", tag),
       subtitle = sprintf("AUC = %.3f", auc_pooled),
       x = "1 - Specificity (false positive rate)",
       y = "Sensitivity (true positive rate)") +
  theme_minimal()

ggsave(file.path(mint_dir, paste0("mint_splsda_", tag, "_LOGOCV_roc.svg")),
       p_roc, width = 7.5, height = 7.5, bg = "white")
write.xlsx(cv_scores, file.path(mint_dir, paste0("cv_scores_mint_splsda_", tag, ".xlsx")))