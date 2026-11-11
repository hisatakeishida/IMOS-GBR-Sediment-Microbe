# ==============================================================================
# 6.08_loro_validation.R
#
# Supporting analysis (not reported): leave-one-reef-out (LORO) validation of
# zoning status prediction, with the MINT sPLS-DA parameters (6.03).
#   1. Setup and input
#   2. LORO with sPLS-DA
#   3. Summary
#   4. Figures
#
# MINT models need at least one sample of every study in training and cannot
# predict a reef held out on its own, so LORO uses single-study sPLS-DA with
# the ncomp and keepX of the final MINT sPLS-DA model. Each held-out reef
# carries a single zoning status, so per-reef accuracy is the proportion of its
# samples predicted correctly; the overall summary also gives balanced accuracy
# over all held-out samples.
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.08_loro_validation.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01); mint_splsda_final_<HAB>_<METHOD>.rds (6.03)
# Output: per-reef results (xlsx, rds) and figure (png)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(parallel)
  library(patchwork)
  library(openxlsx)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
out_dir <- file.path(in_dir, "8_LORO")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col <- "FishingStatus"
reef_col <- "ReefName"
n_cores  <- max(1, detectCores() - 2)

inp  <- readRDS(file.path(in_dir, paste0("zoning_input_", tag, ".rds")))
X    <- inp$clr
Y    <- inp$meta[[zone_col]]
reef <- inp$meta[[reef_col]]

final_mint <- readRDS(file.path(in_dir, "2_MINT_sPLSDA", paste0("mint_splsda_final_", tag, ".rds")))
ncomp <- final_mint$ncomp
keepX <- final_mint$keepX[seq_len(ncomp)]


# ==============================================================================
# 2. LORO WITH sPLS-DA
# ==============================================================================

set.seed(123)
loro <- mclapply(unique(reef), function(r) {
  test <- reef == r
  out  <- tryCatch({
    fit  <- splsda(X[!test, , drop = FALSE], Y[!test], ncomp = ncomp, keepX = keepX, scale = TRUE)
    pred <- predict(fit, newdata = X[test, , drop = FALSE], dist = "centroids.dist")$class$centroids.dist[, ncomp]
    tibble(Sample = rownames(X)[test], Pred = pred, Status = "Success")
  }, error = function(e) tibble(Sample = rownames(X)[test], Pred = NA_character_,
                                Status = conditionMessage(e)))
  mutate(out, Reef = r, Truth = as.character(Y[test]), .before = 1)
}, mc.cores = n_cores) %>%
  bind_rows()

print(count(distinct(loro, Reef, Status), Status))


# ==============================================================================
# 3. SUMMARY
# ==============================================================================

ok <- filter(loro, Status == "Success")

per_reef <- ok %>%
  group_by(Reef, Truth) %>%
  summarise(N_Samples = n(), Correct = sum(Pred == Truth), Accuracy = Correct / N_Samples, .groups = "drop")

overall <- tibble(
  reefs                      = n_distinct(per_reef$Reef),
  mean_reef_accuracy         = mean(per_reef$Accuracy),
  sd_reef_accuracy           = sd(per_reef$Accuracy),
  median_reef_accuracy       = median(per_reef$Accuracy),
  sample_balanced_accuracy   = mean(tapply(ok$Pred == ok$Truth, ok$Truth, mean))
)
print(overall)

confusion <- count(ok, Truth, Pred)
print(confusion)

saveRDS(list(per_sample = loro, per_reef = per_reef, overall = overall),
        file.path(out_dir, paste0("LORO_splsda_", tag, ".rds")))
write.xlsx(list(overall = overall, per_reef = per_reef, confusion = confusion, per_sample = loro),
           file.path(out_dir, paste0("LORO_splsda_", tag, ".xlsx")))


# ==============================================================================
# 4. FIGURES
# ==============================================================================

p_hist <- ggplot(per_reef, aes(100 * Accuracy)) +
  geom_histogram(bins = 15, fill = "steelblue", colour = "black", alpha = 0.7) +
  geom_vline(xintercept = 100 * overall$mean_reef_accuracy, linetype = "dashed", colour = "red3") +
  geom_vline(xintercept = 50, linetype = "dotted", colour = "grey50") +
  labs(x = "Accuracy per held-out reef (%)", y = "Reefs",
       title = paste("Leave-one-reef-out sPLS-DA,", tag)) +
  theme_minimal()

p_conf <- ggplot(confusion, aes(Truth, Pred, fill = n)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = n), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue4", name = "Samples") +
  labs(x = "Actual", y = "Predicted",
       subtitle = sprintf("Balanced accuracy over samples: %.1f%%", 100 * overall$sample_balanced_accuracy)) +
  theme_minimal()

ggsave(file.path(out_dir, paste0("LORO_splsda_", tag, ".png")), p_hist / p_conf, width = 7.5, height = 10)
