# ==============================================================================
# 6.06_plsda_splsda.R
#
# Supporting analysis (not reported): single-study PLS-DA and sPLS-DA of
# zoning status, with sample-level and leave-one-reef-out (LORO) validation.
#   1. Setup and input
#   2. PLS-DA: number of components and final two-component model
#   3. sPLS-DA: tuning and final model
#   4. Performance, AUC, and feature stability
#   5. Sample plots and clustered image map (CIM)
#   6. Leave-one-reef-out validation of PLS-DA and sPLS-DA
#
# Sample-level M-fold cross-validation places replicate samples of a reef in
# both training and test sets, so its error rates are optimistic for a
# reef-level attribute; section 6 holds out whole reefs instead.
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.06_plsda_splsda.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01)
# Output: models and performance (rds), figures (png), LORO results (xlsx)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(BiocParallel)
  library(parallel)
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
out_dir <- file.path(in_dir, "6_PLSDA_sPLSDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col <- "FishingStatus"
reef_col <- "ReefName"
trip_col <- "TripDescription"
fish_cols <- c(Closed = "#41ce80", Open = "#37648b")
param   <- MulticoreParam(workers = 8)
n_cores <- max(1, detectCores() - 2)

inp  <- readRDS(file.path(profile_dir, paste0("zoning_input_", tag, ".rds")))
X    <- inp$clr
Y    <- inp$meta[[zone_col]]
reef <- inp$meta[[reef_col]]
trip <- factor(inp$meta[[trip_col]])

balanced_accuracy <- function(truth, pred) {
  truth <- as.character(truth)
  mean(tapply(as.character(pred) == truth, truth, mean))
}

save_png <- function(file, expr, width = 7.5, height = 7.5) {
  png(file.path(out_dir, file), width = width, height = height, units = "in", res = 300)
  on.exit(dev.off())
  force(expr)
}


# ==============================================================================
# 2. PLS-DA
# ==============================================================================

set.seed(123)
plsda_10   <- plsda(X, Y, ncomp = 10)
perf_pls10 <- perf(plsda_10, validation = "Mfold", folds = 4, nrepeat = 50,
                   BPPARAM = param, progressBar = FALSE)
save_png(paste0("PLSDA_", tag, "_error_rate.png"), plot(perf_pls10, sd = TRUE, legend.position = "horizontal"))
print(perf_pls10$error.rate$BER)

ncomp_pls  <- 2
plsda_fin  <- plsda(X, Y, ncomp = ncomp_pls)
perf_pls   <- perf(plsda_fin, validation = "Mfold", folds = 4, nrepeat = 50,
                   BPPARAM = param, progressBar = FALSE)
print(perf_pls$error.rate$BER)
print(perf_pls$error.rate.class$max.dist)


# ==============================================================================
# 3. sPLS-DA
# ==============================================================================

set.seed(123)
tune_spls <- tune.splsda(X = X, Y = Y, ncomp = 8, validation = "Mfold", folds = 10,
                         nrepeat = 50, test.keepX = c(1:10, seq(20, 150, 10)),
                         BPPARAM = param, progressBar = FALSE)
save_png(paste0("sPLSDA_", tag, "_tuning.png"), plot(tune_spls, sd = TRUE))

ncomp_spls <- tune_spls$choice.ncomp$ncomp
keepX_spls <- tune_spls$choice.keepX[seq_len(ncomp_spls)]
cat(sprintf("sPLS-DA: ncomp = %d; keepX = %s\n", ncomp_spls, paste(keepX_spls, collapse = ", ")))

splsda_fin <- splsda(X = X, Y = Y, ncomp = ncomp_spls, keepX = keepX_spls)


# ==============================================================================
# 4. PERFORMANCE, AUC, AND FEATURE STABILITY
# ==============================================================================

set.seed(123)
perf_spls <- perf(splsda_fin, validation = "Mfold", folds = 4, nrepeat = 50, auc = TRUE,
                  BPPARAM = param, progressBar = FALSE)
print(perf_spls$error.rate)
print(perf_spls$error.rate.class)
auc_spls <- auroc(splsda_fin, roc.comp = 2, print = FALSE)   # training data, not cross-validated

saveRDS(list(plsda = plsda_fin, perf_plsda = perf_pls, tune_splsda = tune_spls,
             splsda = splsda_fin, perf_splsda = perf_spls),
        file.path(out_dir, paste0("plsda_splsda_", tag, ".rds")))

save_png(paste0("sPLSDA_", tag, "_feature_stability.png"), width = 10, height = 5, {
  par(mfrow = c(1, 2))
  barplot(perf_spls$features$stable$comp1, xlab = "Features selected across CV folds",
          ylab = "Stability frequency", main = "Component 1")
  barplot(perf_spls$features$stable$comp2, xlab = "Features selected across CV folds",
          ylab = "Stability frequency", main = "Component 2")
})

sel1 <- selectVar(splsda_fin, comp = 1)
stab <- tibble(Feature = sel1$name, loading = sel1$value$value.var,
               stability = as.numeric(perf_spls$features$stable$comp1[sel1$name]))
write.xlsx(stab, file.path(out_dir, paste0("sPLSDA_", tag, "_comp1_features.xlsx")))


# ==============================================================================
# 5. SAMPLE PLOTS AND CIM
# ==============================================================================

for (m in list(list(fit = plsda_fin, name = "PLSDA"), list(fit = splsda_fin, name = "sPLSDA"))) {
  p_zone <- plotIndiv(m$fit, comp = c(1, 2), group = Y, pch = trip, col = unname(fish_cols[levels(Y)]),
                      ind.names = FALSE, legend = TRUE, ellipse = TRUE,
                      legend.title = "Zoning status", legend.title.pch = "Sampling trip",
                      title = paste(m$name, tag))
  ggsave(file.path(out_dir, paste0(m$name, "_", tag, "_by_zone.png")), p_zone$graph, width = 7.5, height = 7.5)
}

tax <- as.data.frame(as(tax_table(inp$ps), "matrix"))[colnames(X), , drop = FALSE]
feature_labels <- with(tax, paste(Class, Order, Family, Genus, Species, sep = "; "))

save_png(paste0("sPLSDA_", tag, "_CIM_by_zone.png"), width = 14, height = 16, {
  cim(splsda_fin, comp = 1:2, col.names = feature_labels, row.names = TRUE,
      margins = c(37, 25), keysize = c(1, 0.4),
      row.sideColors = unname(fish_cols[as.character(Y)]),
      legend = list(legend = levels(Y), col = unname(fish_cols[levels(Y)]),
                    title = "Zoning status", cex = 0.8),
      title = "sPLS-DA, components 1 and 2", xlab = "Selected features", ylab = "Samples")
})


# ==============================================================================
# 6. LEAVE-ONE-REEF-OUT VALIDATION
# ==============================================================================
# Train on all reefs but one, predict the held-out reef (centroids distance)

run_loro <- function(ncomp, keepX) {
  res <- mclapply(unique(reef), function(r) {
    test <- reef == r
    fit  <- splsda(X[!test, , drop = FALSE], Y[!test], ncomp = ncomp, keepX = keepX[seq_len(ncomp)])
    pr   <- predict(fit, X[test, , drop = FALSE], dist = "centroids.dist")
    tibble(Sample = rownames(X)[test], Reef = r, Truth = as.character(Y[test]),
           Pred = pr$class$centroids.dist[, ncomp])
  }, mc.cores = n_cores)
  bind_rows(res)
}

loro_spls <- run_loro(ncomp_spls, keepX_spls)
loro_pls  <- run_loro(ncomp_pls, rep(ncol(X), ncomp_pls))   # keep all features = PLS-DA

loro_summary <- tibble(
  model = c("PLS-DA", "sPLS-DA"),
  balanced_accuracy = c(balanced_accuracy(loro_pls$Truth, loro_pls$Pred),
                        balanced_accuracy(loro_spls$Truth, loro_spls$Pred))
)
print(loro_summary)

write.xlsx(list(summary = loro_summary, PLSDA = loro_pls, sPLSDA = loro_spls),
           file.path(out_dir, paste0("LORO_", tag, ".xlsx")))
