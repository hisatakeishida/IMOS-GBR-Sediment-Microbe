# ==============================================================================
# 6.05_ordination_dbrda.R
#
# Supporting analysis (not reported): ordination of community composition by
# zoning status, GBR sector, sampling trip, and season, and distance-based
# redundancy analysis (dbRDA).
#   1. Setup and input
#   2. Composition by sampling trip (family level)
#   3. PCA of CLR values coloured by design variables
#   4. dbRDA: distance ~ zoning status + sampling trip + GBR sector
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.05_ordination_dbrda.R <HAB> <METHOD>
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01)
# Output: figures (png) and dbRDA ANOVA table (xlsx)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(vegan)
  library(openxlsx)
  library(patchwork)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
profile_dir <- file.path(in_dir, "0_profiles")   # output of 6.01
out_dir <- file.path(in_dir, "5_ordination_dbRDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col   <- "FishingStatus"
sector_col <- "GBRSector"
trip_col   <- "TripDescription"
n_perm     <- 9999
n_cores    <- max(1, parallel::detectCores() - 2)

fish_cols <- c(Closed = "#41ce80", Open = "#37648b")
trip_cols <- c("Trip 1 | Nov-Dec 2019" = "#f2b09e", "Trip 2 | January 2020" = "#5b5b5b",
               "Trip 3 | February 2020" = "#f3c789", "Trip 4 | July 2020" = "slateblue")
season_cols <- c("Trip 1 | Nov-Dec 2019" = "indianred", "Trip 2 | January 2020" = "indianred4",
                 "Trip 3 | February 2020" = "red3", "Trip 4 | July 2020" = "slateblue")

inp  <- readRDS(file.path(profile_dir, paste0("zoning_input_", tag, ".rds")))
meta <- inp$meta
ps   <- ps_mutate(inp$ps, FishingStatus = factor(FishingStatus, levels = levels(meta[[zone_col]])))


# ==============================================================================
# 2. COMPOSITION BY SAMPLING TRIP
# ==============================================================================

ps_trip <- ps %>%
  tax_transform("compositional", rank = "Family") %>%
  ps_get() %>%
  speedyseq::merge_samples2(trip_col, fun_otu = mean)
p_bar <- comp_barplot(ps_trip, tax_level = "Family", n_taxa = 15) + coord_flip()
ggsave(file.path(out_dir, paste0("composition_by_trip_", tag, ".png")), p_bar, width = 9, height = 5)


# ==============================================================================
# 3. PCA
# ==============================================================================
# PCA on CLR values (pseudocount 1) at the lowest rank (MAG or mOTU)

ord <- ps %>% tax_transform("clr", rank = "unique") %>% ord_calc(method = "PCA")

plot_pca <- function(var, cols = NULL, shape = NULL) {
  p <- ord_plot(ord, fill = var, shape = if (is.null(shape)) var else shape,
                size = 3, alpha = 1, plot_taxa = FALSE) +
    stat_chull(aes(group = .data[[var]], colour = .data[[var]]), geom = "polygon", alpha = 0) +
    scale_shape_manual(values = c(21:25, 7, 12)) +
    theme_bw() + theme(aspect.ratio = 1)
  if (!is.null(cols)) p <- p + scale_fill_manual(values = cols) + scale_colour_manual(values = cols)
  p
}

p_pca <- (plot_pca(zone_col, fish_cols) | plot_pca(sector_col, shape = zone_col)) /
         (plot_pca(trip_col, trip_cols)  | plot_pca(trip_col, season_cols))
ggsave(file.path(out_dir, paste0("PCA_", tag, ".png")), p_pca, width = 15, height = 13)


# ==============================================================================
# 4. dbRDA
# ==============================================================================

dist_aitchison <- dist(inp$clr)
fit_dbrda <- dbrda(reformulate(c(zone_col, trip_col, sector_col), "dist_aitchison"), data = meta)
print(RsquareAdj(fit_dbrda))

set.seed(123)
anova_model <- anova(fit_dbrda, permutations = n_perm)
anova_terms <- anova(fit_dbrda, by = "terms", permutations = n_perm, parallel = n_cores)
print(anova_model)
print(anova_terms)

write.xlsx(list(model = as.data.frame(anova_model), terms = as.data.frame(anova_terms),
                R2 = as.data.frame(RsquareAdj(fit_dbrda))),
           file.path(out_dir, paste0("dbRDA_ANOVA_", tag, ".xlsx")), rowNames = TRUE)

db_pct <- round(100 * fit_dbrda$CCA$eig[1:2] / fit_dbrda$tot.chi, 1)
p_dbrda <- as.data.frame(scores(fit_dbrda, display = "sites", choices = 1:2)) %>%
  mutate(Zone = meta[[zone_col]]) %>%
  ggplot(aes(dbRDA1, dbRDA2, colour = Zone, fill = Zone)) +
  stat_chull(geom = "polygon", alpha = 0.15) +
  geom_point(shape = 21, size = 2, alpha = 0.8) +
  scale_colour_manual(values = fish_cols) + scale_fill_manual(values = fish_cols) +
  labs(x = sprintf("dbRDA1 (%.1f%%)", db_pct[1]), y = sprintf("dbRDA2 (%.1f%%)", db_pct[2]),
       colour = "Zoning status", fill = "Zoning status") +
  theme_bw()
ggsave(file.path(out_dir, paste0("dbRDA_", tag, ".png")), p_dbrda, width = 7.5, height = 7.5)
