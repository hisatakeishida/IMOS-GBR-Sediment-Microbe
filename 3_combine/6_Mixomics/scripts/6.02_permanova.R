# ==============================================================================
# 6.02_permanova.R
#
# Sequential PERMANOVA of community composition (Fig. S3a).
#   1. Setup and input
#   2. PERMANOVA: zoning status alone, and sequential model
#   3. Variance explained per term (figure and table)
#   4. Supporting: reef-level PERMANOVA
#
# Aitchison distances (Euclidean distances of CLR values, pseudocount 1) with
# 9,999 permutations; sequential (type I) terms in the order
#   distance ~ zoning status + sampling trip + GBR sector + reef site
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.02_permanova.R <HAB> <METHOD>
#         HAB = SED or SEA; METHOD = MAG or SINGLEM
# Input : zoning_input_<HAB>_<METHOD>.rds (6.01)
# Output: PERMANOVA tables (xlsx) and variance-explained bar chart (png, svg)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(openxlsx)
  library(tidyverse)
})

args   <- commandArgs(trailingOnly = TRUE)
hab    <- if (length(args) >= 1) args[1] else "SED"
method <- if (length(args) >= 2) args[2] else "MAG"
tag    <- paste0(hab, "_", method)

in_dir  <- "<path/to/6_Mixomics>"
out_dir <- file.path(in_dir, "1_PERMANOVA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_col   <- "FishingStatus"
sector_col <- "GBRSector"
reef_col   <- "ReefName"
trip_col   <- "TripDescription"
n_perm     <- 9999

inp  <- readRDS(file.path(in_dir, paste0("zoning_input_", tag, ".rds")))
meta <- inp$meta
dist_aitchison <- dist(inp$clr)   # Euclidean distance of CLR = Aitchison distance


# ==============================================================================
# 2. PERMANOVA
# ==============================================================================

set.seed(123)
permanova_zone <- adonis2(reformulate(zone_col, "dist_aitchison"),
                          data = meta, permutations = n_perm)
print(permanova_zone)

# Sequential model, zoning status first. Samples are pseudo-replicates of reefs
# for the zoning term, so its F and p are optimistic (see section 4).
set.seed(123)
permanova_full <- adonis2(reformulate(c(zone_col, trip_col, sector_col, reef_col), "dist_aitchison"),
                          data = meta, permutations = n_perm, by = "terms")
print(permanova_full)

write.xlsx(list(zoning_only = as.data.frame(permanova_zone),
                sequential  = as.data.frame(permanova_full)),
           file.path(out_dir, paste0("PERMANOVA_", tag, "_", n_perm, ".xlsx")),
           rowNames = TRUE)


# ==============================================================================
# 3. VARIANCE EXPLAINED PER TERM
# ==============================================================================

plot_permanova_r2 <- function(aov_tab, highlight = zone_col) {
  as.data.frame(aov_tab) %>%
    rownames_to_column("Term") %>%
    filter(Term != "Total") %>%
    mutate(Percent = R2 * 100, Term = fct_reorder(Term, Percent)) %>%
    ggplot(aes(Percent, Term, fill = Term == highlight)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", Percent)), hjust = -0.1, size = 3.5) +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey70"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Variance explained (%)", y = NULL, title = tag) +
    theme_classic()
}

p_permanova <- plot_permanova_r2(permanova_full)
for (ext in c("png", "svg")) {
  ggsave(file.path(out_dir, paste0("PERMANOVA_", tag, "_", n_perm, ".", ext)),
         p_permanova, width = 7.5, height = 4)
}


# ==============================================================================
# 4. SUPPORTING: REEF-LEVEL PERMANOVA
# ==============================================================================
# Mean CLR per reef, so each reef counts once; permutations within GBR sector.
# Not reported in the manuscript.

reef_n    <- table(meta[[reef_col]])
reef_clr  <- rowsum(inp$clr, meta[[reef_col]])
reef_clr  <- reef_clr / as.vector(reef_n[rownames(reef_clr)])
reef_meta <- meta[!duplicated(meta[[reef_col]]), ]
rownames(reef_meta) <- reef_meta[[reef_col]]
reef_meta <- reef_meta[rownames(reef_clr), ]

set.seed(123)
permanova_reef <- adonis2(reformulate(c(zone_col, sector_col, trip_col), "dist(reef_clr)"),
                          data = reef_meta, permutations = n_perm, by = "terms",
                          strata = reef_meta[[sector_col]])
print(permanova_reef)

write.xlsx(as.data.frame(permanova_reef),
           file.path(out_dir, paste0("PERMANOVA_reeflevel_", tag, "_", n_perm, ".xlsx")),
           rowNames = TRUE)
