# ==============================================================================
# 2.15_protein_ordination.R
#
# PCoA of normalised RPKM on Aitchison distances, sediment vs seawater.
#   1. Setup
#   2. Ordination function
#   3. KOs (Fig. 1e): each KO-annotated protein as a feature
#   4. CAZymes (Fig. 1f): proteins summed by CAZyme family
#
# Input : ps_ko.RDS and ps_dbcan.RDS (2.14)
# Output: Fig. 1e and 1f (svg)
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

mem.maxVSize(64000); gc()

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(tidyverse)
})

rds_dir <- "<path/to/3_read_mapping_RPKM/RPKM>"   # output of 2.14
out_dir <- "<path/to/3_read_mapping_RPKM/ordination>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

habitat_col  <- "Source"
trip_col     <- "TripDescriotion"   # column name as spelt in the metadata
habitat_cols <- c(Seawater = "#1f78b4", Sediment = "#33a02c")


# ==============================================================================
# 2. ORDINATION FUNCTION
# ==============================================================================
# rank = "unique": each protein is a feature; otherwise proteins are summed by
# that rank before the distance is calculated

plot_pcoa <- function(ps, rank) {
  ps %>%
    tax_transform(trans = "identity", rank = rank) %>%
    dist_calc("aitchison") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(size = 6, fill = habitat_col, shape = trip_col, plot_taxa = 0) +
    stat_ellipse(aes(linetype = .data[[habitat_col]], colour = .data[[habitat_col]]),
                 linewidth = 0.2) +
    scale_colour_manual(values = habitat_cols) +
    scale_fill_manual(values = habitat_cols) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))
}


# ==============================================================================
# 3. KOs (FIG. 1E)
# ==============================================================================

p_ko <- plot_pcoa(readRDS(file.path(rds_dir, "ps_ko.RDS")), rank = "unique")
p_ko
ggsave(file.path(out_dir, "Fig1e_KO_pcoa.svg"), p_ko, width = 7.5, height = 7.5)


# ==============================================================================
# 4. CAZymes (FIG. 1F)
# ==============================================================================

p_cazy <- plot_pcoa(readRDS(file.path(rds_dir, "ps_dbcan.RDS")), rank = "Gene")
p_cazy
ggsave(file.path(out_dir, "Fig1f_CAZyme_pcoa.svg"), p_cazy, width = 7.5, height = 7.5)

# ==============================================================================
# 5. greening db and hydrogenase 
# ==============================================================================

p_greening_hyddb <- plot_pcoa(readRDS(file.path(rds_dir, "ps_greening_hyddb.RDS")), rank = "Gene")
p_greening_hyddb
ggsave(file.path(out_dir, "Fig_GreeningDb_hydrogenase_pcoa.svg"), p_greening_hyddb, width = 7.5, height = 7.5)


# ==============================================================================
sessionInfo()