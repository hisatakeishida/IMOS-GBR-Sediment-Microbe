# ==============================================================================
# 4.11_novelty_upset.R
#
# UpSet plot of species clusters (95% ANI) shared between GBR sediment and
# seawater MAGs and public genome collections (Fig. 2d).
#
# Input : <prefix>.cluster_summary.tsv from 4.10_novelty_cluster.sh
#         (one row per cluster, number of genomes from each collection)
# Output: Fig. 2d (svg)
# ==============================================================================

suppressPackageStartupMessages({
  library(ComplexUpset)
  library(tidyverse)
})

summary_path <- "<path/to/4_MAG_recovery/9_novelty>/sed_sea_public_ANI95_AF15.cluster_summary.tsv"
out_file     <- "<path/to/4_MAG_recovery/9_novelty>/Fig2d_novelty_upset.svg"

# Presence (1) or absence (0) of each collection in each cluster
df_binary <- read_tsv(summary_path, show_col_types = FALSE) %>%
  column_to_rownames("Rep") %>%
  mutate(across(everything(), ~ as.integer(.x > 0)))

colSums(df_binary)

# Intersections shown, in this order
intersections <- list(
  c("gbr_sediment"),
  c("gbr_sediment", "gbr_seawater"),
  c("gbr_sediment", "china_sea_reef_sediment"),
  c("gbr_sediment", "omd"),
  c("gbr_sediment", "coastal_sediment"),
  c("gbr_sediment", "gbr_seawater", "omd"),
  c("gbr_sediment", "gbr_seawater", "gems", "omd"),
  c("gbr_sediment", "gems"),
  c("gbr_seawater"),
  c("gbr_seawater", "omd"),
  c("gbr_seawater", "gems")
)

p_upset <- upset(df_binary, colnames(df_binary),
                 intersections       = intersections,
                 set_sizes           = FALSE,
                 sort_intersections  = FALSE,
                 sort_sets           = FALSE)
p_upset

ggsave(out_file, p_upset, width = 7.5, height = 7.5)

sessionInfo()