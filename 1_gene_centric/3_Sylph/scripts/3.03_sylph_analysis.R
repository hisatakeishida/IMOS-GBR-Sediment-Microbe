#!/usr/bin/env Rscript
# =============================================================================
# 3.03_sylph_analysis.R
#
# Species-level sylph profiles for IMOS (seawater and sediment):
#   1. Build a phyloseq object from the sylph species table and SingleM
#      environmental metadata
#   2. PCA of CLR-transformed profiles
#   3. Alpha diversity (Shannon) by habitat
#   4. Composition barplots at several ranks
#
# Input:
#   sylph_taxonomy_abundances_sequence_abundance.xlsx  (from sylph_to_excel.R)
#   IMOS_10Mseqs_SingleM_FINAL.xlsx, sheet "seased_env_RS226"
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(vegan)
  library(Maaslin2)
  library(MicrobiomeStat)   # linda()
  library(data.table)       # used by merge_samples2()
  library(openxlsx)
  library(tidyverse)        # loads stringr, dplyr, ggplot2
  library(patchwork)
})

# ---------------------------------------------------------------------------
# 0.  Paths and parameters
# ---------------------------------------------------------------------------
singlem_xlsx <- "<path/to/IMOS_10Mseqs_SingleM_FINAL.xlsx>"
sylph_xlsx   <-"<path/to/sylph_taxonomy_abundances_sequence_abundance.xlsx>"
out_dir      <-"<path/to/3_Sylph>"

rank_used       <- "Family"
min_prevalence  <- 0.1    # DA: taxa detected in >= 10% of all samples
core_prevalence <- 0.8    # core taxa: detected in > 80% of samples
q_cutoff        <- 0.05
n_top_taxa      <- 15     # composition barplots
bar_ranks       <- c("Phylum", "Class", "Order", "Family", "Genus")

habitat_col    <- "Source"
reef_col       <- "ReefName"
trip_col       <- "TripDescription"
habitat_levels <- c("Seawater", "Sediment")   # reference level first
habitat_cols   <- c(Seawater = "#1f78b4", Sediment = "#33a02c")

save_fig <- function(p, file, width = 8, height = 6) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, file), p, width = width, height = height)
}

# ---------------------------------------------------------------------------
# 1.  Read inputs
# ---------------------------------------------------------------------------
sylph_meta <- read.xlsx(singlem_xlsx, sheet = "seased_env_RS226", rowNames = TRUE)
sylph_otu  <- read.xlsx(sylph_xlsx,   sheet = "taxonomy_species", rowNames = TRUE)

sylph_meta[[habitat_col]] <- factor(sylph_meta[[habitat_col]], levels = habitat_levels)

# ---------------------------------------------------------------------------
# 2.  Build phyloseq
# ---------------------------------------------------------------------------
# OTU table (species x samples)
otu_mat <- as.matrix(sylph_otu[rownames(sylph_otu) != "Unassigned", , drop = FALSE])
storage.mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0

# Taxonomy table: split "d__;p__;c__;o__;f__;g__;s__" lineages
tax_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sylph_tax <- str_split_fixed(rownames(otu_mat), ";", length(tax_ranks))
sylph_tax[] <- str_remove(sylph_tax, "^[a-z]__")   # drop rank prefixes
sylph_tax[sylph_tax == ""] <- NA                   # empty ranks (e.g. "g__") to NA
colnames(sylph_tax) <- tax_ranks

# Short taxa names: species name (unique in GTDB; make.unique as a guard)
taxa_ids <- make.unique(sylph_tax[, "Species"])
rownames(otu_mat) <- rownames(sylph_tax) <- taxa_ids

# Match samples to metadata
shared <- intersect(colnames(otu_mat), rownames(sylph_meta))
message(length(shared), " shared samples; ",
        length(setdiff(colnames(otu_mat), rownames(sylph_meta))), " sylph-only; ",
        length(setdiff(rownames(sylph_meta), colnames(otu_mat))), " metadata-only")

ps_sylph <- phyloseq(
  otu_table(otu_mat[, shared, drop = FALSE], taxa_are_rows = TRUE),
  tax_table(sylph_tax),
  sample_data(sylph_meta[shared, , drop = FALSE])
)
ps_sylph <- ps_sylph %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  tax_fix()   # fill NA ranks so microViz can aggregate at any rank
ps_sylph

# ---------------------------------------------------------------------------
# 3.  PCA of CLR profiles at rank_used
# ---------------------------------------------------------------------------
p_habitat_pca <- ps_sylph %>%
  tax_transform(trans = "clr", rank = rank_used) %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    size = 6, fill = habitat_col, shape = trip_col, plot_taxa = 10:1,
    tax_vec_style_all = vec_constraint(linewidth = 0, arrow = NULL)
  ) +
  stat_ellipse(aes(colour = .data[[habitat_col]]), linewidth = 0.2) +
  scale_fill_manual(values = habitat_cols) +
  scale_colour_manual(values = habitat_cols) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  theme_bw() +
  theme(axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "none")
p_habitat_pca
save_fig(p_habitat_pca, paste0("PCA_habitat_", tolower(rank_used), ".pdf"))

# ---------------------------------------------------------------------------
# 4.  Alpha diversity: Shannon index at rank_used
# ---------------------------------------------------------------------------
shannon_col <- paste0("shannon_", rank_used)

alpha_tbl <- ps_sylph %>%
  ps_calc_diversity(rank = rank_used, index = "shannon") %>%
  samdat_tbl()

alpha_summary <- alpha_tbl %>%
  group_by(.data[[habitat_col]]) %>%
  summarise(mean_Shannon = mean(.data[[shannon_col]], na.rm = TRUE),
            sd_Shannon   = sd(.data[[shannon_col]],   na.rm = TRUE))
alpha_summary

wilcox.test(reformulate(habitat_col, shannon_col), data = alpha_tbl)

p_alpha <- ggplot(alpha_tbl, aes(y = .data[[habitat_col]], x = .data[[shannon_col]],
                                 fill = .data[[habitat_col]])) +
  geom_boxplot(width = 0.3) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  scale_fill_manual(values = habitat_cols) +
  labs(x = paste0("Shannon index (", tolower(rank_used), " level)"), y = NULL) +
  theme_bw() +
  theme(legend.position = "none")
p_alpha
save_fig(p_alpha, paste0("Shannon_habitat_", tolower(rank_used), ".pdf"), height = 3)

# ---------------------------------------------------------------------------
# 5.  Composition barplots
# ---------------------------------------------------------------------------
p_comp <- list()

for (rank in bar_ranks) {
  p <- ps_sylph %>%
    ps_arrange(.data[[habitat_col]], .data[[trip_col]]) %>%
    comp_barplot(tax_level = rank, n_taxa = n_top_taxa, merge_other = FALSE,
                 sample_order = "asis", label = NULL, bar_width = 1) +
    facet_grid(rows = vars(.data[[habitat_col]], .data[[trip_col]]),
               scales = "free", space = "free") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0))

  print(p)
  p_comp[[rank]] <- p
  save_fig(p, paste0("composition_", tolower(rank), ".pdf"), width = 8, height = 12)
}