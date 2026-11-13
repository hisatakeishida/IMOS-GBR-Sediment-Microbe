# ==============================================================================
# 2.15_protein_ordination.R
#
# Functional profiles of sediment and seawater: PCoA of normalised RPKM on
# Aitchison distances, and Shannon diversity.
#   1. Setup
#   2. Functions
#   3. KOs (Fig. 1e): each KO-annotated protein as a feature
#   4. CAZymes (Fig. 1f): proteins summed by CAZyme family
#   5. Supporting: metabolic marker genes and HydDB subgroups
#   6. Supporting: Shannon diversity by habitat
#
# Input : ps_ko.RDS, ps_dbcan.RDS, ps_greening_hyddb.RDS (2.14)
# Output: Fig. 1e and 1f and supporting PCoA (svg); Shannon diversity by
#         habitat (svg, tsv)
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

mem.maxVSize(64000)   # macOS only: raises the vector memory limit (MB) for the KO object

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

save_svg <- function(plot, name, width = 7.5, height = 7.5) {
  ggsave(file.path(out_dir, name), plot, width = width, height = height, bg = "white")
}


# ==============================================================================
# 2. FUNCTIONS
# ==============================================================================

# PCoA on Aitchison distances. rank = "unique": each protein is a feature;
# otherwise proteins are summed by that rank first
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

# Shannon index per sample at the given rank
shannon_by_habitat <- function(ps, rank, profile) {
  ps %>%
    ps_calc_diversity(rank = rank, index = "shannon", varname = "Shannon") %>%
    samdat_tbl() %>%
    transmute(.sample_name, Habitat = .data[[habitat_col]], Shannon, Profile = profile)
}


# ==============================================================================
# 3. KOs (Fig. 1e)
# ==============================================================================

ps_ko <- readRDS(file.path(rds_dir, "ps_ko.RDS"))
save_svg(plot_pcoa(ps_ko, rank = "unique"), "Fig1e_KO_pcoa.svg")


# ==============================================================================
# 4. CAZymes (Fig. 1f)
# ==============================================================================

ps_cazy <- readRDS(file.path(rds_dir, "ps_dbcan.RDS"))
save_svg(plot_pcoa(ps_cazy, rank = "Gene"), "Fig1f_CAZyme_pcoa.svg")


# ==============================================================================
# 5. SUPPORTING: METABOLIC MARKER GENES AND HydDB SUBGROUPS
# ==============================================================================

ps_greening_hyddb <- readRDS(file.path(rds_dir, "ps_greening_hyddb.RDS"))
save_svg(plot_pcoa(ps_greening_hyddb, rank = "Gene"), "GreeningDB_hydrogenase_pcoa.svg")


# ==============================================================================
# 6. SUPPORTING: SHANNON DIVERSITY BY HABITAT
# ==============================================================================
# KOs: each KO-annotated protein as a feature; CAZymes and marker genes:
# proteins summed by family or gene. Wilcoxon rank-sum test, sediment vs seawater.

alpha_tbl <- bind_rows(
  shannon_by_habitat(ps_ko,             "unique", "KO proteins"),
  shannon_by_habitat(ps_cazy,           "Gene",   "CAZyme families"),
  shannon_by_habitat(ps_greening_hyddb, "Gene",   "Marker genes and HydDB subgroups")
) %>%
  mutate(Profile = fct_inorder(Profile))

alpha_summary <- alpha_tbl %>%
  group_by(Profile) %>%
  summarise(
    mean_sediment = mean(Shannon[Habitat == "Sediment"]), sd_sediment = sd(Shannon[Habitat == "Sediment"]),
    mean_seawater = mean(Shannon[Habitat == "Seawater"]), sd_seawater = sd(Shannon[Habitat == "Seawater"]),
    wilcox_p      = wilcox.test(Shannon ~ Habitat)$p.value,
    .groups = "drop"
  )
print(alpha_summary)
write_tsv(alpha_summary, file.path(out_dir, "Shannon_habitat_summary.tsv"))

p_alpha <- ggplot(alpha_tbl, aes(Shannon, Habitat, fill = Habitat)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(position = position_jitter(height = 0.2, seed = 1), alpha = 0.5) +
  facet_wrap(~ Profile, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = habitat_cols) +
  labs(x = "Shannon index", y = NULL) +
  theme_bw() +
  theme(legend.position = "none")
save_svg(p_alpha, "Shannon_habitat.svg", width = 6, height = 6)

sessionInfo()
