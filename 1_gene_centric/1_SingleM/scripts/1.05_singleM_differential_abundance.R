# ==============================================================================
# 1.5_singleM_differential_abundance.R
#
# Family-level abundance, prevalence, and differential abundance of sediment and
# seawater communities (mOTU-based).
#   1. Setup
#   2. Input: phyloseq object and family-level relative abundance
#   3. Abundance and prevalence per habitat (Tables S5 and S6)
#   4. Composition of the top 15 classes (Fig. 1c)
#   5. Differential abundance: MaAsLin2
#   6. Differential abundance: LinDA
#   7. Consensus table (Table S4)
#   8. Plot consensus families
#
# A family is considered differentially abundant if MaAsLin2 and LinDA both give
# adjusted p <= 0.05 with the same direction of change (reef site as a random
# effect in both).
#
# Input : singlem_prok_SEASED.rds (prokaryotic SingleM phyloseq object, both
#         habitats; written by 1.4_singlem_community_analysis.R)
# Output: Tables S4 to S6, Fig. 1c, and per-method results in out_dir
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

# BiocManager::install("Maaslin2")
# install.packages("MicrobiomeStat")   # LinDA

library(phyloseq)
library(microViz)
library(Maaslin2)
library(MicrobiomeStat)   # linda()
library(openxlsx)
library(tidyverse)

select <- dplyr::select

in_rds  <- "<path/to/1_SingleM_ANALYSIS>/singlem_prok_SEASED.rds"
out_dir <- "<path/to/1_SingleM_ANALYSIS/abundance_DA>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rank_used       <- "Family"
min_prevalence  <- 0.1    # DA: families detected in >= 10% of all samples
core_prevalence <- 0.8    # Tables S5/S6: families detected in > 80% of samples
q_cutoff        <- 0.05
n_top_class     <- 15     # Fig. 1c
habitat_col     <- "Source"
reef_col        <- "ReefName"
trip_col        <- "TripDescription"
habitat_levels  <- c("Seawater", "Sediment")   # reference level first
habitat_cols    <- c(Seawater = "#1f78b4", Sediment = "#33a02c")


# ==============================================================================
# 2. INPUT
# ==============================================================================

ps <- readRDS(in_rds) %>% phyloseq_validate()
sample_data(ps)[[habitat_col]] <- factor(sample_data(ps)[[habitat_col]],
                                         levels = habitat_levels)

# Family-level relative abundance (samples x families) for one habitat
family_relabu <- function(habitat) {
  ps %>%
    ps_filter(.data[[habitat_col]] == habitat) %>%
    tax_transform(trans = "compositional", rank = rank_used) %>%
    otu_get() %>%
    as("matrix")
}


# ==============================================================================
# 3. ABUNDANCE AND PREVALENCE PER HABITAT (TABLES S5 AND S6)
# ==============================================================================

abundance_prevalence <- function(habitat) {
  ra <- family_relabu(habitat)
  tibble(Family         = colnames(ra),
         MeanRelAbu_pct = 100 * colMeans(ra),
         SD_pct         = 100 * apply(ra, 2, sd),
         Prevalence_pct = 100 * colMeans(ra > 0)) %>%
    mutate(Prevalent_80pct = Prevalence_pct > 100 * core_prevalence) %>%
    arrange(desc(MeanRelAbu_pct))
}

abund_sed <- abundance_prevalence("Sediment")
abund_sea <- abundance_prevalence("Seawater")

write.xlsx(abund_sed, file.path(out_dir, "TableS5_family_abundance_prevalence_sediment_mOTU.xlsx"))
write.xlsx(abund_sea, file.path(out_dir, "TableS6_family_abundance_prevalence_seawater_mOTU.xlsx"))

# Families prevalent in > 80% of samples in each habitat
message("Sediment: ", sum(abund_sed$Prevalent_80pct), " families in > 80% of samples")
message("Seawater: ", sum(abund_sea$Prevalent_80pct), " families in > 80% of samples")


# ==============================================================================
# 4. COMPOSITION OF THE TOP 15 CLASSES (FIG. 1C)
# ==============================================================================
# Samples ordered by habitat and sampling trip

p_class <- ps %>%
  ps_arrange(.data[[habitat_col]], .data[[trip_col]]) %>%
  comp_barplot(tax_level = "Class", n_taxa = n_top_class, merge_other = FALSE,
               sample_order = "asis", label = NULL, bar_width = 1) +
  facet_grid(rows = vars(.data[[habitat_col]], .data[[trip_col]]),
             scales = "free", space = "free") +
  coord_flip() +
  theme(axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0))
p_class

ggsave(file.path(out_dir, "Fig1c_class_composition.pdf"), p_class, width = 8, height = 12)


# ==============================================================================
# 5. DIFFERENTIAL ABUNDANCE: MaAsLin2
# ==============================================================================

ps_fam <- ps %>%
  tax_transform(trans = "compositional", rank = rank_used) %>%
  tax_filter(min_prevalence = min_prevalence, undetected = 0, verbose = FALSE)

# Family names contain characters ("-", spaces, leading digits) that MaAsLin2
# and LinDA rewrite. Use safe IDs for the models and map back afterwards.
fam_ids <- tibble(Family = taxa_names(ps_fam),
                  ID     = sprintf("F%04d", seq_along(taxa_names(ps_fam))))

abund <- as(otu_table(ps_fam), "matrix")
if (taxa_are_rows(ps_fam)) abund <- t(abund)      # samples x families
colnames(abund) <- fam_ids$ID[match(colnames(abund), fam_ids$Family)]

meta <- as(sample_data(ps_fam), "data.frame")
stopifnot(identical(rownames(abund), rownames(meta)))

message(nrow(abund), " samples; ", ncol(abund), " families after prevalence filter")

# Default normalisation (TSS) and transformation (LOG) on relative abundance
fit_maaslin <- Maaslin2(
  input_data      = as.data.frame(abund),
  input_metadata  = meta,
  output          = file.path(out_dir, "maaslin2_family_habitat_reef"),
  fixed_effects   = habitat_col,
  random_effects  = reef_col,
  reference       = paste0(habitat_col, ",", habitat_levels[1]),
  min_prevalence  = 0,     # already filtered above
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE
)

maaslin_res <- fit_maaslin$results %>%
  filter(metadata == habitat_col) %>%
  transmute(ID = feature,
            maaslin2_coef = coef,
            maaslin2_q    = qval)


# ==============================================================================
# 6. DIFFERENTIAL ABUNDANCE: LinDA
# ==============================================================================

fit_linda <- linda(
  feature.dat      = t(abund),             # families x samples
  meta.dat         = meta,
  formula          = paste0("~", habitat_col, " + (1|", reef_col, ")"),
  feature.dat.type = "proportion"
)

linda_term <- paste0(habitat_col, habitat_levels[2])   # e.g. "SourceSediment"
linda_res <- fit_linda$output[[linda_term]] %>%
  rownames_to_column("ID") %>%
  transmute(ID,
            linda_log2FC = log2FoldChange,
            linda_padj   = padj)


# ==============================================================================
# 7. CONSENSUS TABLE (TABLE S4)
# ==============================================================================

da_table <- fam_ids %>%
  left_join(maaslin_res, by = "ID") %>%
  left_join(linda_res,   by = "ID") %>%
  left_join(abund_sed %>% select(Family, RelAbu_Sediment_pct = MeanRelAbu_pct), by = "Family") %>%
  left_join(abund_sea %>% select(Family, RelAbu_Seawater_pct = MeanRelAbu_pct), by = "Family") %>%
  mutate(
    maaslin2_sig = maaslin2_q <= q_cutoff,
    linda_sig    = linda_padj <= q_cutoff,
    same_sign    = sign(maaslin2_coef) == sign(linda_log2FC),
    consensus    = maaslin2_sig & linda_sig & same_sign,
    enriched_in  = case_when(!consensus        ~ NA_character_,
                             linda_log2FC > 0  ~ habitat_levels[2],
                             TRUE              ~ habitat_levels[1])
  ) %>%
  select(-ID) %>%
  arrange(desc(consensus), desc(abs(linda_log2FC)))

write.xlsx(da_table, file.path(out_dir, "TableS4_DA_family_sediment_vs_seawater.xlsx"))

table(da_table$enriched_in, useNA = "ifany")


# ==============================================================================
# 8. PLOT CONSENSUS FAMILIES
# ==============================================================================
# Families resolved only to a higher rank (e.g. "X Order") are excluded

consensus_fams <- da_table %>%
  filter(consensus, !grepl(" (Phylum|Class|Order)$", Family)) %>%
  arrange(desc(linda_log2FC)) %>%
  pull(Family)

plot_df <- ps %>%
  tax_transform(trans = "compositional", rank = rank_used) %>%
  tax_select(consensus_fams, ranks_searched = rank_used, strict_matches = TRUE) %>%
  psmelt() %>%
  mutate(OTU = factor(OTU, levels = consensus_fams))

p_da <- ggplot(plot_df, aes(x = .data[[habitat_col]], y = 100 * Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = .data[[habitat_col]]), height = 0, width = 0.2, alpha = 0.5) +
  scale_colour_manual(values = habitat_cols) +
  labs(x = NULL, y = "Relative abundance (%)") +
  facet_wrap(~ OTU, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")
p_da

ggsave(file.path(out_dir, "DA_consensus_families.pdf"), p_da, width = 12, height = 10)


# ==============================================================================
sessionInfo()