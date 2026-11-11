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

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")   # edit as needed

p_list <- list()

for (rank in ranks) {
  p <- ps %>%
    ps_arrange(.data[[habitat_col]], .data[[trip_col]]) %>%
    comp_barplot(tax_level = rank, n_taxa = n_top_class, merge_other = FALSE,
                 sample_order = "asis", label = NULL, bar_width = 1) +
    facet_grid(rows = vars(.data[[habitat_col]], .data[[trip_col]]),
               scales = "free", space = "free") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0))

  print(p)
  p_list[[rank]] <- p

  ggsave(file.path(out_dir, paste0("Fig1c_",
                                   tolower(rank), "_composition.pdf")),
         p, width = 8, height = 12)
  dev.off()
}

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
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS 15.7.9

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# time zone: Australia/Brisbane
# tzcode source: internal

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] MicrobiomeStat_1.2 Maaslin2_1.16.0    patchwork_1.3.1    lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
#  [8] purrr_1.1.0        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       ggplot2_3.5.2      tidyverse_2.0.0    openxlsx_4.2.8    
# [15] vegan_2.7-1        permute_0.9-8      microViz_0.12.7    phyloseq_1.46.0    data.table_1.17.8 

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3          farver_2.1.2            nloptr_2.2.1           
#   [7] zlibbioc_1.48.2         ragg_1.4.0              vctrs_0.6.5             multtest_2.58.0         minqa_1.2.8             RCurl_1.98-1.17        
#  [13] Rhdf5lib_1.24.2         rhdf5_2.46.1            plyr_1.8.9              commonmark_2.0.0        igraph_2.1.4            lifecycle_1.0.4        
#  [19] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.6-5            R6_2.6.1                GenomeInfoDbData_1.2.11 rbibutils_2.3          
#  [25] clue_0.3-66             digest_0.6.37           numDeriv_2016.8-1.1     colorspace_2.1-1        spatial_7.3-16          S4Vectors_0.40.2       
#  [31] textshaping_1.0.1       labeling_0.4.3          timechange_0.3.0        mgcv_1.8-42             compiler_4.3.1          withr_3.0.2            
#  [37] DBI_1.2.3               biglm_0.9-3             MASS_7.3-60             biomformat_1.30.0       fBasics_4041.97         optparse_1.7.5         
#  [43] tools_4.3.1             ape_5.8-1               zip_2.3.3               glue_1.8.0              stabledist_0.7-2        nlme_3.1-162           
#  [49] rhdf5filters_1.14.1     gridtext_0.1.5          grid_4.3.1              Rtsne_0.17              cluster_2.1.8.1         reshape2_1.4.4         
#  [55] ade4_1.7-23             generics_0.1.4          gtable_0.3.6            microbiome_1.24.0       tzdb_0.5.0              hms_1.1.3              
#  [61] xml2_1.3.8              utf8_1.2.6              XVector_0.42.0          rmutil_1.1.10           BiocGenerics_0.48.1     ggrepel_0.9.6          
#  [67] foreach_1.5.2           pillar_1.11.0           markdown_2.0            logging_0.10-108        robustbase_0.99-4-1     splines_4.3.1          
#  [73] getopt_1.20.4           ggtext_0.1.2            lattice_0.22-7          survival_3.5-5          tidyselect_1.2.1        pbapply_1.7-4          
#  [79] Biostrings_2.70.3       reformulas_0.4.1        litedown_0.7            IRanges_2.36.0          stats4_4.3.1            xfun_0.52              
#  [85] Biobase_2.62.0          statmod_1.5.0           timeDate_4041.110       matrixStats_1.5.0       DEoptimR_1.1-4          stringi_1.8.7          
#  [91] boot_1.3-28.1           codetools_0.2-19        timeSeries_4041.111     hash_2.2.6.3            cli_3.6.5               rpart_4.1.19           
#  [97] systemfonts_1.2.3       Rdpack_2.6.4            dichromat_2.0-0.1       Rcpp_1.1.0              GenomeInfoDb_1.38.8     stable_1.1.6           
# [103] parallel_4.3.1          modeest_2.4.0           bitops_1.0-9            lme4_1.1-37             mvtnorm_1.3-3           lmerTest_3.1-3         
# [109] scales_1.4.0            statip_0.2.3            pcaPP_2.0-5             crayon_1.5.3            rlang_1.1.6            