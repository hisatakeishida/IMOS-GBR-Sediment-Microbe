# ==============================================================================
# 2.16_protein_differential_abundance.R
#
# Differential abundance of metabolic marker genes (including HydDB subgroups)
# between sediment and seawater (Fig. 3c).
#   1. Setup
#   2. Input tables (genes as features, relative abundance)
#   3. MaAsLin2
#   4. LinDA
#   5. Combined table and consensus
#   6. Run
#
# MaAsLin2: TSS normalisation, AST transformation, fixed effects habitat and
#           season, reef as random effect, features in >= 5% of samples
# LinDA   : proportions, ~ habitat + season + (1 | reef), features in >= 5% of samples
# Sediment is the reference level, so positive coefficients mean enrichment in
# seawater. A gene is differentially abundant if both methods give adjusted
# p <= 0.05 for habitat with the same sign.
#
# Input : ps_greening_hyddb.RDS (2.14)
# Output: MAASLIN_LINDA_greeningDB.xlsx
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(Maaslin2)
  library(MicrobiomeStat)   # linda()
  library(openxlsx)
  library(tidyverse)
})

rds_dir <- "<path/to/3_read_mapping_RPKM/RPKM>"   # output of 2.14
out_dir <- "<path/to/4_differential_abundance>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rank_DA        <- "Gene"
habitat_col    <- "Source"
season_col     <- "TripSeason"
random_col     <- "Reef"
habitat_levels <- c("Sediment", "Seawater")   # reference level first
min_prev       <- 0.05
q_cutoff       <- 0.05


# ==============================================================================
# 2. INPUT TABLES
# ==============================================================================

get_metadata <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  meta[[habitat_col]] <- factor(as.character(meta[[habitat_col]]), levels = habitat_levels)
  meta[[season_col]]  <- factor(as.character(meta[[season_col]]))
  meta
}

# Samples x genes, relative abundance within the gene set (MaAsLin2 takes
# samples x features; LinDA takes features x samples, see run_linda)
gene_table <- function(ps) {
  ps %>%
    tax_transform(trans = "compositional", rank = rank_DA) %>%
    otu_table() %>% t() %>% as.data.frame()
}

# Mean relative abundance per habitat (feature names as returned by MaAsLin2)
mean_by_habitat <- function(tab, meta) {
  tab %>%
    rownames_to_column("SampleID") %>%
    mutate(Habitat = meta[SampleID, habitat_col]) %>%
    pivot_longer(-c(SampleID, Habitat), names_to = "feature", values_to = "RelAbundance") %>%
    group_by(feature, Habitat) %>%
    summarise(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Habitat, values_from = MeanRelAbundance,
                names_prefix = "MeanRelAbund_", values_fill = 0) %>%
    mutate(feature = make.names(feature))
}


# ==============================================================================
# 3. MaAsLin2
# ==============================================================================

run_maaslin <- function(tab, meta, label) {
  fit <- Maaslin2(
    input_data      = tab,
    input_metadata  = meta,
    output          = file.path(out_dir, paste0(label, "_TSS_AST_", habitat_col, "_",
                                                season_col, "_", random_col, "asRandom")),
    fixed_effects   = c(habitat_col, season_col),
    random_effects  = random_col,
    normalization   = "TSS",
    transform       = "AST",
    min_prevalence  = min_prev,
    heatmap_first_n = 100
  )
  res <- fit$results

  habitat_res <- res %>%
    filter(metadata == habitat_col, value == habitat_levels[2]) %>%
    select(feature, maaslin_SedSea_coef = coef, maaslin_SedSea_p = pval, maaslin_SedSea_q = qval)

  season_res <- res %>%
    filter(metadata == season_col) %>%
    select(feature, maaslin_TripSeason_coef = coef, maaslin_TripSeason_p = pval,
           maaslin_TripSeason_q = qval)

  full_join(habitat_res, season_res, by = "feature")
}


# ==============================================================================
# 4. LinDA
# ==============================================================================

run_linda <- function(tab, meta) {
  fit <- linda(t(tab), meta,
               formula          = paste0("~", habitat_col, " + ", season_col, " + (1|", random_col, ")"),
               feature.dat.type = "proportion",
               prev.filter      = min_prev)

  habitat_term <- paste0(habitat_col, habitat_levels[2])          # e.g. SourceSeawater
  season_term  <- grep(paste0("^", season_col), names(fit$output), value = TRUE)

  tibble(feature                = make.names(rownames(fit$output[[habitat_term]])),
         linda_SedSea_coef      = fit$output[[habitat_term]]$log2FoldChange,
         linda_SedSea_p         = fit$output[[habitat_term]]$pvalue,
         linda_SedSea_q         = fit$output[[habitat_term]]$padj,
         linda_TripSeason_coef  = fit$output[[season_term]]$log2FoldChange,
         linda_TripSeason_p     = fit$output[[season_term]]$pvalue,
         linda_TripSeason_q     = fit$output[[season_term]]$padj)
}


# ==============================================================================
# 5. COMBINED TABLE AND CONSENSUS
# ==============================================================================

run_da <- function(ps, label) {
  meta <- get_metadata(ps)
  tab  <- gene_table(ps)

  run_maaslin(tab, meta, label) %>%
    left_join(mean_by_habitat(tab, meta), by = "feature") %>%
    full_join(run_linda(tab, meta), by = "feature") %>%
    mutate(
      consensus   = maaslin_SedSea_q <= q_cutoff & linda_SedSea_q <= q_cutoff &
                    sign(maaslin_SedSea_coef) == sign(linda_SedSea_coef),
      enriched_in = case_when(!consensus             ~ NA_character_,
                              linda_SedSea_coef > 0  ~ habitat_levels[2],
                              TRUE                   ~ habitat_levels[1])
    ) %>%
    arrange(desc(consensus), maaslin_SedSea_q)
}


# ==============================================================================
# 6. RUN
# ==============================================================================

# Metabolic marker genes and HydDB subgroups (Fig. 3c)
da_greening <- run_da(readRDS(file.path(rds_dir, "ps_greening_hyddb.RDS")), "greeningDB")
write.xlsx(da_greening, file.path(out_dir, "MAASLIN_LINDA_greeningDB.xlsx"))

table(da_greening$enriched_in, useNA = "ifany")


# ==============================================================================
sessionInfo()