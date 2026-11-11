# ==============================================================================
# 5.08_marker_gene_enrichment.R
#
# Copy numbers of metabolic marker genes and HydDB hydrogenase subgroups in the
# 796 highly complete MAGs, compared between seawater and sediment (Fig. 3,
# Table S16).
#   1. Setup and input
#   2. Count matrix (MAGs x genes), with database names converted to gene names
#   3. Mann-Whitney U tests with Benjamini-Hochberg correction
#   4. Percentage of MAGs with at least one hydrogenase
#
# The generic [NiFe] and [FeFe] hydrogenase marker databases are removed;
# HydDB subgroups are used instead. Genes with no copies in any MAG are not
# tested.
#
# Input : marker_genes.tsv and hyddb.tsv from 5.06 (Genome, Annotation);
#         gene name lookup (GeneName = database or HydDB name, Gene,
#         GeneMetabolism); MAG table with SOURCE and EnrichM_In (TRUE for the
#         796 MAGs with CheckM2 completeness >= 85%)
# Output: enrichment results, count matrix, and MAG table used (csv)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
})

annot_dir     <- "<path/to/5_MAG_annotation/6_enrichment_tests/annotations>"   # output of 5.06
lookup_file   <- "<path/to/GeneNameUpdater.xlsx>"
metadata_file <- "<path/to/MAG_taxonomy_table.xlsx>"
metadata_sheet <- "curatedtaxonomy_Sed1201_Sea876"
out_dir       <- "<path/to/5_MAG_annotation/6_enrichment_tests>"

g1 <- "Seawater"
g2 <- "Sediment"

lookup <- read.xlsx(lookup_file)
gene_lookup       <- lookup %>% distinct(GeneName, Gene)
metabolism_lookup <- lookup %>% distinct(GeneMetabolism, Gene)

dupe <- gene_lookup %>% count(GeneName) %>% filter(n > 1)
if (nrow(dupe) > 0) warning("GeneName -> Gene is not 1:1 for: ", paste(dupe$GeneName, collapse = ", "))

mag_info <- read.xlsx(metadata_file, sheet = metadata_sheet, rowNames = TRUE) %>%
  filter(as.logical(EnrichM_In)) %>%
  rownames_to_column("MAG")
table(mag_info$SOURCE, useNA = "always")


# ==============================================================================
# 2. COUNT MATRIX
# ==============================================================================

annotations <- bind_rows(
  read_tsv(file.path(annot_dir, "marker_genes.tsv"), show_col_types = FALSE),
  read_tsv(file.path(annot_dir, "hyddb.tsv"),        show_col_types = FALSE)
) %>%
  rename(MAG = Genome, GeneName = Annotation) %>%
  filter(!GeneName %in% c("NiFe_hydrogenase", "FeFe_hydrogenase"))

count_long <- annotations %>%
  filter(MAG %in% mag_info$MAG) %>%
  left_join(gene_lookup, by = "GeneName") %>%
  mutate(Gene = coalesce(Gene, GeneName)) %>%     # keep the original name if not in the lookup
  count(MAG, Gene, name = "count") %>%
  complete(MAG = mag_info$MAG, Gene, fill = list(count = 0)) %>%   # MAGs without a gene have 0 copies
  left_join(mag_info %>% select(MAG, SOURCE), by = "MAG")

count_mat <- count_long %>%
  select(MAG, Gene, count) %>%
  pivot_wider(names_from = Gene, values_from = count, values_fill = 0) %>%
  column_to_rownames("MAG")

cat(g1, ": ", sum(mag_info$SOURCE == g1), " MAGs; ", g2, ": ", sum(mag_info$SOURCE == g2), " MAGs\n", sep = "")


# ==============================================================================
# 3. MANN-WHITNEY U TESTS
# ==============================================================================

mw <- function(x, y) suppressWarnings(wilcox.test(x, y, exact = FALSE))

results <- count_long %>%
  group_by(Gene) %>%
  filter(sum(count) > 0) %>%
  summarise(
    g1_mean = mean(count[SOURCE == g1]),
    g2_mean = mean(count[SOURCE == g2]),
    wt      = list(mw(count[SOURCE == g1], count[SOURCE == g2])),
    .groups = "drop"
  ) %>%
  mutate(
    score    = map_dbl(wt, ~ unname(.x$statistic)),
    pvalue   = map_dbl(wt, ~ .x$p.value),
    `FoldChange_g1/g2` = log2((g1_mean + 1e-10) / (g2_mean + 1e-10)),
    `FoldChange_g2/g1` = -`FoldChange_g1/g2`,
    enriched_in = case_when(g1_mean > g2_mean ~ g1,
                            g1_mean < g2_mean ~ g2,
                            TRUE              ~ NA_character_),
    corrected_pvalue = p.adjust(pvalue, method = "BH")
  ) %>%
  left_join(metabolism_lookup, by = "Gene") %>%
  transmute(gene = Gene,
            GeneMetabolism = coalesce(GeneMetabolism, "Unassigned"),
            group_1 = g1, group_2 = g2, enriched_in,
            group_1_mean = g1_mean, group_2_mean = g2_mean,
            `FoldChange_g1/g2`, `FoldChange_g2/g1`,
            score, pvalue, corrected_pvalue) %>%
  arrange(corrected_pvalue)

write.csv(results,   file.path(out_dir, "marker_genes_mannwhitney.csv"), row.names = FALSE)
write.csv(count_mat, file.path(out_dir, "marker_genes_countmat_UsedInEnrichment.csv"))
write.csv(mag_info,  file.path(out_dir, "marker_genes_maginfo_UsedInEnrichment.csv"), row.names = FALSE)


# ==============================================================================
# 4. MAGs WITH AT LEAST ONE HYDROGENASE
# ==============================================================================
# HydDB subgroups: gene names starting with "[NiFe]", "[FeFe]", or "[Fe]"

hyd_genes <- unique(count_long$Gene[str_detect(count_long$Gene, "^\\[(NiFe|FeFe|Fe)\\]")])

count_long %>%
  filter(Gene %in% hyd_genes) %>%
  group_by(MAG, SOURCE) %>%
  summarise(any_hyd = any(count > 0), .groups = "drop") %>%
  group_by(SOURCE) %>%
  summarise(MAGs = n(), with_hydrogenase = sum(any_hyd), pct = 100 * with_hydrogenase / MAGs)

sessionInfo()
