# ==============================================================================
# 2.14_protein_annotation_objects.R
#
# Builds phyloseq objects of normalised RPKM for each annotation, with a common
# taxonomy (GeneMetabolism, Gene), from the KO object of 2.13.
#   1. Setup
#   2. Input: KO RPKM object and updated sample metadata
#   3. Annotation tables
#   4. Phyloseq objects
#
# Marker genes: HydDB subgroups replace the generic [NiFe] and [FeFe]
# hydrogenase marker genes, and HydDB takes priority for proteins hit by both.
#
# Input : ko.RPKM.rds (2.13); annotation tables (2.8); marker gene lookup
#         (database name -> Gene, GeneMetabolism); sample metadata
# Output: ps_ko.RDS, ps_greening_hyddb.RDS, ps_dbcan.RDS
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(openxlsx)
  library(tidyverse)
})

rpkm_rds      <- "<path/to/3_read_mapping_RPKM/RPKM>/ko.RPKM.rds"
annot_dir     <- "<path/to/2_annotation/summary>"
marker_lookup <- "<path/to/marker_gene_lookup.tsv>"   # columns: marker_gene, Gene, GeneMetabolism
metadata_path <- "<path/to/genecentric_kegg_cazy.xlsx>"
out_dir       <- "<path/to/3_read_mapping_RPKM/RPKM>"


# ==============================================================================
# 2. INPUT
# ==============================================================================

kos_ps <- readRDS(rpkm_rds)

new_meta <- sample_data(read.xlsx(metadata_path, rowNames = TRUE))
stopifnot(all(sample_names(kos_ps) %in% sample_names(new_meta)))
sample_data(kos_ps) <- new_meta[sample_names(kos_ps), ]


# ==============================================================================
# 3. ANNOTATION TABLES (GeneID, GeneMetabolism, Gene)
# ==============================================================================

# KOs
ko_annot <- as.data.frame(tax_table(kos_ps)) %>%
  rownames_to_column("GeneID") %>%
  filter(ko_id != "Unassigned") %>%
  transmute(GeneID, GeneMetabolism = "KO", Gene = ko_id)

# Metabolic marker genes, without the generic hydrogenase databases
annot_greening <- read_tsv(file.path(annot_dir, "marker_genes_annot.tsv"), show_col_types = FALSE) %>%
  left_join(read_tsv(marker_lookup, show_col_types = FALSE), by = "marker_gene") %>%
  filter(!Gene %in% c("NiFe_hydrogenase", "FeFe_hydrogenase")) %>%
  select(GeneID, GeneMetabolism, Gene)

# HydDB subgroups
annot_hyddb <- read_tsv(file.path(annot_dir, "hyddb_annot.tsv"), show_col_types = FALSE) %>%
  transmute(GeneID, GeneMetabolism = "Trace gas metabolism", Gene = hyddb_id)

# Marker genes and HydDB, HydDB first
annot_greening_hyddb <- bind_rows(
  annot_greening %>% filter(!GeneID %in% annot_hyddb$GeneID),
  annot_hyddb
)

# CAZyme families (".hmm" suffix removed)
annot_cazy <- read_tsv(file.path(annot_dir, "dbcan_annot.tsv"), show_col_types = FALSE) %>%
  transmute(GeneID, GeneMetabolism = "CAZyme", Gene = sub("\\.hmm$", "", cazy_id))


# ==============================================================================
# 4. PHYLOSEQ OBJECTS
# ==============================================================================

# Keep proteins in `annot`, with taxonomy (GeneMetabolism, Gene)
build_filtered_ps <- function(ps, annot) {
  annot_sub <- annot %>%
    filter(GeneID %in% taxa_names(ps)) %>%
    distinct(GeneID, .keep_all = TRUE)

  ps_sub <- prune_taxa(annot_sub$GeneID, ps)

  tax_df <- annot_sub[match(taxa_names(ps_sub), annot_sub$GeneID),
                      c("GeneMetabolism", "Gene")] %>% as.data.frame()
  rownames(tax_df) <- taxa_names(ps_sub)

  tax_table(ps_sub) <- tax_table(as.matrix(tax_df))
  ps_sub
}

ps_ko             <- build_filtered_ps(kos_ps, ko_annot)
ps_greening_hyddb <- build_filtered_ps(kos_ps, annot_greening_hyddb)
ps_dbcan          <- build_filtered_ps(kos_ps, annot_cazy)

sapply(list(ps_ko = ps_ko, ps_greening_hyddb = ps_greening_hyddb, ps_dbcan = ps_dbcan), ntaxa)

saveRDS(ps_ko,             file.path(out_dir, "ps_ko.RDS"))
saveRDS(ps_greening_hyddb, file.path(out_dir, "ps_greening_hyddb.RDS"))
saveRDS(ps_dbcan,          file.path(out_dir, "ps_dbcan.RDS"))


# ==============================================================================
sessionInfo()