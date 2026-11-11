# ==============================================================================
# 2.14_protein_annotation_objects.R
#
# Phyloseq objects of normalised RPKM (2.13) for KOs, metabolic marker genes
# with HydDB subgroups, and CAZyme families, each with a common taxonomy
# (GeneMetabolism, Gene).
#   1. Setup and input paths
#   2. Normalised RPKM object and updated sample metadata
#   3. Annotation tables
#   4. Phyloseq objects
#
# Marker genes: the generic [NiFe] and [FeFe] hydrogenase marker genes are
# replaced by HydDB subgroups, and HydDB takes priority for proteins hit by both.
#
# Input : dram_ko_e6.RPKM.rds (2.13); marker gene hits with updated gene names
#         (GeneID, Gene, GeneName, GeneMetabolism); HydDB classification
#         (GeneID, hyddb_group); run_dbcan families (GeneID, cazy_id); sample
#         metadata (xlsx, sample IDs as row names)
# Output: ps_greening_hyddb.RDS, ps_ko.RDS, ps_dbcan.RDS
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT PATHS
# ==============================================================================
mem.maxVSize(64000); gc()

suppressPackageStartupMessages({
  library(phyloseq)
  library(openxlsx)
  library(tidyverse)
})

rpkm_rds       <- "<path/to/3_read_mapping_RPKM/RPKM>/dram_ko_e6.RPKM.rds"
metadata_path  <- "<path/to/genecentric_kegg_cazy.xlsx>"
greening_annot <- "<path/to/2_annotation/summary>/1_greeningdb_annot.tsv"
hyddb_annot    <- "<path/to/2_annotation/summary>/3_hyddb_annot.tsv"
dbcan_annot    <- "<path/to/2_annotation/summary>/2_dbcan_annot.tsv"
out_dir        <- "<path/to/3_read_mapping_RPKM/RPKM>"

# ==============================================================================
# 2. NORMALISED RPKM OBJECT AND SAMPLE METADATA
# ==============================================================================

kos_ps   <- readRDS(rpkm_rds)
new_meta <- sample_data(read.xlsx(metadata_path, rowNames = TRUE))

stopifnot(all(sample_names(kos_ps) %in% sample_names(new_meta)))
sample_data(kos_ps) <- new_meta[sample_names(kos_ps), ]


# ==============================================================================
# 3. ANNOTATION TABLES
# ==============================================================================

# KOs, from the tax_table of the RPKM object
ko_annot <- as.data.frame(tax_table(kos_ps)) %>%
  rownames_to_column("GeneID") %>%
  filter(ko_id != "Unassigned") %>%
  transmute(GeneID, GeneMetabolism = "KO", Gene = ko_id)

# Metabolic marker genes, without the generic hydrogenase markers
annot_greening <- read.delim(greening_annot) %>%
  filter(!Gene %in% c("NiFe_hydrogenase", "FeFe_hydrogenase")) %>%
  select(-GeneName)

# HydDB subgroups
annot_hyddb <- read.delim(hyddb_annot) %>%
  transmute(GeneID, GeneMetabolism = "Trace gas metabolism", Gene = hyddb_group)

# Proteins in both tables are taken from HydDB
cat(sprintf("Proteins with both a marker gene and a HydDB hit: %d\n",
            sum(annot_greening$GeneID %in% annot_hyddb$GeneID)))
annot_greening_hyddb <- bind_rows(
  filter(annot_greening, !GeneID %in% annot_hyddb$GeneID),
  annot_hyddb
)

# CAZyme families
annot_cazy <- read.delim(dbcan_annot) %>%
  transmute(GeneID, GeneMetabolism = "CAZyme", Gene = sub("\\.hmm$", "", cazy_id))


# ==============================================================================
# 4. PHYLOSEQ OBJECTS
# ==============================================================================
# Keep only proteins in `annot`, with tax_table (GeneMetabolism, Gene)

build_filtered_ps <- function(ps, annot) {
  annot_sub <- annot %>%
    filter(GeneID %in% taxa_names(ps)) %>%
    distinct(GeneID, .keep_all = TRUE)
  ps_sub <- prune_taxa(annot_sub$GeneID, ps)
  tax_df <- annot_sub[match(taxa_names(ps_sub), annot_sub$GeneID), c("GeneMetabolism", "Gene")] %>%
    as.data.frame()
  rownames(tax_df) <- taxa_names(ps_sub)
  tax_table(ps_sub) <- tax_table(as.matrix(tax_df))
  ps_sub
}

ps_ko             <- build_filtered_ps(kos_ps, ko_annot)
ps_greening_hyddb <- build_filtered_ps(kos_ps, annot_greening_hyddb)
ps_cazy           <- build_filtered_ps(kos_ps, annot_cazy)

print(sapply(list(ps_ko = ps_ko, ps_greening_hyddb = ps_greening_hyddb, ps_cazy = ps_cazy), ntaxa))

saveRDS(ps_greening_hyddb, file.path(out_dir, "ps_greening_hyddb.RDS"))
saveRDS(ps_ko,             file.path(out_dir, "ps_ko.RDS"))
saveRDS(ps_cazy,           file.path(out_dir, "ps_dbcan.RDS"))
