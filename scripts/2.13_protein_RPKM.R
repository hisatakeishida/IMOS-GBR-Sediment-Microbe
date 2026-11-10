# ==============================================================================
# 2.13_protein_RPKM.R
#
# Normalised protein abundance (RPKM scaled by SingleM marker genes) for one
# annotation (KO, CAZyme, HydDB, or marker genes).
#   1. Setup
#   2. Protein counts: low-abundance filter and count matrix
#   3. Annotation and sample metadata
#   4. Count phyloseq object and library size
#   5. RPKM for catalogue proteins
#   6. RPKM for SingleM marker genes
#   7. Normalise protein RPKM by mean SingleM marker gene RPKM
#
# RPKM = count / (mapped reads in sample / 1e6) / (gene length / 1e3)
# Normalised RPKM = protein RPKM / mean RPKM of all SingleM marker genes
#
# Usage: Rscript scripts/2.13_protein_RPKM.R ko      # or cazy, hyddb, marker
#
# Input : long-format counts and lengths (2.10, 2.12), annotation tables (2.8),
#         sample metadata
# Output: <annotation>.COUNT.rds and <annotation>.RPKM.rds (phyloseq objects,
#         proteins as taxa), library size plot
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(data.table)
  library(openxlsx)
  library(tidyverse)
})

args       <- commandArgs(trailingOnly = TRUE)
annotation <- if (length(args) > 0) args[1] else "ko"

in_dir        <- "<path/to/3_read_mapping_RPKM>"
annot_dir     <- "<path/to/2_annotation/summary>"
metadata_path <- "<path/to/genecentric_kegg_cazy.xlsx>"
out_dir       <- file.path("<path/to/3_read_mapping_RPKM/RPKM>", annotation)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

count_path         <- file.path(in_dir, "diamond_long.catalogue.tsv")
singlem_count_path <- file.path(in_dir, "diamond_long.singlem.tsv")
protein_len_path   <- file.path(in_dir, "comp_protein_length.fixed.tsv")
singlem_len_path   <- file.path(in_dir, "singlem_marker_length.tsv")

annot_files <- c(ko     = "ko_annot.tsv",
                 cazy   = "dbcan_annot.tsv",
                 hyddb  = "hyddb_annot.tsv",
                 marker = "marker_genes_annot.tsv")
if (!annotation %in% names(annot_files)) stop("Unknown annotation: ", annotation)
annot_path <- file.path(annot_dir, annot_files[[annotation]])

min_relabu   <- 1e-6   # protein counts < 0.0001% of a sample's mapped reads are removed
habitat_cols <- c(Seawater = "#1f78b4", Sediment = "#33a02c")

# Long format (Sample, GeneID, Count) to genes x samples matrix
long_to_matrix <- function(long) {
  wide <- dcast(as.data.table(long), GeneID ~ Sample, value.var = "Count", fill = 0)
  mat  <- as.matrix(wide[, -1])
  rownames(mat) <- wide$GeneID
  mat
}

# RPKM from a genes x samples count matrix and a two-column length file (bp)
rpkm <- function(counts, len_path) {
  len    <- fread(len_path, header = FALSE, col.names = c("GeneID", "Length"))
  len_bp <- len$Length[match(rownames(counts), len$GeneID)]
  if (anyNA(len_bp)) stop("Gene IDs missing from ", len_path)
  rpm <- sweep(counts, 2, colSums(counts) / 1e6, FUN = "/")
  sweep(rpm, 1, len_bp / 1e3, FUN = "/")
}


# ==============================================================================
# 2. PROTEIN COUNTS
# ==============================================================================

counts_long <- fread(count_path) %>%
  group_by(Sample) %>%
  mutate(rel_abundance = Count / sum(Count)) %>%
  ungroup() %>%
  filter(rel_abundance >= min_relabu) %>%
  select(Sample, GeneID, Count)

counts <- long_to_matrix(counts_long)
rm(counts_long)
message(nrow(counts), " proteins x ", ncol(counts), " samples after filtering")


# ==============================================================================
# 3. ANNOTATION AND SAMPLE METADATA
# ==============================================================================

tax <- fread(annot_path, header = TRUE, colClasses = "character")
tax <- tax[tax[[1]] %in% rownames(counts)]
tax_mat <- as.matrix(tax[, -1, with = FALSE])
rownames(tax_mat) <- tax[[1]]
tax_mat[is.na(tax_mat) | tax_mat == ""] <- "Unassigned"
message(nrow(tax_mat), " mapped proteins with ", annotation, " annotation table entries")

meta <- sample_data(read.xlsx(metadata_path, rowNames = TRUE))

# Proteins in the annotation table only; RPKM (sections 5 to 7) uses all mapped
# proteins so that sequencing depth is calculated from all mapped reads
make_ps <- function(mat) {
  mat <- mat[rownames(mat) %in% rownames(tax_mat), , drop = FALSE]
  phyloseq(otu_table(mat, taxa_are_rows = TRUE), meta, tax_table(tax_mat)) %>%
    tax_fix(min_length = 1) %>%
    phyloseq_validate()
}


# ==============================================================================
# 4. COUNT PHYLOSEQ OBJECT AND LIBRARY SIZE
# ==============================================================================

ps_count <- make_ps(counts)
ps_count
saveRDS(ps_count, file.path(out_dir, paste0(annotation, ".COUNT.rds")))

lib_df <- tibble(Sample      = sample_names(ps_count),
                 LibrarySize = sample_sums(ps_count),
                 Source      = sample_data(ps_count)$Source) %>%
  mutate(Sample = fct_reorder(Sample, LibrarySize))

p_lib <- ggplot(lib_df, aes(x = Sample, y = LibrarySize, fill = Source)) +
  geom_col() +
  scale_fill_manual(values = habitat_cols) +
  labs(title = "Mapped reads per sample", x = "Sample", y = "Reads") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top")

ggsave(file.path(out_dir, "numreadsmapped.png"), p_lib, width = 7.5, height = 7.5)


# ==============================================================================
# 5. RPKM FOR CATALOGUE PROTEINS
# ==============================================================================

rpkm_prot <- rpkm(counts, protein_len_path)


# ==============================================================================
# 6. RPKM FOR SINGLEM MARKER GENES
# ==============================================================================

singlem_counts <- long_to_matrix(fread(singlem_count_path))
singlem_counts <- singlem_counts[, colnames(singlem_counts) %in% colnames(counts), drop = FALSE]

rpkm_singlem      <- rpkm(singlem_counts, singlem_len_path)
singlem_mean_rpkm <- colMeans(rpkm_singlem, na.rm = TRUE)
summary(singlem_mean_rpkm)


# ==============================================================================
# 7. NORMALISE PROTEIN RPKM BY MEAN SINGLEM MARKER GENE RPKM
# ==============================================================================

stopifnot(all(colnames(rpkm_prot) %in% names(singlem_mean_rpkm)))
rpkm_norm <- sweep(rpkm_prot, 2, singlem_mean_rpkm[colnames(rpkm_prot)], FUN = "/")

ps_rpkm <- make_ps(rpkm_norm)
ps_rpkm
saveRDS(ps_rpkm, file.path(out_dir, paste0(annotation, ".RPKM.rds")))


# ==============================================================================
sessionInfo()
