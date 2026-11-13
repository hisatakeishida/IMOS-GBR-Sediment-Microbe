# ==============================================================================
# 2.13_protein_RPKM.R
#
# Normalised abundance of catalogue proteins annotated with DRAM KOfam
# (E-value 1e-6): RPKM divided by the mean RPKM of SingleM marker genes.
#   1. Setup and input paths
#   2. Protein counts: low-abundance filter and count matrix
#   3. Sample metadata and KO annotation
#   4. Count phyloseq object and library size
#   5. RPKM for catalogue proteins
#   6. RPKM for SingleM marker genes
#   7. Normalise protein RPKM by mean SingleM marker gene RPKM
#
# Low-abundance filter: protein counts < 1e-6 (0.0001%) of a sample's mapped
# reads are removed.
# RPKM = count / (mapped reads in sample / 1e6) / (gene length in nt / 1e3)
# Normalised RPKM = protein RPKM / mean RPKM of all SingleM marker genes in the sample
#
# Input : long-format read counts (Sample, GeneID, Count) for catalogue proteins
#         and SingleM marker genes (2.12); gene lengths (2.10); KO annotation
#         (GeneID, ko_id; 2.08); sample metadata (xlsx, sample IDs as row names)
# Output: dram_ko_e6.COUNT.rds and dram_ko_e6.RPKM.rds (phyloseq objects,
#         proteins as taxa); numreadsmapped.png (mapped reads per sample)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT PATHS
# ==============================================================================

# On the HPC, packages were installed in a project library:
# .libPaths("<path/to/R/Rlibs>")

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(data.table)
  library(openxlsx)
  library(tidyverse)
})

mode <- "dram_ko_e6"

countmat_path       <- "<path/to/3_read_mapping_RPKM/summary>/gene_catalogue_counts.longformat.tsv"
singlem_otumat_path <- "<path/to/3_read_mapping_RPKM/summary>/singlem_counts.longformat.tsv"
protein_len         <- "<path/to/3_read_mapping_RPKM/summary>/gene_catalogue_protein_length.tsv"
singlem_len         <- "<path/to/3_read_mapping_RPKM/summary>/singlem_protein_length.tsv"

genename_path       <- "<path/to/2_annotation/summary>/gene_catalogue_ko_annotations.tsv"
metadata_path       <- "<path/to/genecentric_kegg_cazy.xlsx>"
out_dir             <- "<path/to/3_read_mapping_RPKM/RPKM>"
fig_dir             <- file.path(out_dir, "FIG", mode)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

count_rds <- file.path(out_dir, paste0(mode, ".COUNT.rds"))
rpkm_rds  <- file.path(out_dir, paste0(mode, ".RPKM.rds"))
min_rel_abundance <- 1e-6

# Long-format counts to a numeric matrix (genes x samples)
long_to_matrix <- function(long) {
  wide <- dcast(as.data.table(long), GeneID ~ Sample, value.var = "Count", fill = 0)
  mat  <- as.data.frame(wide[, -1])
  rownames(mat) <- wide$GeneID
  mat[] <- lapply(mat, as.numeric)
  mat
}

# RPKM of a genes x samples count matrix, given lengths in nucleotides
rpkm <- function(counts, len_path) {
  len <- fread(len_path, header = FALSE, sep = "\t", col.names = c("GeneID", "Length"))
  len <- len[match(rownames(counts), len$GeneID), ]
  if (anyNA(len$GeneID) || !all(rownames(counts) == len$GeneID)) {
    stop("Gene IDs in ", basename(len_path), " do not match the count matrix.")
  }
  rpm <- sweep(counts, 2, colSums(counts) / 1e6, FUN = "/")
  sweep(rpm, 1, len$Length / 1e3, FUN = "/")
}


# ==============================================================================
# 2. PROTEIN COUNTS
# ==============================================================================

otu_long <- fread(countmat_path, header = TRUE, sep = "\t")

otu_long <- otu_long %>%
  group_by(Sample) %>%
  mutate(rel_abundance = Count / sum(Count)) %>%
  ungroup() %>%
  filter(rel_abundance >= min_rel_abundance) %>%
  select(Sample, GeneID, Count)

otu <- long_to_matrix(otu_long)
rm(otu_long); gc()
cat(sprintf("Proteins: %d; samples: %d\n", nrow(otu), ncol(otu)))


# ==============================================================================
# 3. SAMPLE METADATA AND KO ANNOTATION
# ==============================================================================

sample_data_in <- sample_data(read.xlsx(metadata_path, rowNames = TRUE))

tax_rel <- as.matrix(fread(genename_path, header = TRUE, sep = "\t"))
rownames(tax_rel) <- tax_rel[, 1]
tax_rel <- tax_rel[rownames(tax_rel) %in% rownames(otu), -1, drop = FALSE]
tax_rel[is.na(tax_rel) | tax_rel == ""] <- "Unassigned"
tax_in <- tax_table(tax_rel)


# ==============================================================================
# 4. COUNT PHYLOSEQ OBJECT AND LIBRARY SIZE
# ==============================================================================

otu_annot <- otu[rownames(otu) %in% rownames(tax_rel), ]

pseq <- phyloseq(otu_table(otu_annot, taxa_are_rows = TRUE), sample_data_in, tax_in) %>%
  tax_fix(min_length = 1) %>%
  phyloseq_validate()
saveRDS(pseq, count_rds)

lib_df <- tibble(Sample = sample_names(pseq),
                 LibrarySize = sample_sums(pseq),
                 Source = sample_data(pseq)$Source) %>%
  mutate(Sample = fct_reorder(Sample, LibrarySize))

p_lib <- ggplot(lib_df, aes(Sample, LibrarySize, fill = Source)) +
  geom_col() +
  scale_fill_manual(values = c(Seawater = "#1f78b4", Sediment = "#33a02c")) +
  labs(title = "Mapped reads per sample (KO-annotated proteins)", x = "Sample", y = "Reads") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
ggsave(file.path(fig_dir, "numreadsmapped.png"), p_lib, width = 7.5, height = 7.5)


# ==============================================================================
# 5. RPKM FOR CATALOGUE PROTEINS
# ==============================================================================
# Per-million scaling uses all proteins passing the filter, not only the
# KO-annotated ones

rpkm_prot <- rpkm(otu, protein_len)


# ==============================================================================
# 6. RPKM FOR SingleM MARKER GENES
# ==============================================================================

singlem_otu <- long_to_matrix(fread(singlem_otumat_path, header = TRUE, sep = "\t"))
singlem_otu <- singlem_otu[, colnames(singlem_otu) %in% colnames(otu), drop = FALSE]

rpkm_singlem      <- rpkm(singlem_otu, singlem_len)
singlem_mean_rpkm <- colMeans(rpkm_singlem, na.rm = TRUE)
print(summary(singlem_mean_rpkm))


# ==============================================================================
# 7. NORMALISED RPKM
# ==============================================================================

stopifnot(all(colnames(rpkm_prot) %in% names(singlem_mean_rpkm)))
rpkm_prot_scaled <- sweep(rpkm_prot, 2, singlem_mean_rpkm[colnames(rpkm_prot)], FUN = "/")

pseq_rpkm <- phyloseq(otu_table(rpkm_prot_scaled, taxa_are_rows = TRUE), sample_data_in, tax_in) %>%
  tax_fix(min_length = 1) %>%
  phyloseq_validate()
pseq_rpkm

saveRDS(pseq_rpkm, rpkm_rds)
