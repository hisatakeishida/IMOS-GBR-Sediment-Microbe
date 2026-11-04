# **********************************************************************************
# Loading packages 
# **********************************************************************************

library("data.table")
library("phyloseq")
library("microbiome")
library("ComplexHeatmap")
library("microViz")
library("ggtext")
library("ggraph")
library("DT")
library("corncob")
library("openxlsx")
library("dplyr")
library("ggrepel")

# **********************************************************************************
# Loading datasets 
# **********************************************************************************
#datasets
countmat_path="ko_diamond_long_format.tsv"
metadata_path="genecentric_kegg_cazy.xlsx"
genename_path="0_genename-KO.tsv"
protein_len="1_genelength-KO.tsv"
singlem_otumat_path="singlem_countmatrix.tsv"
singlem_len="2_genelength-SingleM.tsv"
figout="kos"
RDSOUT="KO.RPKM.rds"

#protein information
# Prepare OTU table
otu_raw <- fread(countmat_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Get relative abundance for filtering
otu_raw <- otu_raw %>%
  group_by(Sample) %>%
  mutate(total_count = sum(Count),
         rel_abundance = Count / total_count) %>%
  ungroup()

# Actual filtering
otu_raw <- otu_raw %>%
  filter(rel_abundance >= 1e-6) %>%
  select(Sample, GeneID, Count)

print("Long to Matrix")
otu_raw <- as.data.table(otu_raw)
otu_mat <- dcast(otu_raw, GeneID ~ Sample, value.var = "Count", fill = 0)
otu <- as.data.frame(otu_mat[, -1])
rownames(otu) <- otu_mat$GeneID
otu[] <- lapply(otu, as.numeric)

#sample information
sample_data <- sample_data(read.xlsx(metadata_path, sheet = "kegg_env", rowNames=TRUE))
rownames(sample_data) <- gsub("\\.kegg$", "", rownames(sample_data))
tax_rel <- as.matrix(fread(genename_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
rownames(tax_rel) <- tax_rel[,1]
tax_rel <- tax_rel[, -1, drop = FALSE]
tax_rel_filtered <- tax_rel[rownames(tax_rel) %in% rownames(otu), , drop = FALSE]
otu2 <- otu[rownames(otu) %in% rownames(tax_rel_filtered), ]

#into phyloseq object
sample_data_in <- sample_data(sample_data)
otu_in <- otu_table(otu2, taxa_are_rows = TRUE)
tax_in <- tax_table(tax_rel_filtered)
pseq = phyloseq(otu_table = otu_in, sample_data=sample_data_in, taxonomyTable =tax_in)
pseq

pseq <- tax_fix(pseq, min_length = 1)
pseq <- phyloseq_validate(pseq)
pseq@sam_data$SampleID <- rownames(sample_data(pseq))
pseq

# **********************************************************************************
# Analsysis based on count (library size and rarefraction curve)
# **********************************************************************************
# looking at library size 
lib.size <- apply(t(otu_table(pseq)), 1, sum)
plot_df <- data.frame(
  Sample = names(lib.size),
  LibrarySize = lib.size,
  Source = sample_data(pseq)$Source
)
plot_df$Sample <- factor(plot_df$Sample, levels = plot_df$Sample[order(plot_df$LibrarySize)])

p <- ggplot(plot_df, aes(x = Sample, y = LibrarySize, fill = Source)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(
    title = "Library Size per Sample",
    x = "Sample",
    y = "Total Reads"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),  # Hide sample labels if too many
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
# jac_path1 <- file.path(figout, "numreadsmapped.png")
# ggsave(file=jac_path1, plot=p, width=7.5, height=7.5)

# looking at rarefraction curve 
mat <- t(otu_table(pseq))
class(mat) <- "matrix"
class(mat)
dim(mat)
raremax <- min(rowSums(mat))
meta_rarecur <-  sample_data(pseq)
dim(meta_rarecur)

meta_rarecur$Site <- rownames(meta_rarecur)
x <- vegan::rarecurve(mat, step = 100000, sample = raremax, col = "blue", label = FALSE,  tidy = TRUE)  %>% left_join(meta_rarecur)

p <- ggplot(x) +
  geom_line(aes(x = Sample, y = Species, group = Site, colour = Source)) +
  theme_bw() +
  scale_color_manual(values = c("#1f78b4", "#33a02c"))

# jac_path1 <- file.path(figout, "rarecurve.ko.100000.png")
# ggsave(file=jac_path1, plot=p, width=7.5, height=7.5)

# **********************************************************************************
# Analsysis based on count (ordination)
# **********************************************************************************

agg_levels <- c("ko_id", "unique")
ordmethods <- c("PCA")
tranmodes <- c("identity", "compositional", "hellinger", "clr")

for (transform in tranmodes) {
  for (agg in agg_levels) {
    for (ord in ordmethods) {
      message("Running: agg = ", agg, ", ord = ", ord, ", transformation = ", transform)
      p <- pseq %>%
        tax_transform(trans = transform, rank = agg) %>%
        # dist_calc("bray") %>%
        ord_calc(method = ord) %>%
        ord_plot(
          size = 4,
          fill = "Source", shape = "TripDescriotion",
          plot_taxa = 10:1,
          tax_vec_length = 0.3,
          tax_vec_style_all = vec_constraint(linewidth = 0, arrow = NULL),
          tax_lab_style = tax_lab_style(
            type = "label", max_angle = 90, aspect_ratio = 1,
            size = 3, alpha = 0.9, fontface = "bold",
            label.r = unit(0, "mm")
          )
        ) +
        # geom_text_repel(aes(label = SampleID), vjust = -1, hjust = 0.5, size = 3, colour = "black") +
        stat_ellipse(aes(linetype = Source, colour = Source), linewidth = 0.2) + theme_bw() +
        theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
        scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
        scale_shape_manual(values = c(21, 22, 23, 24)) +
        scale_color_manual(values = c("Seawater" = "#1f78b4", "Sediment" = "#33a02c"))
      # filename <- paste0("ko_count-", agg, "_", ord, "_", transform, ".png")
      # image_path <- file.path(figout, filename)
      # ggsave(file = image_path, plot = p, width = 7.5, height = 7.5)
    }
  }
}

# **********************************************************************************
# Analsysis based on RPKM 
# **********************************************************************************

# **********************************************************************************
# #RPKM for protein db 
# **********************************************************************************
protein_len_df <- fread(protein_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(protein_len_df) <- c("GeneID", "Length")
protein_len_df <- protein_len_df[match(rownames(otu), protein_len_df$GeneID), ]

## Count up the total mapped reads in a sample and divide that number by 1,000,000. 
## This is our “per million” scaling factor
per_million_scaling <- colSums(otu) / 1e6 

## Divide the read counts by the “per million” scaling factor. 
## This normalizes for sequencing depth, giving you reads per million (RPM)
rpm <- sweep(otu, 2, per_million_scaling, FUN = "/")

## Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
gene_length_kb <- protein_len_df$Length / 1000

## Sanity check: ensure same order between rownames(otu) and gene IDs
if (!all(rownames(otu) == protein_len_df$GeneID)) {
  stop("Gene IDs in gene length file and OTU matrix do not match in order.")
}

## RPKM = RPM / gene length in kb
rpkm_prot <- sweep(rpm, 1, gene_length_kb, FUN = "/")  # divide each row by gene length

# **********************************************************************************
# #RPKM for singlem db 
# **********************************************************************************
singlem_otu_raw <- fread(singlem_otumat_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
singlem_otu_mat <- as.data.frame(singlem_otu_raw)
rownames(singlem_otu_mat) <- singlem_otu_mat[[1]]   # use Sample column as row names
singlem_otu_mat[[1]] <- NULL                # remove the Sample column
singlem_otu <- t(singlem_otu_mat)
singlem_otu <- as.data.frame(singlem_otu)
singlem_otu[] <- lapply(singlem_otu, as.numeric)
colnames(singlem_otu) <- gsub("\\.singlem$", "", colnames(singlem_otu))
singlem_otu_filtered <- singlem_otu[colnames(singlem_otu) %in% colnames(otu), , drop = FALSE]

singlem_len_df <- fread(singlem_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(singlem_len_df) <- c("GeneID", "Length")
singlem_len_df <- singlem_len_df[match(rownames(singlem_otu_filtered), singlem_len_df$GeneID), ]

## Count up the total mapped reads in a sample and divide that number by 1,000,000. 
## This is our “per million” scaling factor
per_million_scaling <- colSums(singlem_otu_filtered) / 1e6 

## Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
rpm <- sweep(singlem_otu_filtered, 2, per_million_scaling, FUN = "/")

## Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
gene_length_kb <- singlem_len_df$Length / 1000

## Sanity check: ensure same order between rownames(otu) and gene IDs
if (!all(rownames(singlem_otu_filtered) == singlem_len_df$GeneID)) {
  stop("Gene IDs in gene length file and OTU matrix do not match in order.")
}

# RPKM = RPM / gene length in kb
rpkm_singlem <- sweep(rpm, 1, gene_length_kb, FUN = "/")  # divide each row by gene length
dim(rpkm_singlem)
head(rpkm_singlem[1:3, 1:3])
rpkm_singlem_mean_persamp <- colMeans(rpkm_singlem, na.rm = TRUE)

# **********************************************************************************
# #Normalize protein RPKM with Singlem RPKM
# **********************************************************************************
# stopifnot(all(colnames(rpkm_prot) %in% names(rpkm_singlem_mean_persamp)))

# Match the means to the order of rpkm_prot columns
scaling_factors <- rpkm_singlem_mean_persamp[colnames(rpkm_prot)]

# Divide each column in rpkm_prot by the corresponding scaling factor
rpkm_prot_scaled <- sweep(rpkm_prot, 2, scaling_factors, FUN = "/")

rpkm_prot_scaled_int <- rpkm_prot_scaled %>% mutate_all(~ as.integer(round(.)))
dim(rpkm_prot_scaled_int)

#into phyloseq object
otu_in_rpkm <- otu_table(rpkm_prot_scaled_int, taxa_are_rows = TRUE)
pseq_rpkm = phyloseq(otu_table = otu_in_rpkm, sample_data=sample_data_in, taxonomyTable =tax_in)
pseq_rpkm <- tax_fix(pseq_rpkm, min_length = 1)
pseq_rpkm <- phyloseq_validate(pseq_rpkm)
pseq_rpkm

#"SAVING AS RDS for later 
saveRDS(pseq_rpkm, RDSOUT)

# # **********************************************************************************
# # Analsysis based on RPKM 
# # **********************************************************************************
# pseq_rpkm <- readRDS(RDSOUT)

lib.size <- apply(t(otu_table(pseq_rpkm)), 1, sum)

# Build data frame with sample IDs and Source
plot_df <- data.frame(
  Sample = names(lib.size),
  LibrarySize = lib.size,
  Source = sample_data(pseq_rpkm)$Source
)
plot_df$Sample <- factor(plot_df$Sample, levels = plot_df$Sample[order(plot_df$LibrarySize)])

p <- ggplot(plot_df, aes(x = Sample, y = LibrarySize, fill = Source)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(
    title = "Library Size per Sample",
    x = "Sample",
    y = "Total Reads"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),  # Hide sample labels if too many
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
# jac_path1 <- file.path(figout, "rpkm.numreadsmapped.png")
# ggsave(file=jac_path1, plot=p, width=7.5, height=7.5)

mat <- t(otu_table(pseq_rpkm))
class(mat) <- "matrix"
class(mat)
dim(mat)
raremax <- min(rowSums(mat))
meta_rarecur <-  sample_data(pseq_rpkm)
meta_rarecur$Site <- rownames(meta_rarecur)
x <- vegan::rarecurve(mat, step = 10000, sample = raremax, col = "blue", label = FALSE,  tidy = TRUE)  %>% left_join(meta_rarecur)

p <- ggplot(x) +
  geom_line(aes(x = Sample, y = Species, group = Site, colour = Source)) +
  theme_bw() +
  scale_color_manual(values = c("#1f78b4", "#33a02c"))

# jac_path1 <- file.path(figout, "rpkm.rarecurve.ko.10000.png")
# ggsave(file=jac_path1, plot=p, width=7.5, height=7.5)

# # **********************************************************************************
# # Distance based ordination
# # **********************************************************************************

agg_levels <- c("ko_id", "unique")
ordmethods <- c("NMDS")
pseq_rpkm@sam_data$SampleID <- rownames(sample_data(pseq_rpkm))

for (agg in agg_levels) {
  for (ord in ordmethods) {
    message("Running: agg = ", agg, ", ord = ", ord)
    p <- pseq_rpkm %>%
      ps_filter(SampleID != "lynchs_1__21S000275_R1_allRunsCombined")  %>%
      tax_agg(agg) %>%
      tax_transform("identity") %>%
      dist_calc("aitchison") %>%
      ord_calc(method = ord) %>%
      ord_plot(
        size = 4,
        fill = "Source", shape = "TripDescriotion",  # fix if typo
        plot_taxa = 10:1,
        tax_vec_length = 0.3,
        tax_vec_style_all = vec_constraint(linewidth = 0, arrow = NULL),
        tax_lab_style = tax_lab_style(
          type = "label", max_angle = 90, aspect_ratio = 1,
          size = 3, alpha = 0.9, fontface = "bold",
          label.r = unit(0, "mm")
        )
      ) +
      geom_text_repel(aes(label = SampleID), vjust = -1, hjust = 0.5, size = 3, colour = "black") +
      stat_ellipse(aes(linetype = Source, colour = Source), linewidth = 0.2) + theme_bw() +
      theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
      scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
      scale_shape_manual(values = c(21, 22, 23, 24)) +
      scale_color_manual(values = c("Seawater" = "#1f78b4", "Sediment" = "#33a02c"))
    # filename <- paste0("ko_rpkm_withname-removedlynch", "aitchison_", agg, "_", ord, ".png")
    # image_path <- file.path(figout, filename)
    # ggsave(file = image_path, plot = p, width = 7.5, height = 7.5)
  }
}
