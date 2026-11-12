# ==============================================================================
# 5.15_cgc_substrate_heatmap.R
#
# Polysaccharide substrates of CAZyme gene clusters (CGCs) in sediment and
# seawater MAGs, by order (Fig. 4).
#   1. Setup and input
#   2. CGC substrate counts per MAG, summed per habitat and order
#   3. Substrates shown and row order
#   4. Side bars: substrate counts per habitat
#   5. Top bar and column annotation: MAGs per order, habitat, and phylum
#   6. Relative abundance of each order (SingleM)
#   7. Heatmaps
#
# Cells: number of CGCs predicted (run_dbcan easy_substrate, 5.02) to act on
# each substrate, summed over the MAGs of an order (columns, split by habitat).
# Side bars: values per habitat from sheet substrate_count_MAGs_detail. Rows are
# ordered by the habitat with the higher value (sediment first), then by that
# value. Column order: habitat, then phylum (class for Pseudomonadota).
#
# Input : MAG metadata (sheet derep: Source, Phylum, Phylum_Class, Order, Family);
#         substrate_summary.xlsx (sheets substrate_matrix_detail_cgcs: substrate
#         x MAG CGC counts; substrate_count_MAGs_detail: Source (substrate),
#         Sediment, Seawater); coverm_seased_order.xlsx (sheet
#         SeaSedOrder_singlem: relative abundance per order)
# Output: derep85_mag_heatmap.substrate.svg (heatmap) and
#         derep85_mag_heatmap.substrate.combined.svg (with the relative
#         abundance heatmap above)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(ComplexHeatmap)
  library(circlize)
  library(openxlsx)
  library(tidyverse)
})

# merge_samples2() from speedyseq (M. McLaren): remotes::install_github("mikemc/speedyseq")
merge_samples2 <- speedyseq::merge_samples2

mag_metadata_xlsx <- "<path/to/Sea__Sed_indiass_coass.xlsx>"
substrate_xlsx    <- "<path/to/5_MAG_annotation/2_CAZymes_CGCs/substrate_summary.xlsx>"
relabund_xlsx     <- "<path/to/coverm_seased_order.xlsx>"
out_dir           <- "<path/to/5_MAG_annotation/2_CAZymes_CGCs/figures>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_level <- "Source_Phylum_Class_Order"   # columns: habitat x order
bar_digits  <- 1                             # decimals of the side bar labels

substrates_shown <- c(
  "Beta-glucan_degradation", "Alpha-glucan_degradation", "Glycogen_biosynthesis",
  "Starch_degradation", "Xylan_degradation", "Glycogen_degradation", "Pectin_degradation",
  "Fucoidan_degradation", "Alginate_degradation", "Beta-mannan_degradation",
  "Host glycan_degradation", "Mucin_degradation", "N-glycan_degradation",
  "Alpha-galactan_degradation", "Raffinose_degradation", "Alginate_biosynthesis",
  "Starch_biosynthesis", "Alpha-glucan_biosynthesis"
)

source_colors <- c(Sediment = "#33a02c", Seawater = "#1f78b4")
phylum_colors <- c(
  "Asgardarchaeota" = "#A6CEE3", "Thermoplasmatota" = "#B2DF8A", "Thermoproteota" = "#33A02B",
  "4484-113" = "#FB9A99", "Acidobacteriota" = "#E31A1C", "Actinomycetota" = "#FDBF6F",
  "Alphaproteobacteria" = "#FF7F00", "Bacillota" = "#6A3D9A", "Bacteroidota" = "#FFFF99",
  "Bacteroidota_A" = "#B15928", "Bdellovibrionota" = "#1ff8ff", "Bdellovibrionota_B" = "#1B9E77",
  "Bipolaricaulota" = "#D95F02", "Campylobacterota" = "#7570B3", "Chlamydiota" = "#E7298A",
  "Chloroflexota" = "#66A61E", "Cyanobacteriota" = "#E6AB02", "Desulfobacterota" = "#A6761D",
  "Desulfobacterota_D" = "#666666", "Gammaproteobacteria" = "#b249d5", "Gemmatimonadota" = "#4b6a53",
  "Hydrogenedentota" = "#7edc45", "Krumholzibacteriota" = "#5c47b8", "Latescibacterota" = "#cfd251",
  "Margulisbacteria" = "#ff69b4", "Marinisomatota" = "#cd3e50", "Myxococcota" = "#69c86c",
  "Myxococcota_A" = "#da6130", "Nitrospirota" = "#5e79b2", "Patescibacteriota" = "#c29545",
  "Planctomycetota" = "#532a5a", "SAR324" = "#5f7b35", "Spirochaetota" = "#c497cf",
  "Tectomicrobia" = "#773a27", "Verrucomicrobiota" = "#7cb9cb", "Nanobdellota" = "#1F78B1",
  "Babelota" = "#CAB2D6", "Muiribacteriota" = "#83d5af"
)


# ==============================================================================
# 2. CGC SUBSTRATE COUNTS, SUMMED PER HABITAT AND ORDER
# ==============================================================================

mag_meta <- read.xlsx(mag_metadata_xlsx, sheet = "derep", rowNames = TRUE) %>%
  mutate(Source_Phylum_Class       = paste(Source, Phylum_Class, sep = "_"),
         Source_Phylum_Class_Order = paste(Source_Phylum_Class, Order, sep = "_"))

counts_raw <- read.xlsx(substrate_xlsx, sheet = "substrate_matrix_detail_cgcs")
counts <- counts_raw[, -1] %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
rownames(counts) <- counts_raw[[1]]

tax <- matrix(rownames(counts), dimnames = list(rownames(counts), "ID"))

ps_mag <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(mag_meta), tax_table(tax))
cat(sprintf("MAGs with substrate predictions and metadata: %d\n", nsamples(ps_mag)))

ps_order <- merge_samples2(ps_mag, group_level, fun_otu = sum)
ps_order <- prune_samples(sample_sums(ps_order) > 0, ps_order)
ps_order <- prune_taxa(taxa_sums(ps_order) > 0, ps_order)


# ==============================================================================
# 3. SUBSTRATES SHOWN AND ROW ORDER
# ==============================================================================

habitat_counts <- read.xlsx(substrate_xlsx, sheet = "substrate_count_MAGs_detail") %>%
  column_to_rownames("Source")

shown <- intersect(substrates_shown, intersect(taxa_names(ps_order), rownames(habitat_counts)))
missing <- setdiff(substrates_shown, shown)
if (length(missing)) message("Substrates not found and left out: ", paste(missing, collapse = ", "))

# Higher in sediment first, then higher in seawater, then ties; within each,
# by the larger value
habitat_counts <- habitat_counts[shown, ] %>%
  mutate(group    = case_when(Sediment > Seawater ~ 1, Seawater > Sediment ~ 2, TRUE ~ 3),
         sort_key = pmax(Sediment, Seawater)) %>%
  arrange(group, desc(sort_key))
row_ids <- rownames(habitat_counts)

ps_order <- prune_taxa(row_ids, ps_order)

samdf <- as(sample_data(ps_order), "data.frame") %>%
  mutate(Source = as.character(Source), Phylum_Class = as.character(Phylum_Class))
col_order <- order(factor(samdf$Source, levels = names(source_colors)), samdf$Phylum_Class)
samdf     <- samdf[col_order, ]

otu_mat <- as(otu_table(ps_order), "matrix")
if (!taxa_are_rows(ps_order)) otu_mat <- t(otu_mat)
otu_mat <- otu_mat[row_ids, rownames(samdf), drop = FALSE]


# ==============================================================================
# 4. SIDE BARS: SUBSTRATE COUNTS PER HABITAT
# ==============================================================================

bar_sed  <- habitat_counts$Sediment
bar_sea  <- -habitat_counts$Seawater
ylim_bar <- 1.2 * max(abs(c(bar_sed, bar_sea)), na.rm = TRUE)

bar_anno <- rowAnnotation(
  Source_Sum_seawater = anno_barplot(bar_sea, gp = gpar(fill = source_colors[["Seawater"]], col = "black"),
                                     width = unit(2, "cm"), bar_width = 0.8, border = TRUE,
                                     ylim = c(-ylim_bar, 0)),
  Source_Sum_sediment = anno_barplot(bar_sed, gp = gpar(fill = source_colors[["Sediment"]], col = "black"),
                                     width = unit(2, "cm"), bar_width = 0.8, border = TRUE,
                                     ylim = c(0, ylim_bar))
)

# Value labels at the end of each bar (rows are not clustered, so row i is row_ids[i])
label_bars <- function(fontsize) {
  n <- length(row_ids)
  y <- unit(1 - (seq_len(n) - 0.5) / n, "npc")
  sed_frac <- pmin(1, pmax(0, coalesce(bar_sed / ylim_bar, 0)))
  sea_frac <- pmin(1, pmax(0, coalesce(abs(bar_sea) / ylim_bar, 0)))
  decorate_annotation("Source_Sum_sediment", {
    grid.text(formatC(bar_sed, format = "f", digits = bar_digits),
              x = unit(sed_frac + 0.02, "npc"), y = y, just = "left", gp = gpar(fontsize = fontsize))
  })
  decorate_annotation("Source_Sum_seawater", {
    grid.text(formatC(abs(bar_sea), format = "f", digits = bar_digits),
              x = unit(1 - sea_frac - 0.02, "npc"), y = y, just = "right", gp = gpar(fontsize = fontsize))
  })
}


# ==============================================================================
# 5. TOP BAR AND COLUMN ANNOTATION
# ==============================================================================

genome_counts <- table(sample_data(ps_mag)[[group_level]])[colnames(otu_mat)]

top_anno <- HeatmapAnnotation(
  Genome_Count = anno_barplot(as.numeric(genome_counts), gp = gpar(fill = "grey30", col = "black"),
                              bar_width = 0.9, border = TRUE),
  annotation_name_side = "left", annotation_name_rot = 90
)

bottom_anno <- HeatmapAnnotation(
  Source = samdf$Source, Phylum_Class = samdf$Phylum_Class,
  col = list(Source = source_colors, Phylum_Class = phylum_colors),
  show_legend = FALSE
)


# ==============================================================================
# 6. RELATIVE ABUNDANCE OF EACH ORDER (SingleM)
# ==============================================================================
# Orders missing from the relative abundance table are shown as 0

rel_abu <- 100 * as.matrix(read.xlsx(relabund_xlsx, sheet = "SeaSedOrder_singlem", rowNames = TRUE))
missing_orders <- setdiff(colnames(otu_mat), rownames(rel_abu))
if (length(missing_orders)) {
  rel_abu <- rbind(rel_abu, matrix(0, length(missing_orders), ncol(rel_abu),
                                   dimnames = list(missing_orders, colnames(rel_abu))))
}
rel_abu <- t(rel_abu[colnames(otu_mat), , drop = FALSE])

rel_max <- max(rel_abu, na.rm = TRUE)
col_fun_rel <- colorRamp2(c(0, min(rel_abu[rel_abu > 0]), rel_max / 2, rel_max),
                          c("#d3d3d3", "#3b4798", "#faf6c0", "#981f2b"))


# ==============================================================================
# 7. HEATMAPS
# ==============================================================================

max_val <- max(otu_mat, na.rm = TRUE)
col_fun <- colorRamp2(c(0, 0.001, max_val / 2, max_val), c("#d3d3d3", "white", "#40347c", "#be3831"))

ht_rel <- Heatmap(rel_abu, name = "Relative abundance (%)", col = col_fun_rel,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  row_names_side = "left", show_column_names = FALSE,
                  column_split = samdf$Source, border = TRUE, column_gap = unit(0.5, "mm"))

ht <- Heatmap(
  otu_mat, name = "CGCs", col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 0.1, fill = fill))
  },
  top_annotation = top_anno, bottom_annotation = bottom_anno, right_annotation = bar_anno,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  column_split = samdf$Source, border = TRUE,
  row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"),
  use_raster = TRUE, show_heatmap_legend = FALSE,
  width = unit(30, "cm"), height = unit(20, "cm")
)

heatmap_legend <- Legend(title = "CGCs (summed per order)", col_fun = col_fun, direction = "horizontal")
source_legend  <- Legend(title = "Source", at = names(source_colors),
                         legend_gp = gpar(fill = source_colors), ncol = 2)

draw_figure <- function(file, ht_list, fontsize) {
  svg(file.path(out_dir, file), width = 18, height = 15, bg = "white")
  draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       heatmap_legend_list = list(heatmap_legend), annotation_legend_list = list(source_legend))
  label_bars(fontsize)
  dev.off()
}

draw_figure("derep85_mag_heatmap.substrate.svg", ht, fontsize = 10)
draw_figure("derep85_mag_heatmap.substrate.combined.svg", ht_rel %v% ht, fontsize = 8)
