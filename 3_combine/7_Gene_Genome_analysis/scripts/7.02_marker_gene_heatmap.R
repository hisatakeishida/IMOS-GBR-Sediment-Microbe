# ==============================================================================
# 7.02_marker_gene_heatmap.R
#
# Metabolic marker genes in sediment and seawater microbes from the
# genome-centric and gene-centric approaches (Fig. 3a-c).
#   1. Setup and input
#   2. Marker genes and HydDB subgroups grouped by metabolism
#   3. Genome-centric: % of MAGs per phylum (class for Pseudomonadota) encoding
#      each function (Fig. 3a)
#   4. Genome-centric: mean copy number and Mann-Whitney U log2 fold change (Fig. 3b)
#   5. Gene-centric: mean relative abundance and MaAsLin2/LinDA consensus (Fig. 3c)
#   6. Heatmap
#
# A MAG encodes a function (row) if it has at least one hit to any of the
# function's genes. Genome-centric values are for the 796 MAGs with CheckM2
# completeness >= 85% used in 5.08. For functions with several genes, side
# panels show: summed mean copy number; mean log2 fold change and lowest
# adjusted p (Mann-Whitney U); mean relative abundance; and significant if any
# gene is supported by both MaAsLin2 and LinDA (adjusted p <= 0.05, same sign).
# Side bar lengths are on a square-root scale.
#
# Input : marker_genes_countmat_UsedInEnrichment.csv, marker_genes_maginfo_UsedInEnrichment.csv,
#         marker_genes_mannwhitney.csv (5.08); MAASLIN_LINDA_greeningDB.xlsx (2.16)
# Output: Fig. 3a-c heatmap (svg, with and without legends); table of plotted values (xlsx)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(svglite)
  library(openxlsx)
  library(tidyverse)
})
select <- dplyr::select

enrich_dir <- "<path/to/5_MAG_annotation/6_enrichment_tests>"     # 5.08 output
da_file    <- "<path/to/2_Protein/DA/MAASLIN_LINDA_greeningDB.xlsx>" # 2.16 output
out_dir    <- "<path/to/7_Gene_Genome_analysis/2_metabolism>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q_cutoff <- 0.05

# Gene names pass through make.names() so they match across tables
# (e.g. "[NiFe]_1a" becomes "X.NiFe._1a")
genome_genetable <- read.csv(file.path(enrich_dir, "countmat_UsedInEnrichment.csv"),
                             row.names = 1, check.names = TRUE)
genome_metadata  <- read.csv(file.path(enrich_dir, "maginfo_UsedInEnrichment.csv"))
manw_result      <- read.csv(file.path(enrich_dir, "manw_enrichment.csv")) %>%
  mutate(gene = make.names(gene))
gene_enrichment  <- read.xlsx(da_file) %>%
  rename(gene = feature)

stopifnot(all(rownames(genome_genetable) %in% genome_metadata$MAG))  # MAGs as rows


# ==============================================================================
# 2. FUNCTIONS (ROWS) AND THEIR GENES
# ==============================================================================

metab_groups <- list(
  "Aerobic respiration" = list(
    "Oxidative phosphorylation (atpA)"        = "atpA",
    "NADH oxidation (nuoF)"                   = "nuoF",
    "Succinate oxidation (sdhA)"              = "sdhA",
    "Aerobic respiration High O2 (coxA,cyoA)" = c("coxA", "cyoA"),
    "Aerobic respiration Low O2 (cydA,ccoN)"  = c("cydA", "ccoN")
  ),
  "Alternative electron acceptor" = list(
    "Fumarate reduction (frdA)"       = "frdA",
    "Reductive dehalogenation (rdhA)" = "rdhA",
    "Iron reduction (mtrB,omcB)"      = c("mtrB", "omcB"),
    "Selenate reduction (ygfK)"       = "ygfK"
  ),
  "Alternative electron donor" = list(
    "Formate oxidation (fdhA)" = "fdhA",
    "Arsenite oxidation (ARO)" = "ARO",
    "Iron oxidation (cyc2)"    = "cyc2"
  ),
  "Nitrogen metabolism" = list(
    "Nitrate reduction (narG,napA)"                  = c("narG", "napA"),  # NO3- to NO2-
    "Denitrification:Nitrite reduction (nirK,nirS)"  = c("nirK", "nirS"),  # NO2- to NO
    "Denitrification:Nitric oxide reduction (norB)"  = "norB",             # NO to N2O
    "Denitrification:Nitrous oxide reduction (nosZ)" = "nosZ",             # N2O to N2
    "DNRA (nrfA)"                                    = "nrfA",             # NO2- to NH4+
    "anammox (hzsA)"                                 = "hzsA",
    "Nitrogen fixation (nifH)"                       = "nifH",
    "Nitrification:Ammonia oxidation (amoA)"         = "amoA",
    "Nitrification:Nitrite oxidation (nxrA)"         = "nxrA"              # NO2- to NO3-
  ),
  "Sulfur metabolism" = list(
    "Sulfate reduction (dsrA)"        = "dsrA",
    "Sulfide oxidation (sqr,FCC)"     = c("sqr", "FCC"),
    "Sulfur disproportionation (sor)" = "sor",
    "Thiosulfate oxidation (soxB)"    = "soxB"
  ),
  "Trace gas metabolism" = list(
    "Aerobic CO oxidation (coxL)"   = "coxL",
    "Anaerobic CO oxidation (cooS)" = "cooS",
    "[NiFe] G1 (H2 uptake)"         = c("X.NiFe._1a", "X.NiFe._1b", "X.NiFe._1c", "X.NiFe._1d", "X.NiFe._1e",
                                        "X.NiFe._1f", "X.NiFe._1h", "X.NiFe._1k", "X.NiFe._1l"),
    "[NiFe] G2 (Sensing/H2 uptake)" = c("X.NiFe._2a", "X.NiFe._2b", "X.NiFe._2c", "X.NiFe._2d"),
    "[NiFe] G3 (H2 fermentation)"   = c("X.NiFe._3a", "X.NiFe._3b", "X.NiFe._3c", "X.NiFe._3d"),
    "[NiFe] G4 (H2 evolving)"       = c("X.NiFe._4a", "X.NiFe._4b", "X.NiFe._4c", "X.NiFe._4d",
                                        "X.NiFe._4e", "X.NiFe._4f", "X.NiFe._4g"),
    "[FeFe] Group A (Fermentative/Bifurcating)" = c("X.FeFe._A1", "X.FeFe._A2", "X.FeFe._A3"),
    "[FeFe] Group B (Colonic)"                  = "X.FeFe._B",
    "[FeFe] Group C (Sensing)"                  = c("X.FeFe._C1", "X.FeFe._C3")
  ),
  "Carbon metabolism" = list(
    "rTCA cycle (aclB)"              = "aclB",
    "Wood-Ljungdahl pathway (acsB)"  = "acsB",
    "4-hydroxybutyrate cycle (hbsT)" = "hbsT",
    "3-hydroxypropionate cycle (mcr)" = "mcr",
    "RuBisCO/CBB cycle (rbcL)"       = "rbcL"
  ),
  "Phototrophy" = list(
    "Photosystem-based phototrophy (psaA)" = "psaA",
    "Photosystem-based phototrophy (psbA)" = "psbA",
    "Rhodopsin-based phototrophy (rho)"    = "RHO"
  ),
  "Others" = list(
    "Arsenate reduction (arsC)" = "arsC"
  ),
  "Volatile organic compound metabolism" = list(
    "Aerobic hydrocarbon oxidation (hmoA)" = "hmoA",
    "Aerobic alkane oxidation (pxmA)"      = "pxmA",
    "Soluble diiron monooxygenase (smoX)"  = "smoX"
  )
)

# Functions searched but not shown in Fig. 3
drop_labels <- c("Arsenate reduction (arsC)", "anammox (hzsA)",
                 "Aerobic hydrocarbon oxidation (hmoA)", "Aerobic alkane oxidation (pxmA)",
                 "Soluble diiron monooxygenase (smoX)")

label_gene_map <- imap_dfr(metab_groups, function(cat_list, cat_name) {
  imap_dfr(cat_list, function(genes, label) tibble(metabolism = cat_name, label = label, gene = genes))
}) %>%
  filter(!label %in% drop_labels)

missing_genes <- setdiff(label_gene_map$gene, colnames(genome_genetable))
if (length(missing_genes)) message("Not in the MAG count matrix: ", paste(missing_genes, collapse = ", "))


# ==============================================================================
# 3. % OF MAGs ENCODING EACH FUNCTION, PER HABITAT AND PHYLUM (Fig. 3a)
# ==============================================================================

mag_groups <- genome_metadata %>%
  select(MAG, Source = SOURCE, Phylum_Class) %>%
  filter(MAG %in% rownames(genome_genetable)) %>%
  mutate(group = paste(Source, Phylum_Class, sep = " | "))

label_presence <- label_gene_map %>%
  filter(gene %in% colnames(genome_genetable)) %>%
  group_by(metabolism, label) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  mutate(pres = map(genes, function(gs) {
    tibble(MAG = rownames(genome_genetable),
           pres = as.integer(rowSums(genome_genetable[, gs, drop = FALSE]) > 0))
  })) %>%
  select(metabolism, label, pres) %>%
  unnest(pres)

presence_long <- label_presence %>%
  inner_join(mag_groups, by = "MAG") %>%
  group_by(metabolism, label, group, Source, Phylum_Class) %>%
  summarise(pct = 100 * mean(pres), .groups = "drop")

# Rows in metab_groups order; columns sediment first, then by phylum
row_info <- label_gene_map %>%
  distinct(metabolism, label) %>%
  filter(label %in% presence_long$label)

group_meta <- presence_long %>%
  distinct(group, Source, Phylum_Class) %>%
  mutate(Source = factor(Source, levels = c("Sediment", "Seawater"))) %>%
  arrange(Source, Phylum_Class)

pct_mat <- presence_long %>%
  select(label, group, pct) %>%
  pivot_wider(names_from = group, values_from = pct, values_fill = 0) %>%
  column_to_rownames("label") %>%
  as.matrix()
pct_mat <- pct_mat[row_info$label, group_meta$group, drop = FALSE]

row_metab    <- factor(row_info$metabolism, levels = names(metab_groups))
group_counts <- count(mag_groups, group)$n[match(group_meta$group, count(mag_groups, group)$group)]


# ==============================================================================
# 4. GENOME-CENTRIC: COPY NUMBER AND MANN-WHITNEY U (Fig. 3b)
# ==============================================================================
# 5.08: group 1 = seawater, group 2 = sediment. log2 fold change shown as
# log2(seawater / sediment), so negative values indicate enrichment in sediments.

manw_label <- label_gene_map %>%
  left_join(manw_result %>%
              select(gene, group_1_mean, group_2_mean, corrected_pvalue) %>%
              mutate(logFC = log2((group_1_mean + 1e-3) / (group_2_mean + 1e-3))),
            by = "gene") %>%
  group_by(metabolism, label) %>%
  summarise(logFC         = mean(logFC, na.rm = TRUE),
            padj          = suppressWarnings(min(corrected_pvalue, na.rm = TRUE)),
            mean_seawater = sum(group_1_mean, na.rm = TRUE),
            mean_sediment = sum(group_2_mean, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(logFC = replace(logFC, !is.finite(logFC), 0),
         padj  = replace(padj, !is.finite(padj), 1))
manw_label <- manw_label[match(row_info$label, manw_label$label), ]


# ==============================================================================
# 5. GENE-CENTRIC: RELATIVE ABUNDANCE AND MaAsLin2/LinDA (Fig. 3c)
# ==============================================================================

gene_da <- gene_enrichment %>%
  mutate(consensus_sig = maaslin_SedSea_q <= q_cutoff & linda_SedSea_q <= q_cutoff &
           sign(maaslin_SedSea_coef) == sign(linda_SedSea_coef)) %>%
  select(gene, MeanRelAbund_Sediment, MeanRelAbund_Seawater, consensus_sig)

label_da <- label_gene_map %>%
  left_join(gene_da, by = "gene") %>%
  group_by(metabolism, label) %>%
  summarise(relabund_sediment = mean(MeanRelAbund_Sediment, na.rm = TRUE),
            relabund_seawater = mean(MeanRelAbund_Seawater, na.rm = TRUE),
            consensus_sig     = any(consensus_sig, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(relabund_sediment, relabund_seawater), ~ replace(.x, !is.finite(.x), 0)))
label_da <- label_da[match(row_info$label, label_da$label), ]

write.xlsx(list(pct_MAGs = rownames_to_column(as.data.frame(pct_mat), "label"),
                genome_centric = manw_label,
                gene_centric = label_da),
           file.path(out_dir, "Fig3a-c_values.xlsx"))


# ==============================================================================
# 6. HEATMAP
# ==============================================================================

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

col_fun      <- colorRamp2(c(0, 0.001, 50, 100), c("#d3d3d3", "white", "#40347c", "#be3831"))
pval_col_fun <- colorRamp2(c(0, 0.05, 1), c("#F0F921FF", "#0D0887FF", "white"))
sig_cols     <- c(`FALSE` = "#f0f0f0", `TRUE` = "#be3831")

metab_colors <- microViz::distinct_palette(n = nlevels(droplevels(row_metab)), pal = "kelly", add = NA)
names(metab_colors) <- levels(droplevels(row_metab))

# Side bars (square-root scale; sediment drawn to the left)
copy_sed  <- -sqrt(manw_label$mean_sediment)
copy_sea  <-  sqrt(manw_label$mean_seawater)
abund_sed <- -sqrt(label_da$relabund_sediment)
abund_sea <-  sqrt(label_da$relabund_seawater)
ylim_abu  <- max(abs(c(abund_sed, abund_sea)))
sig_fill  <- sig_cols[as.character(label_da$consensus_sig)]

bar_anno <- rowAnnotation(
  copies_sediment = anno_barplot(copy_sed, gp = gpar(fill = "#33a02c", col = "black"),
                                 width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
                                 ylim = c(min(copy_sed), 0)),
  copies_seawater = anno_barplot(copy_sea, gp = gpar(fill = "#1f78b4", col = "black"),
                                 width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
                                 ylim = c(0, max(copy_sea))),
  log2FC_MWU      = anno_barplot(manw_label$logFC, gp = gpar(fill = pval_col_fun(manw_label$padj), col = "black"),
                                 width = unit(2.5, "cm"), bar_width = 0.8, border = TRUE,
                                 ylim = range(manw_label$logFC)),
  relabund_sediment = anno_barplot(abund_sed, gp = gpar(fill = sig_fill, col = "black"),
                                   width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
                                   ylim = c(-ylim_abu, 0)),
  relabund_seawater = anno_barplot(abund_sea, gp = gpar(fill = sig_fill, col = "black"),
                                   width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
                                   ylim = c(0, ylim_abu))
)

left_anno <- rowAnnotation(Metabolism = row_metab, col = list(Metabolism = metab_colors),
                           show_legend = FALSE, annotation_name_side = "top")

top_anno <- HeatmapAnnotation(
  MAG_count    = anno_barplot(group_counts, gp = gpar(fill = "grey30", col = "black"),
                              bar_width = 0.9, border = TRUE),
  Source       = group_meta$Source,
  Phylum_Class = group_meta$Phylum_Class,
  col = list(Source = source_colors, Phylum_Class = phylum_colors),
  annotation_name_side = "left", show_legend = FALSE
)

ht <- Heatmap(
  pct_mat, name = "% MAGs", col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 0.1, fill = fill))
  },
  top_annotation = top_anno, right_annotation = bar_anno, left_annotation = left_anno,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12),
  row_split = row_metab, column_split = group_meta$Source,
  border = TRUE, row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
  show_heatmap_legend = FALSE
)

shown_phyla <- intersect(names(phylum_colors), unique(group_meta$Phylum_Class))
legends <- list(
  Legend(title = "Phylum (class for Pseudomonadota)", at = shown_phyla,
         legend_gp = gpar(fill = phylum_colors[shown_phyla]), ncol = 8),
  Legend(title = "Mann-Whitney U adjusted p", col_fun = pval_col_fun, at = c(0, 0.05, 1),
         direction = "horizontal", legend_width = unit(4, "cm")),
  Legend(title = "MaAsLin2 and LinDA (adjusted p <= 0.05)", at = c("FALSE", "TRUE"),
         labels = c("Not significant", "Both significant"),
         legend_gp = gpar(fill = sig_cols), direction = "horizontal", ncol = 2),
  Legend(title = "Metabolism", at = names(metab_colors),
         legend_gp = gpar(fill = metab_colors), direction = "horizontal", ncol = 5)
)
heatmap_legend <- Legend(title = "% MAGs", col_fun = col_fun, direction = "horizontal",
                         legend_width = unit(4, "cm"))

svglite(file.path(out_dir, "Fig3a-c_marker_gene_heatmap.svg"), width = 20, height = 10)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     heatmap_legend_list = list(heatmap_legend))
dev.off()

svglite(file.path(out_dir, "Fig3a-c_marker_gene_heatmap.withlegend.svg"), width = 20, height = 10)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     heatmap_legend_list = list(heatmap_legend), annotation_legend_list = legends)
dev.off()

