#!/usr/bin/env Rscript
# =============================================================================
# 5.08_mag_gene_enrichment_heatmap.R
#
# Genome-centric metabolic marker genes in seawater and sediment MAGs
#   1. Mann-Whitney U test per gene (Seawater vs Sediment MAGs) with BH FDR
#   2. Heatmap of the percentage of MAGs carrying each metabolism, by
#      Source x Phylum_Class, with side barplots of:
#        - mean gene count per MAG (sediment, seawater)
#        - Mann-Whitney log2 fold change (coloured by adjusted p)
#        - gene-centric mean relative abundance (sediment, seawater), coloured
#          by MaAsLin2 + LinDA consensus significance
#
# Input:
#   2077seased_MAGs_customdb.hyddb.xlsx   MAG x gene count matrix (rows = MAGs)
#   table_S5.xlsx                         MAG metadata (curated taxonomy, SOURCE,
#                                         Phylum_Class, EnrichM_In)
#   GeneNameUpdater.xlsx                  GeneName -> Gene and GeneMetabolism
#   MAASLIN_LINDA_greeningDB_June30.xlsx  gene-centric MaAsLin2 / LinDA results
#
# Output (out_dir):
#   manw_enrichment.csv, countmat_UsedInEnrichment.csv, maginfo_UsedInEnrichment.csv
#   greening_mag_heatmap.svg, greening_mag_heatmap.withlegend.svg
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(svglite)
  library(microViz)   # distinct_palette()
})

# -----------------------------------------------------------------------------
# 0.  Paths and parameters
# -----------------------------------------------------------------------------
imos_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS"
comp_dir <- file.path(imos_dir, "0_comp_prot_newdb")
out_dir  <- file.path(comp_dir, "genome_centric")

count_file      <- file.path(comp_dir, "genome_centric", "2077seased_MAGs_customdb.hyddb.xlsx")
metadata_file   <- file.path(imos_dir, "imosmanuscript", "3_Sup_Table", "ARCH", "table_S5.xlsx")
metadata_sheet  <- "curatedtaxonomy_Sed1201_Sea876"
gene_names_file <- file.path(comp_dir, "genome_centric", "GeneNameUpdater.xlsx")
gene_da_file    <- file.path(comp_dir, "gene_centric", "MAASLIN_LINDA_greeningDB_June30.xlsx")

outfile_enrichment <- file.path(out_dir, "manw_enrichment.csv")
outfile_countmat   <- file.path(out_dir, "countmat_UsedInEnrichment.csv")
outfile_maginfo    <- file.path(out_dir, "maginfo_UsedInEnrichment.csv")

g1 <- "Seawater"
g2 <- "Sediment"
q_cutoff <- 0.05

source_colors <- c(Sediment = "#33a02c", Seawater = "#1f78b4")

phylum_colors <- c(
  "Asgardarchaeota"="#A6CEE3","Thermoplasmatota"="#B2DF8A","Thermoproteota"="#33A02B",
  "4484-113"="#FB9A99","Acidobacteriota"="#E31A1C","Actinomycetota"="#FDBF6F",
  "Alphaproteobacteria"="#FF7F00","Bacillota"="#6A3D9A","Bacteroidota"="#FFFF99",
  "Bacteroidota_A"="#B15928","Bdellovibrionota"="#1ff8ff","Bdellovibrionota_B"="#1B9E77",
  "Bipolaricaulota"="#D95F02","Campylobacterota"="#7570B3","Chlamydiota"="#E7298A",
  "Chloroflexota"="#66A61E","Cyanobacteriota"="#E6AB02","Desulfobacterota"="#A6761D",
  "Desulfobacterota_D"="#666666","Gammaproteobacteria"="#b249d5","Gemmatimonadota"="#4b6a53",
  "Hydrogenedentota"="#7edc45","Krumholzibacteriota"="#5c47b8","Latescibacterota"="#cfd251",
  "Margulisbacteria"="#ff69b4","Marinisomatota"="#cd3e50","Myxococcota"="#69c86c",
  "Myxococcota_A"="#da6130","Nitrospirota"="#5e79b2","Patescibacteriota"="#c29545",
  "Planctomycetota"="#532a5a","SAR324"="#5f7b35","Spirochaetota"="#c497cf",
  "Tectomicrobia"="#773a27","Verrucomicrobiota"="#7cb9cb","Nanobdellota"="#1F78B1",
  "Babelota"="#CAB2D6","Muiribacteriota"="#83d5af"
)

# Metabolism labels and their marker genes. Gene names use make.names()
# form, e.g. "[NiFe]_1a" is "X.NiFe._1a".
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
    "Denitrification:Nitrite oxide reduction (nosZ)" = "nosZ",             # N2O to N2
    "DNRA (NrfA)"                                    = "nrfA",             # NO2- to NH4+
    "anammox (hzsA)"                                 = "hzsA",             # NH4+ to N2
    "Nitrogen fixation (nifH)"                       = "nifH",
    "Nitrification:Ammonia oxidation (amoA)"         = "amoA",
    "Nitrification:Nitrite oxidation (NxrA)"         = "nxrA"              # NO2- to NO3-
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
    "[NiFe] G1 (H2 uptake)"         = c("X.NiFe._1a", "X.NiFe._1b", "X.NiFe._1c", "X.NiFe._1d",
                                        "X.NiFe._1e", "X.NiFe._1f", "X.NiFe._1h", "X.NiFe._1k",
                                        "X.NiFe._1l"),
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
    "3-hydroxybutyrate cycle (mcr)"  = "mcr",
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

# Labels left out of the heatmap
drop_labels <- c(
  "Arsenate reduction (arsC)",
  "anammox (hzsA)",
  "Aerobic hydrocarbon oxidation (hmoA)",
  "Aerobic alkane oxidation (pxmA)",
  "Soluble diiron monooxygenase (smoX)"
)

# =============================================================================
# PART 1. Mann-Whitney U enrichment (Seawater vs Sediment MAGs)
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1  Gene name lookup (Greening DB): GeneName -> Gene, Gene -> GeneMetabolism
# -----------------------------------------------------------------------------
message("Loading count matrix and metadata")

annot_update <- read.xlsx(gene_names_file, rowNames = FALSE)
gene_lookup       <- annot_update %>% distinct(GeneName, Gene)
metabolism_lookup <- annot_update %>% distinct(GeneMetabolism, Gene)

dupe <- gene_lookup %>% count(GeneName) %>% filter(n > 1)
if (nrow(dupe) > 0) {
  warning("GeneName -> Gene mapping is not 1:1; left_join will duplicate rows for: ",
          paste(dupe$GeneName, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 1.2  MAG x gene counts, renamed via the lookup (generic NiFe/FeFe dropped)
# -----------------------------------------------------------------------------
count_mat <- read.xlsx(count_file, rowNames = TRUE) %>%
  rownames_to_column("MAG") %>%
  dplyr::select(-any_of(c("NiFe_hydrogenase", "FeFe_hydrogenase"))) %>%
  pivot_longer(-MAG, names_to = "GeneName", values_to = "count") %>%
  left_join(gene_lookup, by = "GeneName") %>%
  mutate(Gene = coalesce(Gene, GeneName)) %>%   # keep original name if unmatched
  group_by(MAG, Gene) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = count, values_fill = 0) %>%
  column_to_rownames("MAG")

# -----------------------------------------------------------------------------
# 1.3  MAG metadata (EnrichM_In == TRUE) and alignment with count_mat
# -----------------------------------------------------------------------------
mag_info <- read.xlsx(metadata_file, sheet = metadata_sheet, rowNames = TRUE) %>%
  filter(as.logical(EnrichM_In)) %>%
  rownames_to_column("MAG")

stopifnot("SOURCE" %in% names(mag_info))
print(table(mag_info$SOURCE, useNA = "always"))

shared_mags   <- intersect(rownames(count_mat), mag_info$MAG)
only_in_count <- setdiff(rownames(count_mat), mag_info$MAG)
only_in_meta  <- setdiff(mag_info$MAG, rownames(count_mat))

message("count_mat MAGs: ", nrow(count_mat),
        " | mag_info MAGs (EnrichM_In == TRUE): ", nrow(mag_info),
        " | shared: ", length(shared_mags))
if (length(only_in_count) > 0) message("In count_mat only (head): ", paste(head(only_in_count, 5), collapse = ", "))
if (length(only_in_meta)  > 0) message("In mag_info only (head): ",  paste(head(only_in_meta,  5), collapse = ", "))

count_mat <- count_mat[shared_mags, , drop = FALSE]
mag_info  <- mag_info %>% filter(MAG %in% shared_mags)

g1_mags <- mag_info$MAG[mag_info$SOURCE == g1]
g2_mags <- mag_info$MAG[mag_info$SOURCE == g2]
message(g1, ": ", length(g1_mags), " MAGs; ", g2, ": ", length(g2_mags), " MAGs")
stopifnot(length(g1_mags) > 0, length(g2_mags) > 0)

# -----------------------------------------------------------------------------
# 1.4  Long format; drop genes that are zero in every MAG
# -----------------------------------------------------------------------------
long_df <- count_mat %>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG, names_to = "gene", values_to = "count") %>%
  mutate(SOURCE = case_when(MAG %in% g1_mags ~ g1,
                            MAG %in% g2_mags ~ g2)) %>%
  filter(!is.na(SOURCE))

genes_keep <- long_df %>%
  group_by(gene) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(gene)

dropped_n <- n_distinct(long_df$gene) - length(genes_keep)
if (dropped_n > 0) message("Dropping ", dropped_n, " all-zero gene(s) from the test")
long_df <- long_df %>% filter(gene %in% genes_keep)

# -----------------------------------------------------------------------------
# 1.5  Mann-Whitney U per gene, log2 fold change, BH correction
#      (all-zero genes were dropped, so every gene has signal in >= 1 group)
# -----------------------------------------------------------------------------
message("Mann-Whitney U test per gene")

mw <- function(x, y) suppressWarnings(wilcox.test(x, y, exact = FALSE))

manw_result <- long_df %>%
  group_by(gene) %>%
  summarise(
    g1_mean = mean(count[SOURCE == g1]) + 1e-10,
    g2_mean = mean(count[SOURCE == g2]) + 1e-10,
    wt      = list(mw(count[SOURCE == g1], count[SOURCE == g2])),
    .groups = "drop"
  ) %>%
  mutate(
    score              = map_dbl(wt, ~ unname(.x$statistic)),
    pvalue             = map_dbl(wt, ~ .x$p.value),
    corrected_pvalue   = p.adjust(pvalue, method = "BH"),
    `FoldChange_g1/g2` = log2(g1_mean / g2_mean),
    `FoldChange_g2/g1` = log2(g2_mean / g1_mean),
    enriched_in        = case_when(g1_mean > g2_mean ~ g1,
                                   g1_mean < g2_mean ~ g2),
    group_1 = g1,
    group_2 = g2
  ) %>%
  left_join(metabolism_lookup, by = c("gene" = "Gene")) %>%
  mutate(GeneMetabolism = coalesce(GeneMetabolism, "Unassigned")) %>%
  dplyr::select(gene, GeneMetabolism, group_1, group_2, enriched_in,
                group_1_mean = g1_mean, group_2_mean = g2_mean,
                `FoldChange_g1/g2`, `FoldChange_g2/g1`,
                score, pvalue, corrected_pvalue)

# -----------------------------------------------------------------------------
# 1.6  Write tables
# -----------------------------------------------------------------------------
write.csv(manw_result, outfile_enrichment, row.names = FALSE)
write.csv(count_mat,   outfile_countmat,   row.names = TRUE)
write.csv(mag_info,    outfile_maginfo,    row.names = FALSE)
message("Wrote: ", outfile_enrichment, "\n       ", outfile_countmat, "\n       ", outfile_maginfo)

# =============================================================================
# PART 2. Heatmap: % MAGs carrying each metabolism, by Source x Phylum_Class
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1  Inputs, with gene names in make.names() form to match metab_groups
# -----------------------------------------------------------------------------
genome_genetable <- count_mat
colnames(genome_genetable) <- make.names(colnames(genome_genetable))

genome_metadata <- mag_info

manw_result <- manw_result %>% mutate(gene = make.names(gene))

gene_enrichment <- read.xlsx(gene_da_file, rowNames = TRUE) %>%
  rownames_to_column("gene") %>%
  mutate(gene = make.names(gene))

# -----------------------------------------------------------------------------
# 2.2  Summary: % of MAGs with any hydrogenase ([NiFe] or [FeFe]) hit
# -----------------------------------------------------------------------------
hyd_genes <- grep("^X\\.(NiFe|FeFe)\\.", colnames(genome_genetable), value = TRUE)
has_hyd   <- rowSums(genome_genetable[, hyd_genes, drop = FALSE]) > 0

genome_metadata %>%
  mutate(has_hydrogenase = has_hyd[MAG]) %>%
  group_by(SOURCE) %>%
  summarise(n_MAGs = n(),
            pct_with_hydrogenase = 100 * mean(has_hydrogenase), .groups = "drop") %>%
  print()

# -----------------------------------------------------------------------------
# 2.3  Label <-> gene lookup and per-MAG label presence (any gene > 0)
# -----------------------------------------------------------------------------
label_gene_map <- imap_dfr(metab_groups, \(cat_list, cat_name) {
  imap_dfr(cat_list, \(genes, label) tibble(metabolism = cat_name,
                                            label      = label,
                                            gene       = genes))
})

missing_genes <- setdiff(unique(label_gene_map$gene), colnames(genome_genetable))
if (length(missing_genes)) {
  message("Missing from genome_genetable: ", paste(missing_genes, collapse = ", "))
}

mag_label_presence <- label_gene_map %>%
  filter(gene %in% colnames(genome_genetable)) %>%
  group_by(metabolism, label) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  mutate(pres = map(genes, \(gs) {
    tibble(MAG  = rownames(genome_genetable),
           pres = as.integer(rowSums(genome_genetable[, gs, drop = FALSE]) > 0))
  })) %>%
  dplyr::select(metabolism, label, pres) %>%
  unnest(pres)

# -----------------------------------------------------------------------------
# 2.4  % MAGs per Source x Phylum_Class group
# -----------------------------------------------------------------------------
mag_groups <- genome_metadata %>%
  dplyr::select(MAG, Source = SOURCE, Phylum_Class) %>%
  filter(MAG %in% rownames(genome_genetable)) %>%
  mutate(group = paste(Source, Phylum_Class, sep = " | "))

presence_long <- mag_label_presence %>%
  inner_join(mag_groups, by = "MAG") %>%
  group_by(metabolism, label, group, Source, Phylum_Class) %>%
  summarise(pct = mean(pres) * 100, .groups = "drop")

mat_wide <- presence_long %>%
  dplyr::select(metabolism, label, group, pct) %>%
  pivot_wider(names_from = group, values_from = pct, values_fill = 0) %>%
  filter(!label %in% drop_labels)

# Rows in metab_groups order
label_order <- label_gene_map %>% distinct(metabolism, label)
mat_wide <- mat_wide %>%
  slice(order(match(label, label_order$label)))

# Columns: Sediment first, then Phylum_Class alphabetically
group_meta <- presence_long %>%
  distinct(group, Source, Phylum_Class) %>%
  mutate(Source = factor(Source, levels = c("Sediment", "Seawater"))) %>%
  arrange(Source, Phylum_Class)
col_order <- group_meta$group

heat_mat <- mat_wide %>%
  dplyr::select(all_of(col_order)) %>%
  as.matrix()
rownames(heat_mat) <- mat_wide$label

row_metab    <- factor(mat_wide$metabolism, levels = names(metab_groups))
group_counts <- mag_groups %>% count(group) %>% slice(match(col_order, group))

# -----------------------------------------------------------------------------
# 2.5  Gene-centric MaAsLin2 + LinDA (Sediment vs Seawater), per label
#      Consensus = both q <= q_cutoff with the same direction; a label is
#      significant if ANY of its genes is. Mean relative abundance is
#      square-root transformed; sediment is plotted to the left (negative).
# -----------------------------------------------------------------------------
gene_da <- gene_enrichment %>%
  mutate(consensus_sig = maaslin_SedSea_q <= q_cutoff & linda_SedSea_q <= q_cutoff &
           sign(maaslin_SedSea_coef) == sign(linda_SedSea_coef)) %>%
  dplyr::select(gene, MeanRelAbund_Sediment, MeanRelAbund_Seawater, consensus_sig)

label_da <- label_gene_map %>%
  left_join(gene_da, by = "gene") %>%
  group_by(metabolism, label) %>%
  summarise(
    MeanRelAbund_Sediment = sqrt(mean(MeanRelAbund_Sediment, na.rm = TRUE)),
    MeanRelAbund_Seawater = sqrt(mean(MeanRelAbund_Seawater, na.rm = TRUE)),
    consensus_sig         = any(consensus_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  slice(match(rownames(heat_mat), label)) %>%
  mutate(across(c(MeanRelAbund_Sediment, MeanRelAbund_Seawater),
                ~ if_else(is.finite(.x), .x, 0)),
         MeanRelAbund_Sediment = -MeanRelAbund_Sediment)

# -----------------------------------------------------------------------------
# 2.6  Genome-centric Mann-Whitney results, per label
#      logFC = log2(Seawater / Sediment) of mean gene count per MAG;
#      p = smallest adjusted p among the label's genes.
# -----------------------------------------------------------------------------
manw_label <- label_gene_map %>%
  left_join(manw_result %>%
              dplyr::select(gene, group_1_mean, group_2_mean, corrected_pvalue) %>%
              mutate(logFC = log2((group_1_mean + 1e-3) / (group_2_mean + 1e-3))),
            by = "gene") %>%
  group_by(metabolism, label) %>%
  summarise(logFC   = mean(logFC, na.rm = TRUE),
            pval    = suppressWarnings(min(corrected_pvalue, na.rm = TRUE)),
            g1_mean = sum(group_1_mean, na.rm = TRUE),   # seawater
            g2_mean = sum(group_2_mean, na.rm = TRUE),   # sediment
            .groups = "drop") %>%
  slice(match(rownames(heat_mat), label)) %>%
  mutate(logFC = if_else(is.finite(logFC), logFC, 0),
         pval  = if_else(is.finite(pval),  pval,  1))

genecounts_sediment <- -sqrt(manw_label$g2_mean)
genecounts_seawater <-  sqrt(manw_label$g1_mean)

# -----------------------------------------------------------------------------
# 2.7  Colour scales
# -----------------------------------------------------------------------------
col_fun      <- colorRamp2(c(0, 0.001, 50, 100), c("#d3d3d3", "white", "#40347c", "#be3831"))
pval_col_fun <- colorRamp2(c(0, 0.05, 1), c("#F0F921FF", "#0D0887FF", "white"))
sig_col_fun  <- colorRamp2(c(0, 1), c("#f0f0f0", "#be3831"))

metab_levels <- levels(droplevels(row_metab))
metab_colors <- setNames(distinct_palette(n = length(metab_levels), pal = "kelly", add = NA),
                         metab_levels)

ylim_abu <- max(abs(c(label_da$MeanRelAbund_Sediment, label_da$MeanRelAbund_Seawater)), na.rm = TRUE)

# -----------------------------------------------------------------------------
# 2.8  Annotations
# -----------------------------------------------------------------------------
row_anno <- rowAnnotation(
  Metabolism = row_metab,
  col = list(Metabolism = metab_colors),
  show_legend = FALSE,
  annotation_name_side = "top"
)

bar_anno <- rowAnnotation(
  genecount_sediment = anno_barplot(
    genecounts_sediment,
    gp = gpar(fill = source_colors[["Sediment"]], col = "black"),
    width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
    ylim = c(min(genecounts_sediment, na.rm = TRUE), 0)),
  genecount_seawater = anno_barplot(
    genecounts_seawater,
    gp = gpar(fill = source_colors[["Seawater"]], col = "black"),
    width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
    ylim = c(0, max(genecounts_seawater, na.rm = TRUE))),
  logFC_MANW = anno_barplot(
    manw_label$logFC,
    gp = gpar(fill = pval_col_fun(manw_label$pval), col = "black"),
    width = unit(2.5, "cm"), bar_width = 0.8, border = TRUE,
    ylim = range(manw_label$logFC, na.rm = TRUE)),
  sed_relabu = anno_barplot(
    label_da$MeanRelAbund_Sediment,
    gp = gpar(fill = sig_col_fun(as.integer(label_da$consensus_sig)), col = "black"),
    width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
    ylim = c(-ylim_abu, 0)),
  sea_relabu = anno_barplot(
    label_da$MeanRelAbund_Seawater,
    gp = gpar(fill = sig_col_fun(as.integer(label_da$consensus_sig)), col = "black"),
    width = unit(1.25, "cm"), bar_width = 0.8, border = TRUE,
    ylim = c(0, ylim_abu))
)

top_anno <- HeatmapAnnotation(
  MAG_count = anno_barplot(group_counts$n,
                           gp = gpar(fill = "grey30", col = "black"),
                           bar_width = 0.9, border = TRUE),
  Source       = group_meta$Source,
  Phylum_Class = group_meta$Phylum_Class,
  col = list(Source = source_colors, Phylum_Class = phylum_colors),
  annotation_name_side = "left", show_legend = FALSE
)

# -----------------------------------------------------------------------------
# 2.9  Heatmap and legends
# -----------------------------------------------------------------------------
ht <- Heatmap(
  matrix = heat_mat,
  name   = "% MAGs with pathway",
  col    = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 0.1, fill = fill))
  },
  top_annotation      = top_anno,
  right_annotation    = bar_anno,
  left_annotation     = row_anno,
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 12),
  column_names_gp     = gpar(fontsize = 12),
  row_split           = row_metab,
  column_split        = group_meta$Source,
  border              = TRUE,
  row_gap             = unit(1, "mm"),
  column_gap          = unit(1, "mm"),
  show_heatmap_legend = FALSE
)

phyla_present <- intersect(names(phylum_colors), unique(group_meta$Phylum_Class))

heatmap_legend <- Legend(title = "% MAGs", col_fun = col_fun,
                         direction = "horizontal", legend_width = unit(4, "cm"))
phylum_legend  <- Legend(title = "Phylum_Class", at = phyla_present,
                         legend_gp = gpar(fill = phylum_colors[phyla_present]), ncol = 8)
pval_legend    <- Legend(title = "MANW adj.p", col_fun = pval_col_fun,
                         at = c(0, 0.05, 1), labels = c("0", "0.05", "1"),
                         direction = "horizontal", legend_width = unit(4, "cm"))
sig_legend     <- Legend(title = paste0("MaAsLin+LINDA consensus (q≤", q_cutoff, ")"),
                         at = c(0, 1), labels = c("Not significant", "Both significant"),
                         legend_gp = gpar(fill = c("#f0f0f0", "#be3831")),
                         direction = "horizontal", ncol = 2)
metab_legend   <- Legend(title = "Metabolism", at = names(metab_colors),
                         legend_gp = gpar(fill = metab_colors),
                         direction = "horizontal", ncol = 5)

# -----------------------------------------------------------------------------
# 2.10 Draw
# -----------------------------------------------------------------------------
svglite(file.path(out_dir, "greening_mag_heatmap.svg"), width = 20, height = 20)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     heatmap_legend_list = list(heatmap_legend))
dev.off()

svglite(file.path(out_dir, "greening_mag_heatmap.withlegend.svg"), width = 20, height = 20)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     heatmap_legend_list    = list(heatmap_legend),
     annotation_legend_list = list(phylum_legend, pval_legend, sig_legend, metab_legend))
dev.off()

message("Done.")