# ==============================================================================
# 5.14_gcf_bigfam_distance.R
#
# Gene cluster families (GCFs) of sediment and seawater MAGs compared with the
# OMD and RMD, and their novelty relative to BiG-FAM (Fig. 5).
#   1. Setup and input
#   2. GCFs per MAG in each database
#   3. GCF composition by BGC class
#   4. Mean distance to BiG-FAM per GCF and novel GCFs (mean distance > 0.4)
#   5. Histograms of mean distance per GCF, coloured by BGC class
#
# GCFs are from clust-o-matic (5.11) on dereplicated datasets, BGCs >= 5,000 bp.
# antiSMASH products are grouped into BGC classes (sheet antismash_class). For
# each GCF, the distance to BiG-FAM is the mean of its BGCs' BiG-SLiCE distances
# (5.13); the GCF's BGC class and phylum are those of its longest BGC.
#
# Input : workbook BGCs_GCFs_IMOS_OMD.xlsx with sheets
#           GCF_sum_derep    : GCFs per MAG per database (Derep, GCFs per MAG (5000 bp))
#           COM_sum_derep    : one row per BGC (SOURCE, BGC_FUNCTION, GCF_ID)
#           antismash_class  : antiSMASH product (Antismash) to BGC class (Summary)
#           gcf_dist_0.4_SED, gcf_dist_0.4_SEA : one row per BGC (clustomatic_gcf,
#             bs_dist, BGC_FUNCTION, len_nt, Phylum)
# Output: Fig. 5 panels (svg); GCF composition and novel GCF counts (xlsx)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(patchwork)
  library(tidyverse)
})

infile_bs <- "<path/to/BGCs_GCFs_IMOS_OMD.xlsx>"
out_dir   <- "<path/to/5_MAG_annotation/7_BGCs/figures>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

novel_cutoff <- 0.4
binlen       <- 0.05

db_levels <- c("MGD_Sediment", "MGD_Seawater", "OMD", "RMD")
db_names  <- c(IMOS_SED = "MGD_Sediment", IMOS_SEA = "MGD_Seawater",
               OMD_DB = "OMD", `RMDv1-antismash_GBK` = "RMD")
db_colors <- c(MGD_Sediment = "#33a02c", MGD_Seawater = "#1f78b4",
               OMD = "#6AABC6", RMD = "#C73C74")

bgc_colors <- c(
  "MIBIG"   = "#474747",
  "NRPS"    = "#B5D465",
  "Others"  = "#FAD753",
  "RiPPs"   = "#EB8F6B",
  "T1PKS"   = "#8EA0C7",
  "T2/3PKS" = "#D790C1",
  "Terpene" = "#82BFA6"
)

rename_db <- function(x) factor(coalesce(db_names[x], x), levels = db_levels)

antismash_category <- read.xlsx(infile_bs, sheet = "antismash_class")

# antiSMASH products grouped into BGC classes
add_bgc_class <- function(df) {
  df %>%
    left_join(antismash_category, by = c("BGC_FUNCTION" = "Antismash")) %>%
    mutate(BGC_FUNCTION = coalesce(Summary, BGC_FUNCTION)) %>%
    select(-Summary)
}

save_svg <- function(plot, name, width, height) {
  ggsave(file.path(out_dir, name), plot, width = width, height = height, bg = "white")
}


# ==============================================================================
# 2. GCFs PER MAG
# ==============================================================================

gcf_summary <- read.xlsx(infile_bs, sheet = "GCF_sum_derep") %>%
  mutate(db = rename_db(Derep))

p_gcf_per_mag <- ggplot(gcf_summary, aes(db, `GCFs.per.MAG.(5000.bp)`, fill = db)) +
  geom_col(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", `GCFs.per.MAG.(5000.bp)`)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = db_colors) +
  labs(x = NULL, y = "Number of GCFs per MAG") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save_svg(p_gcf_per_mag, "GCFs_per_MAG.svg", width = 5, height = 5)


# ==============================================================================
# 3. GCF COMPOSITION BY BGC CLASS
# ==============================================================================
# Each GCF counted once per database and BGC class; MIBiG excluded

gcf_composition <- read.xlsx(infile_bs, sheet = "COM_sum_derep") %>%
  filter(SOURCE != "MIBIG") %>%
  add_bgc_class() %>%
  mutate(SOURCE = rename_db(SOURCE)) %>%
  distinct(SOURCE, BGC_FUNCTION, GCF_ID) %>%
  count(SOURCE, BGC_FUNCTION, name = "count") %>%
  group_by(SOURCE) %>%
  mutate(proportion = 100 * count / sum(count)) %>%
  ungroup()

unmapped <- setdiff(unique(gcf_composition$BGC_FUNCTION), names(bgc_colors))
if (length(unmapped)) warning("BGC classes without a colour: ", paste(unmapped, collapse = ", "))

# Stack order: most GCFs overall at the bottom
class_order <- gcf_composition %>%
  group_by(BGC_FUNCTION) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(BGC_FUNCTION)
gcf_composition <- mutate(gcf_composition, BGC_FUNCTION = factor(BGC_FUNCTION, levels = class_order))

p_composition <- ggplot(gcf_composition, aes(SOURCE, proportion, fill = BGC_FUNCTION)) +
  geom_col(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", proportion)), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = bgc_colors) +
  labs(x = NULL, y = "Proportion (%)", fill = "BGC class") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save_svg(p_composition, "GCF_composition_proportion.svg", width = 5, height = 5)
save_svg(p_composition + theme(legend.position = "right"), "GCF_composition_proportion_withlegend.svg",
         width = 5, height = 5)


# ==============================================================================
# 4. MEAN DISTANCE TO BiG-FAM AND NOVEL GCFs
# ==============================================================================

read_bgcs <- function(sheet) add_bgc_class(read.xlsx(infile_bs, sheet = sheet, rowNames = TRUE))

bgcs <- list(Sediment = read_bgcs("gcf_dist_0.4_SED"),
             Seawater = read_bgcs("gcf_dist_0.4_SEA"))

gcfs <- map(bgcs, ~ .x %>%
  group_by(clustomatic_gcf) %>%
  summarise(mean_bs_dist = mean(bs_dist, na.rm = TRUE),
            BGC_FUNCTION = BGC_FUNCTION[which.max(len_nt)],
            Phylum       = Phylum[which.max(len_nt)],
            .groups      = "drop"))

novel <- imap_dfr(gcfs, ~ tibble(
  Habitat    = .y,
  GCFs       = n_distinct(.x$clustomatic_gcf),
  Novel_GCFs = sum(.x$mean_bs_dist > novel_cutoff),
  Novel_pct  = 100 * Novel_GCFs / GCFs))
print(novel)

write.xlsx(list(GCF_composition = gcf_composition, novel_GCFs = novel),
           file.path(out_dir, "GCF_composition_novelty_BiG-FAM.xlsx"))


# ==============================================================================
# 5. HISTOGRAMS
# ==============================================================================

dist_hist <- function(df) {
  ggplot(df, aes(mean_bs_dist, fill = BGC_FUNCTION)) +
    geom_histogram(binwidth = binlen, position = "stack", colour = "black", linewidth = 0.3) +
    geom_vline(xintercept = novel_cutoff) +
    scale_fill_manual(values = bgc_colors) +
    labs(x = "Distance to BiG-FAM", y = "GCF count", fill = "BGC class") +
    theme_minimal() +
    theme(axis.text    = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.line    = element_line(colour = "black"),
          legend.position = "none")
}

plots <- map(gcfs, dist_hist)

# Same axes for both habitats: x up to the largest mean distance, y up to the
# tallest plotted bar
x_max <- ceiling(max(map_dbl(gcfs, ~ max(.x$mean_bs_dist))) * 100) / 100
y_max <- max(map_dbl(plots, ~ max(layer_data(.x, 1)$ymax)))

plots <- map(plots, ~ .x + coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)))
plots$Seawater <- plots$Seawater + theme(axis.title.y = element_blank())

save_svg(plots$Sediment, "GCF_dist_bigfam_sed.svg", width = 5, height = 5)
save_svg(plots$Seawater, "GCF_dist_bigfam_sea.svg", width = 5, height = 5)

sessionInfo()
