# ==============================================================================
# 5.14_gcf_bigfam_distance.R
#
# Novelty of sediment and seawater gene cluster families (GCFs) relative to
# BiG-FAM (Fig. 5).
#   1. Setup and input
#   2. Mean distance to BiG-FAM per GCF
#   3. Number and percentage of novel GCFs (mean distance > 0.4)
#   4. Histograms of mean distance per GCF, coloured by BGC class
#
# For each GCF (clust-o-matic, 5.11), the distance to BiG-FAM is the mean of
# its BGCs' BiG-SLiCE distances to their closest BiG-FAM GCF; the GCF's BGC
# class and phylum are those of its longest BGC.
#
# Input : workbook with one sheet per habitat (gcf_dist_0.4_SED,
#         gcf_dist_0.4_SEA), one row per BGC, with columns clustomatic_gcf,
#         bs_dist (BiG-SLiCE distance to BiG-FAM), BGC_FUNCTION (antiSMASH
#         product), len_nt, and Phylum; and sheet antismash_class mapping
#         antiSMASH products (Antismash) to BGC classes (Summary)
# Output: Fig. 5 distance histograms (svg); novel GCF counts (xlsx)
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

bgc_colors <- c(
  "MIBIG"   = "#474747",
  "NRPS"    = "#B5D465",
  "Others"  = "#FAD753",
  "RiPPs"   = "#EB8F6B",
  "T1PKS"   = "#8EA0C7",
  "T2/3PKS" = "#D790C1",
  "Terpene" = "#82BFA6"
)

antismash_category <- read.xlsx(infile_bs, sheet = "antismash_class")

# One row per BGC, with antiSMASH products grouped into BGC classes
read_bgcs <- function(sheet) {
  read.xlsx(infile_bs, sheet = sheet, rowNames = TRUE) %>%
    left_join(antismash_category, by = c("BGC_FUNCTION" = "Antismash")) %>%
    mutate(BGC_FUNCTION = if_else(!is.na(Summary), Summary, BGC_FUNCTION)) %>%
    select(-Summary)
}

bgcs <- list(Sediment = read_bgcs("gcf_dist_0.4_SED"),
             Seawater = read_bgcs("gcf_dist_0.4_SEA"))


# ==============================================================================
# 2. MEAN DISTANCE TO BiG-FAM PER GCF
# ==============================================================================

gcfs <- map(bgcs, ~ .x %>%
  group_by(clustomatic_gcf) %>%
  summarise(mean_bs_dist = mean(bs_dist, na.rm = TRUE),
            BGC_FUNCTION = BGC_FUNCTION[which.max(len_nt)],
            Phylum       = Phylum[which.max(len_nt)],
            .groups      = "drop"))


# ==============================================================================
# 3. NOVEL GCFs
# ==============================================================================

novel <- imap_dfr(gcfs, ~ tibble(
  Habitat     = .y,
  GCFs        = n_distinct(.x$clustomatic_gcf),
  Novel_GCFs  = sum(.x$mean_bs_dist > novel_cutoff),
  Novel_pct   = 100 * Novel_GCFs / GCFs))
novel

write.xlsx(novel, file.path(out_dir, "GCF_novelty_BiG-FAM.xlsx"))


# ==============================================================================
# 4. HISTOGRAMS
# ==============================================================================

dist_hist <- function(df) {
  ggplot(df, aes(x = mean_bs_dist, fill = BGC_FUNCTION)) +
    geom_histogram(binwidth = binlen, position = "stack", colour = "black", linewidth = 0.3) +
    geom_vline(xintercept = novel_cutoff) +
    scale_fill_manual(values = bgc_colors) +
    labs(x = "Distance to BiG-FAM", y = "GCF count", fill = "GCF class") +
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

plots$Sediment | plots$Seawater

ggsave(file.path(out_dir, "GCF_dist_bigfam_sed.svg"), plots$Sediment, width = 3, height = 2)
ggsave(file.path(out_dir, "GCF_dist_bigfam_sea.svg"), plots$Seawater, width = 3, height = 2)

sessionInfo()
