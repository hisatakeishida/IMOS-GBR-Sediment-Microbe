# ==============================================================================
# 4.15_genome_features.R
#
# Genome size, GC content, and number of predicted genes of sediment and
# seawater MAGs (Fig. 2e-g), as measured and normalised by completeness.
#   1. Setup and input
#   2. Wilcoxon rank-sum tests and habitat means
#   3. Density plots with habitat means (dashed lines), for both versions
#
# Normalised: genome size and number of predicted genes divided by CheckM2
# completeness, to estimate values for complete genomes. GC content and coding
# density are not normalised.
#
# Input : table of MAG statistics, one row per MAG (1,201 sediment and 876
#         seawater), with columns Source (Sediment or Seawater), Genome_size
#         (bp), GC_Content (%), Coding_Density_CheckM2 (fraction),
#         Total_Coding_Sequences, and Completeness_C2 (%; CheckM2)
# Output: two figures (measured and normalised); test results and means as
#         a table
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggpubr)
  library(tidyverse)
})

stats_path  <- "<path/to/Sea__Sed_indiass_coass.xlsx>"
stats_sheet <- "derep"
out_dir     <- "<path/to/4_MAG_recovery/7_taxonomy>"

habitat_cols <- c(Seawater = "#1f78b4", Sediment = "#33a02c")

magstats <- read.xlsx(stats_path, sheet = stats_sheet, rowNames = TRUE) %>%
  mutate(completeness       = Completeness_C2 / 100,
         # As measured
         Genome_size_Mb     = Genome_size / 1e6,
         Genes_k            = Total_Coding_Sequences / 1000,
         Coding_density     = Coding_Density_CheckM2 * 100,
         # Normalised by completeness
         Genome_size_est_Mb = Genome_size_Mb / completeness,
         Genes_est_k        = Genes_k / completeness)

table(magstats$Source)
summary(magstats$Completeness_C2)


# ==============================================================================
# 2. WILCOXON RANK-SUM TESTS AND HABITAT MEANS
# ==============================================================================

traits <- c("Genome_size_Mb", "Genome_size_est_Mb", "GC_Content",
            "Coding_density", "Genes_k", "Genes_est_k")

stats_table <- map_dfr(traits, function(trait) {
  means <- tapply(magstats[[trait]], magstats$Source, mean, na.rm = TRUE)
  tibble(trait         = trait,
         mean_Sediment = means[["Sediment"]],
         mean_Seawater = means[["Seawater"]],
         wilcox_p      = wilcox.test(as.formula(paste(trait, "~ Source")), data = magstats)$p.value)
})
stats_table

write.xlsx(stats_table, file.path(out_dir, "genome_features_wilcoxon.xlsx"))

habitat_means <- magstats %>%
  group_by(Source) %>%
  summarise(across(all_of(traits), ~ mean(.x, na.rm = TRUE)))


# ==============================================================================
# 3. DENSITY PLOTS
# ==============================================================================

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position  = "none",
        panel.grid       = element_blank(),
        axis.ticks       = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8))

density_plot <- function(trait, x_label) {
  ggplot(magstats, aes(x = .data[[trait]], fill = Source)) +
    geom_density(alpha = 0.7) +
    geom_vline(data = habitat_means, aes(xintercept = .data[[trait]], colour = Source),
               linetype = 2, linewidth = 1) +
    scale_fill_manual(values = habitat_cols) +
    scale_colour_manual(values = habitat_cols) +
    labs(x = x_label, y = "Density") +
    base_theme
}

p_gc <- density_plot("GC_Content", "GC content (%)")

# As measured
p_measured <- ggarrange(
  density_plot("Genome_size_Mb", "Genome size (Mb)"),
  p_gc,
  density_plot("Genes_k", "Predicted genes (thousands)"),
  ncol = 1, nrow = 3, align = "v")
p_measured
ggsave(file.path(out_dir, "Fig2e-g_genome_features_measured.svg"), p_measured, width = 5, height = 9)

# Normalised by completeness
p_normalised <- ggarrange(
  density_plot("Genome_size_est_Mb", "Estimated genome size (Mb)"),
  p_gc,
  density_plot("Genes_est_k", "Estimated predicted genes (thousands)"),
  ncol = 1, nrow = 3, align = "v")
p_normalised
ggsave(file.path(out_dir, "Fig2e-g_genome_features_normalised.svg"), p_normalised, width = 5, height = 9)

sessionInfo()