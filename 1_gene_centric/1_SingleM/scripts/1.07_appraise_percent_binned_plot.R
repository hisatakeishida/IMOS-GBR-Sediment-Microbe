# ==============================================================================
# 1.07_appraise_percent_binned_plot.R
#
# Percentage of SingleM mOTU coverage in each sample assigned to a MAG or a
# GTDB R10-RS226 genome (singlem appraise, 1.06), at species and genus level,
# for sediment and seawater, all reads and 10 M read pairs.
#   1. Setup
#   2. Read appraise tables
#   3. Scatter plot
#
# One point per sample; black bars mark the pooled value over all samples (the
# "total" row written by singlem appraise). The "total" and "average" rows and
# "#" comment lines are excluded from the per-sample values.
#
# Input files are named <NN>_<hab>_<reads>_<ref>.<rank>_appraise.tsv, e.g.
# 00_sed_allreads_mags.species_appraise.tsv. Habitat is taken from the first
# digit of NN (0 = sediment, 1 = seawater); a warning is given if the habitat
# in the file name disagrees.
#
# Input : 4_appraise_novelty/0_results/*_appraise.tsv
# Output: appraise_percent_binned.svg and .png; per-sample values (tsv)
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

results_dir <- "<path/to/1_gene_centric/1_SingleM/4_appraise_novelty/0_results>"
out_dir     <- "<path/to/1_gene_centric/1_SingleM/4_appraise_novelty>"

habitat_cols <- c(Sediment = "#33a02c", Seawater = "#1f78b4")


# ==============================================================================
# 2. READ APPRAISE TABLES
# ==============================================================================

name_re <- "^(\\d)\\d_(sed|sea)_(allreads|10Mreads)_(mags|tkrs226genomes)\\.(species|genus)_appraise\\.tsv$"

read_appraise <- function(path) {
  parts <- str_match(basename(path), name_re)
  if (is.na(parts[1])) {
    message("Skipped (name not recognised): ", basename(path))
    return(NULL)
  }
  habitat <- c("0" = "Sediment", "1" = "Seawater")[[parts[2]]]
  if (c(sed = "Sediment", sea = "Seawater")[[parts[3]]] != habitat) {
    warning(basename(path), " is named '", parts[3], "' but numbered as ", habitat,
            "; using ", habitat, call. = FALSE)
  }
  read_tsv(path, comment = "#", show_col_types = FALSE) %>%
    mutate(file      = basename(path),
           Habitat   = habitat,
           Reads     = c(allreads = "All reads", `10Mreads` = "10 M read pairs")[[parts[4]]],
           Reference = c(mags = "MAGs", tkrs226genomes = "GTDB R10-RS226")[[parts[5]]],
           Rank      = str_to_title(parts[6]),
           Domain    = str_remove(domain, "^d__"))
}

appraise <- list.files(results_dir, pattern = "_appraise\\.tsv$", full.names = TRUE) %>%
  map(read_appraise) %>%
  bind_rows() %>%
  mutate(Habitat   = factor(Habitat, levels = names(habitat_cols)),
         Reads     = factor(Reads, levels = c("All reads", "10 M read pairs")),
         Reference = factor(Reference, levels = c("MAGs", "GTDB R10-RS226")),
         Group     = factor(paste(Domain, Rank, sep = "\n"),
                            levels = c("Archaea\nSpecies", "Archaea\nGenus",
                                       "Bacteria\nSpecies", "Bacteria\nGenus")))

pooled     <- filter(appraise, sample == "total")
per_sample <- filter(appraise, !sample %in% c("total", "average"))

write_tsv(per_sample, file.path(out_dir, "appraise_percent_binned.tsv"))

per_sample %>%
  group_by(Reads, Reference, Habitat, Domain, Rank) %>%
  summarise(samples = n(), mean = mean(percent_binned), median = median(percent_binned),
            .groups = "drop") %>%
  print(n = Inf)


# ==============================================================================
# 3. SCATTER PLOT
# ==============================================================================

dodge <- position_dodge(width = 0.7)

p <- ggplot(per_sample, aes(Group, percent_binned, colour = Habitat)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7, seed = 1),
             size = 0.8, alpha = 0.5) +
  geom_errorbar(data = pooled, aes(ymin = percent_binned, ymax = percent_binned, group = Habitat),
                position = dodge, width = 0.3, colour = "black", linewidth = 0.8) +
  facet_grid(Reads ~ Reference) +
  scale_colour_manual(values = habitat_cols) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "mOTU coverage assigned (%)", colour = NULL,
       caption = "Black bars: pooled over all samples") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "appraise_percent_binned.svg"), p, width = 9, height = 7)
ggsave(file.path(out_dir, "appraise_percent_binned.png"), p, width = 9, height = 7, dpi = 300)
