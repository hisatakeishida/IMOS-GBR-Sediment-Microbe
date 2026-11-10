# ==============================================================================
# 4.16_phylum_counts.R
#
# Number of sediment and seawater MAGs in the 20 most represented phyla
# (Pseudomonadota split into classes), as a diverging bar plot: sediment to
# the left, seawater to the right (Fig. 2b).
#
# Input : curated taxonomy of the 1,201 sediment and 876 seawater MAGs, one row
#         per MAG, with columns SOURCE (Sediment or Seawater) and
#         Phylum_curated (from 4.9)
# Output: Fig. 2b (svg)
# ==============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
})

taxonomy_path  <- "<path/to/MAG_taxonomy_table.xlsx>"
taxonomy_sheet <- "curatedtaxonomy_Sed1201_Sea876"
out_file       <- "<path/to/4_MAG_recovery/7_taxonomy>/Fig2b_phylum_counts.svg"

n_top        <- 20
habitat_cols <- c(Seawater = "#1f78b4", Sediment = "#33a02c")

mag_taxonomy <- read.xlsx(taxonomy_path, sheet = taxonomy_sheet)
count(mag_taxonomy, SOURCE)


# ==============================================================================
# COUNTS PER PHYLUM (TOP 20, OTHERS COMBINED)
# ==============================================================================

top_phyla <- mag_taxonomy %>%
  count(Phylum_curated, sort = TRUE) %>%
  slice_head(n = n_top) %>%
  pull(Phylum_curated)

counts <- mag_taxonomy %>%
  mutate(Phylum_plot = if_else(Phylum_curated %in% top_phyla, Phylum_curated, "Others")) %>%
  count(Phylum_plot, SOURCE) %>%
  mutate(n_signed = if_else(SOURCE == "Sediment", -n, n))   # sediment plotted to the left

# Order by total number of MAGs, with Others at the bottom
phylum_order <- counts %>%
  group_by(Phylum_plot) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(total) %>%
  pull(Phylum_plot)
phylum_order <- c("Others", setdiff(phylum_order, "Others"))
counts$Phylum_plot <- factor(counts$Phylum_plot, levels = phylum_order)


# ==============================================================================
# PLOT
# ==============================================================================

p_phyla <- ggplot(counts, aes(x = Phylum_plot, y = n_signed, fill = SOURCE)) +
  geom_col(position = "identity", width = 0.9, colour = "black") +
  geom_text(aes(label = n), colour = "white", size = 3.5,
            position = position_stack(vjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = habitat_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text       = element_text(size = 8))
p_phyla

ggsave(out_file, p_phyla, width = 3.23, height = 4)

sessionInfo()