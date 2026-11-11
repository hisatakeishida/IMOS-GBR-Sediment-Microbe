# ==============================================================================
# 7.01_taxonomic_comparison.R
#
# Sediment community composition from the gene-centric (SingleM mOTUs) and
# genome-centric (CoverM, MAGs) approaches (Fig. S2, Table S7).
#   1. Setup and input
#   2. Mean relative abundance per sampling trip (tables)
#   3. log2 fold changes, CoverM / SingleM (Fig. S2a-d)
#   4. Rank comparison of the most abundant taxa (Fig. S2e-h)
#   5. Class-level composition per reef, by trip and approach (Fig. S2i)
#   6. Supporting: shared and unique families, CLR agreement, core families
#
# Comparisons use mean profiles per reef site (reads summed over the replicate
# samples of each reef, then converted to relative abundance). Taxa named by
# tax_fix() placeholders (unclassified at that rank; names containing a space)
# are excluded from sections 3 and 4.
#
# Input : profile_SED_MAG.rds and profile_SED_SINGLEM.rds (6.01)
# Output: tables (xlsx) and figures (svg)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(microbiome)   # core_members()
  library(ggrepel)
  library(ggVennDiagram)
  library(patchwork)
  library(openxlsx)
  library(tidyverse)
})

in_dir  <- "<path/to/6_Mixomics>"
out_dir <- "<path/to/7_Gene_Genome_analysis/1_taxonomy>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reef_col <- "ReefName"
trip_col <- "TripDescription"
ranks    <- setNames(c("Phylum", "Class", "Order", "Family"),
                     c("Phylum", "Class", "Order", "Family"))
name_mag <- "MAG_CoverM"
name_sm  <- "mOTU_SingleM"

mag_sed     <- readRDS(file.path(in_dir, "profile_SED_MAG.rds"))
singlem_sed <- readRDS(file.path(in_dir, "profile_SED_SINGLEM.rds"))

# merge_samples2() from speedyseq (M. McLaren) merges samples while keeping
# categorical sample variables: remotes::install_github("mikemc/speedyseq")
merge_samples2 <- speedyseq::merge_samples2

mag_reef     <- merge_samples2(mag_sed,     reef_col, fun_otu = sum)
singlem_reef <- merge_samples2(singlem_sed, reef_col, fun_otu = sum)


# ==============================================================================
# 2. MEAN RELATIVE ABUNDANCE PER SAMPLING TRIP
# ==============================================================================

rel_by_trip <- function(ps, rank) {
  ps %>%
    tax_transform("compositional", rank = rank) %>%
    ps_get() %>%
    merge_samples2(trip_col, fun_otu = mean) %>%
    otu_table() %>%
    { 100 * as(., "matrix") } %>%
    as.data.frame()
}

for (rank in c("Class", "Family")) {
  write.xlsx(list(CoverM  = rel_by_trip(mag_sed, rank),
                  SingleM = rel_by_trip(singlem_sed, rank)),
             file.path(out_dir, paste0(rank, "_relabund_by_trip.xlsx")),
             rowNames = TRUE)
}


# ==============================================================================
# 3. log2 FOLD CHANGES, CoverM / SingleM (Fig. S2a-d)
# ==============================================================================
# Mean relative abundance (%) across reefs per taxon, for taxa detected by both
# approaches; the 10 highest and 10 lowest log2(CoverM / SingleM) are shown.
# Negative values: taxon under-represented by CoverM relative to SingleM.

mean_by_taxon <- function(ps, rank, col) {
  ps %>%
    tax_transform("compositional", rank = rank) %>%
    psmelt() %>%
    group_by(.data[[rank]]) %>%
    summarise("{col}" := mean(Abundance), .groups = "drop")
}

log2fc_taxa <- function(ps_mag, ps_sm, rank, top_n = 10) {
  df <- full_join(mean_by_taxon(ps_mag, rank, "mean_mag"),
                  mean_by_taxon(ps_sm,  rank, "mean_sm"), by = rank) %>%
    replace_na(list(mean_mag = 0, mean_sm = 0)) %>%
    filter(mean_mag > 0, mean_sm > 0, !grepl(" ", .data[[rank]])) %>%
    mutate(mean_mag = 100 * mean_mag, mean_sm = 100 * mean_sm,
           log2FC = log2(mean_mag / mean_sm))

  top <- bind_rows(slice_max(df, log2FC, n = top_n, with_ties = FALSE),
                   slice_min(df, log2FC, n = top_n, with_ties = FALSE))
  top[[rank]] <- factor(top[[rank]], levels = top[[rank]][order(top$log2FC)])
  max_abs <- max(abs(top$log2FC))

  p <- ggplot(top, aes(.data[[rank]], log2FC, fill = mean_sm)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    scale_y_continuous(limits = c(-max_abs, max_abs)) +
    scale_fill_viridis_c(option = "plasma") +
    labs(x = rank, y = paste0("log2 fold change (", name_mag, " / ", name_sm, ")"),
         fill = paste0("Mean abundance (%)\n(", name_sm, ")")) +
    theme_minimal(base_size = 14)

  list(plot = p, data = df)
}

fc_res <- lapply(ranks, function(r) log2fc_taxa(mag_reef, singlem_reef, r))
fig_s2_ad <- wrap_plots(lapply(fc_res, `[[`, "plot"), nrow = 1)
ggsave(file.path(out_dir, "FigS2a-d_log2FC_CoverM_SingleM.svg"), fig_s2_ad, width = 22, height = 6)
write.xlsx(lapply(fc_res, `[[`, "data"), file.path(out_dir, "log2FC_CoverM_SingleM.xlsx"))


# ==============================================================================
# 4. RANK COMPARISON OF THE MOST ABUNDANT TAXA (Fig. S2e-h)
# ==============================================================================
# Taxa among the top 5 by mean relative abundance in either approach; taxa
# missing from one approach get that approach's lowest rank + 1. Spearman
# correlation of ranks over all taxa.

rank_taxa <- function(ps, rank) {
  ps %>%
    tax_transform("compositional", rank = rank) %>%
    psmelt() %>%
    filter(!grepl(" ", .data[[rank]])) %>%
    group_by(.data[[rank]]) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    mutate(rank_pos = row_number())
}

compare_ranks <- function(ps_mag, ps_sm, rank, top_n = 5) {
  r_mag <- rank_taxa(ps_mag, rank)
  r_sm  <- rank_taxa(ps_sm, rank)
  top_taxa <- union(head(r_mag[[rank]], top_n), head(r_sm[[rank]], top_n))

  dominance <- full_join(r_mag, r_sm, by = rank, suffix = c("_mag", "_sm")) %>%
    replace_na(list(mean_abund_mag = 0, mean_abund_sm = 0,
                    rank_pos_mag = max(r_mag$rank_pos) + 1,
                    rank_pos_sm  = max(r_sm$rank_pos) + 1))

  p <- dominance %>%
    filter(.data[[rank]] %in% top_taxa) %>%
    ggplot(aes(rank_pos_mag, rank_pos_sm, label = .data[[rank]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(colour = "steelblue", size = 3) +
    geom_text_repel(max.overlaps = 20, size = 6) +
    labs(x = paste("Rank in", name_mag), y = paste("Rank in", name_sm)) +
    theme_minimal(base_size = 14)

  list(rankplot = p, data = dominance,
       spearman = cor.test(dominance$rank_pos_mag, dominance$rank_pos_sm, method = "spearman"))
}

rank_res <- lapply(ranks, function(r) compare_ranks(mag_reef, singlem_reef, r))
fig_s2_eh <- wrap_plots(lapply(rank_res, `[[`, "rankplot"), nrow = 1)
ggsave(file.path(out_dir, "FigS2e-h_rank_CoverM_SingleM.svg"), fig_s2_eh, width = 22, height = 6)
print(lapply(rank_res, `[[`, "spearman"))


# ==============================================================================
# 5. CLASS-LEVEL COMPOSITION BY TRIP AND APPROACH (Fig. S2i)
# ==============================================================================

make_combined <- function(ps_mag, ps_sm, rank) {
  d1 <- ps_mutate(ps_mag, dataset = name_mag)
  d2 <- ps_mutate(ps_sm,  dataset = name_sm)
  sample_names(d1) <- paste0("MAG_",     sample_names(d1))
  sample_names(d2) <- paste0("SingleM_", sample_names(d2))
  merge_phyloseq(d1 %>% tax_agg(rank) %>% ps_get(),
                 d2 %>% tax_agg(rank) %>% ps_get())
}

combined_class <- make_combined(mag_reef, singlem_reef, "Class")

fig_s2_i <- combined_class %>%
  ps_arrange(.data[[trip_col]], .data[[reef_col]]) %>%
  comp_barplot("Class", merge_other = FALSE, n_taxa = 15,
               sample_order = "asis", label = reef_col) +
  coord_flip() +
  facet_wrap(vars(dataset, .data[[trip_col]]), scales = "free", nrow = 2)
ggsave(file.path(out_dir, "FigS2i_class_composition.svg"), fig_s2_i, width = 22, height = 11)


# ==============================================================================
# 6. SUPPORTING: SHARED FAMILIES, CLR AGREEMENT, CORE FAMILIES
# ==============================================================================

# Families detected by only one approach
clr_mag <- as(otu_get(tax_transform(mag_reef,     "clr", rank = "Family")), "matrix")
clr_sm  <- as(otu_get(tax_transform(singlem_reef, "clr", rank = "Family")), "matrix")

only_mag <- setdiff(colnames(clr_mag), colnames(clr_sm))
only_sm  <- setdiff(colnames(clr_sm), colnames(clr_mag))
cat(sprintf("Families: %d CoverM only, %d SingleM only, %d shared\n",
            length(only_mag), length(only_sm), length(intersect(colnames(clr_mag), colnames(clr_sm)))))

# Agreement of CLR values for shared families and reefs
common_taxa <- intersect(colnames(clr_mag), colnames(clr_sm))
common_reef <- intersect(rownames(clr_mag), rownames(clr_sm))
cat(sprintf("Pearson r of family CLR values (shared reefs and families): %.3f\n",
            cor(as.vector(clr_mag[common_reef, common_taxa]),
                as.vector(clr_sm[common_reef, common_taxa]))))

# Core families: present in >= 80% of samples
core_ps <- make_combined(mag_sed, singlem_sed, "Family") %>%
  tax_transform("compositional", rank = "Family") %>%
  ps_get()
core_list <- lapply(c(CoverM = name_mag, SingleM = name_sm), function(d) {
  core_members(ps_filter(core_ps, dataset == d), detection = 0, prevalence = 0.8)
})
print(lengths(core_list))
print(intersect(core_list$CoverM, core_list$SingleM))
ggsave(file.path(out_dir, "core_families_venn.svg"), ggVennDiagram(core_list), width = 6, height = 5)

write.xlsx(list(CoverM_only  = tibble(Family = only_mag),
                SingleM_only = tibble(Family = only_sm),
                core_CoverM  = tibble(Family = core_list$CoverM),
                core_SingleM = tibble(Family = core_list$SingleM)),
           file.path(out_dir, "families_shared_unique_core.xlsx"))
