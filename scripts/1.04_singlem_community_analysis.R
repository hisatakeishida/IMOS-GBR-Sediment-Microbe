# ==============================================================================
# 1.4_singlem_community_analysis.R
#
# mOTU-based (SingleM) community analysis of GBR sediment and seawater.
#   1. Setup: packages, paths, palettes, helpers
#   2. Load SingleM profiles and scale to estimated prokaryotic read counts
#   3. Sediment vs seawater: PCA (Fig. 1d), Shannon diversity, PERMANOVA
#   4. Within-habitat structure: PCA (Fig. 1g-h) and sequential PERMANOVA
#   5. Within-habitat structure: reef-level PERMANOVA
#   6. QC: sequencing facility check (sediment)
#
# Input : IMOS_10Mseqs_SingleM_FINAL.xlsx with sheets
#           seased_cov_RS226 (coverage matrix from 1.3_singlem_prof_coverage.py,
#                             taxa x samples),
#           seased_tax_RS226 (taxonomy), seased_env_RS226 (sample metadata,
#                             including microbial_fraction output)
# Output: singlem_prok_SED.rds, singlem_prok_SEA.rds (used by downstream scripts),
#         figures and PERMANOVA tables in out_dir
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

# ---- Packages ----------------------------------------------------------------
# install.packages(c("BiocManager", "openxlsx", "tidyverse", "patchwork",
#                    "data.table", "vegan"))
# BiocManager::install("phyloseq")
# install.packages("microViz",
#   repos = c("https://david-barnett.r-universe.dev", getOption("repos")))

library(data.table)   # used by merge_samples2()
library(phyloseq)
library(microViz)
library(vegan)
library(openxlsx)
library(tidyverse)
library(patchwork)

select  <- dplyr::select
tax_agg <- microViz::tax_agg

# ---- Paths and switches ------------------------------------------------------
singlem_xlsx <- "<path/to/IMOS_10Mseqs_SingleM_FINAL.xlsx>"
out_dir      <- "<path/to/1_SingleM_ANALYSIS>"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_outputs <- TRUE   # write figures and tables to out_dir
n_perm       <- 9999   # PERMANOVA permutations
rank_used    <- "Family"

save_fig <- function(plot, file, width = 7.5, height = 7.5) {
  if (save_outputs) ggsave(file.path(out_dir, file), plot = plot, width = width, height = height)
  invisible(plot)
}
save_tab <- function(tab, file) {
  if (save_outputs) write.xlsx(as.data.frame(tab), file.path(out_dir, file), rowNames = TRUE)
  invisible(tab)
}

# ---- Metadata column names ---------------------------------------------------
trip_col    <- "TripDescription"   # sampling trip (Trip 1 to Trip 4)
zone_col    <- "FishingStatus"     # reef zoning status: Closed (NTMR) or Open (fished)
sector_col  <- "GBRSector"
reef_col    <- "ReefName"
habitat_col <- "Source"            # Sediment or Seawater

# ---- Palettes ----------------------------------------------------------------
trip_col_vec <- c(
  "Trip 1 | Nov-Dec 2019"  = "#f2b09e",
  "Trip 2 | January 2020"  = "#5b5b5b",
  "Trip 3 | February 2020" = "#f3c789",
  "Trip 4 | July 2020"     = "slateblue")

source_cols <- c(Seawater = "#1f78b4", Sediment = "#33a02c")
seq_cols    <- c(ACE = "red", MICROBA = "blue", QUT = "green")

# ---- Helper: ordination label style ------------------------------------------
taxa_label_style <- function() {
  tax_lab_style(type = "label", max_angle = 90, aspect_ratio = 1,
                size = 3, alpha = 0.9, fontface = "bold",
                label.r = unit(0, "mm"))
}

# ---- Helper: within-habitat PCA coloured and shaped by sampling trip ---------
plot_trip_pca <- function(ps, rank = rank_used, trans = "clr", aspect = 1) {
  ps %>%
    tax_agg(rank) %>%
    tax_transform(trans) %>%
    ord_calc(method = "PCA") %>%
    ord_plot(size = 4, alpha = 1, fill = trip_col, shape = trip_col, plot_taxa = FALSE) +
    stat_chull(aes(group = .data[[trip_col]], colour = .data[[trip_col]]),
               geom = "polygon", alpha = 0) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_fill_manual(values = trip_col_vec) +
    scale_colour_manual(values = trip_col_vec) +
    theme_bw() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
          aspect.ratio = aspect)
}

# ---- Helper: PERMANOVA R2 bar plot -------------------------------------------
plot_permanova_r2 <- function(aov_tab, highlight = zone_col) {
  as.data.frame(aov_tab) %>%
    rownames_to_column("Term") %>%
    filter(Term != "Total") %>%
    mutate(Percent = R2 * 100, Term = fct_reorder(Term, Percent)) %>%
    ggplot(aes(Percent, Term, fill = Term == highlight)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", Percent)), hjust = -0.1, size = 3.5) +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey70"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Variation explained (%)", y = NULL) +
    theme_classic()
}

# ---- Helper: Aitchison distance at a given rank ------------------------------
aitchison_dist <- function(ps, rank = rank_used) {
  ps %>%
    tax_agg(rank) %>%
    tax_transform("identity") %>%
    dist_calc("aitchison") %>%
    dist_get()
}

# ---- Helper: merge_samples2() ------------------------------------------------
# Adapted from speedyseq (M. McLaren). Like phyloseq::merge_samples(), but
# handles categorical sample variables: abundances are combined with `fun_otu`;
# sample variables collapse to their unique value, otherwise NA.
setGeneric("merge_samples2",
           function(x, group, fun_otu = sum, funs = list(), reorder = FALSE)
             standardGeneric("merge_samples2"))

setMethod("merge_samples2", signature("phyloseq"),
          function(x, group, fun_otu = sum, funs = list(), reorder = FALSE) {
            if (length(group) == 1) {
              stopifnot(group %in% sample_variables(x))
              group <- sample_data(x)[[group]]
            } else {
              stopifnot(identical(length(group), nsamples(x)))
            }
            if (anyNA(group)) {
              warning("`group` has missing values; corresponding samples will be dropped")
              x <- prune_samples(!is.na(group), x)
              group <- group[!is.na(group)]
            }
            otu.merged <- merge_samples2(otu_table(x), group, fun_otu = fun_otu, reorder = reorder)
            sam.merged <- if (!is.null(access(x, "sam_data")))
              merge_samples2(sample_data(x), group, funs = funs) else NULL
            phyloseq(otu.merged, sam.merged,
                     access(x, "tax_table"), access(x, "phy_tree"), access(x, "refseq"))
          })

setMethod("merge_samples2", signature("otu_table"),
          function(x, group, fun_otu = sum, reorder = FALSE) {
            stopifnot(identical(length(group), nsamples(x)))
            needs_flip <- taxa_are_rows(x)
            if (needs_flip) x <- t(x)
            if (anyNA(group)) {
              warning("`group` has missing values; corresponding samples will be dropped")
              x <- x[!is.na(group), ]
              group <- group[!is.na(group)]
            }
            if (identical(fun_otu, sum)) {
              x.merged <- rowsum(x, group, reorder = reorder)
            } else {
              stopifnot(!".group" %in% colnames(x))
              f <- purrr::as_mapper(fun_otu)
              x <- x %>% as("matrix") %>% data.table::as.data.table() %>% cbind(.group = group)
              x.merged <- if (reorder) x[, lapply(.SD, f), keyby = .(.group)]
                          else         x[, lapply(.SD, f), by = .(.group)]
              rns <- x.merged$.group
              x.merged[, .group := NULL]
              x.merged <- as(x.merged, "matrix")
              rownames(x.merged) <- rns
            }
            x.merged <- otu_table(x.merged, taxa_are_rows = FALSE)
            if (needs_flip) x.merged <- t(x.merged)
            x.merged
          })

setMethod("merge_samples2", signature("sample_data"),
          function(x, group, funs = list(), reorder = FALSE) {
            if (length(group) == 1) {
              stopifnot(group %in% sample_variables(x))
              group <- x[[group]]
            } else {
              stopifnot(identical(length(group), nsamples(x)))
            }
            if (anyNA(group)) {
              warning("`group` has missing values; corresponding samples will be dropped")
              x <- x[!is.na(group), ]
              group <- group[!is.na(group)]
            }
            var_in_funs <- names(x) %>% rlang::set_names(. %in% names(funs), .)
            funs <- purrr::map2(var_in_funs, names(var_in_funs),
                                ~ if (.x) purrr::as_mapper(funs[[.y]]) else unique_or_na)
            new_sample_names <- group %>% unique() %>% sort() %>% as.character()
            x.merged <- purrr::map2(x, funs, ~ merge_groups(.x, group = group, f = .y)) %>%
              data.frame() %>%
              vctrs::vec_set_names(new_sample_names)
            if (!reorder) {
              initial_order <- group %>% unique() %>% as.character()
              x.merged <- x.merged[initial_order, , drop = FALSE]
            }
            sample_data_stable(x.merged)
          })

unique_or_na <- function(x) UseMethod("unique_or_na")
unique_or_na.default <- function(x) {
  if (length(unique(x)) == 1) x[[1]]
  else if (is.atomic(x)) as(NA, typeof(x))
  else NA
}
unique_or_na.factor <- function(x) {
  if (length(unique(x)) == 1) x[[1]]
  else factor(NA, levels = levels(x), ordered = is.ordered(x))
}
merge_groups <- function(x, group, f = unique_or_na) {
  f <- purrr::as_mapper(f)
  split(x, group) %>%
    purrr::map(f) %>%
    { vctrs::vec_c(!!!., .name_spec = rlang::zap()) }
}
# Keeps row names as-is (phyloseq::sample_data() resets them to "1", "2", ...)
sample_data_stable <- function(object) {
  stopifnot(identical(class(object), "data.frame"))
  object <- phyloseq:::reconcile_categories(object)
  new("sample_data", object)
}


# ==============================================================================
# 2. LOAD SINGLEM PROFILES AND SCALE TO PROKARYOTIC READ COUNTS
# ==============================================================================

singlem_meta <- read.xlsx(singlem_xlsx, sheet = "seased_env_RS226", rowNames = TRUE)
singlem_otu  <- read.xlsx(singlem_xlsx, sheet = "seased_cov_RS226", rowNames = TRUE)
singlem_tax  <- read.xlsx(singlem_xlsx, sheet = "seased_tax_RS226", rowNames = TRUE) %>%
  as.matrix()

# Excel date serial to Date (original serial kept)
singlem_meta$Sampling_Date <- singlem_meta$Date
singlem_meta$Date <- as.Date(singlem_meta$Date, origin = "1899-12-30")

# Scale relative coverage to estimated prokaryotic read counts using the
# SingleM microbial_fraction output (metagenome_size, read_fraction), and add
# one taxon holding the estimated non-prokaryotic reads per sample
sm <- singlem_meta[colnames(singlem_otu), ]
stopifnot(!anyNA(sm$metagenome_size))
prok_reads    <- sm$metagenome_size * sm$read_fraction    / 100
nonprok_reads <- sm$metagenome_size * sm$nonprok_fraction / 100

singlem_counts <- as.matrix(singlem_otu) %>%
  sweep(2, colSums(.), `/`) %>%
  sweep(2, prok_reads, `*`) %>%
  rbind(SingleM_nonprok = nonprok_reads) %>%
  round()

nonprok_tax <- setNames(rep("", ncol(singlem_tax)), colnames(singlem_tax))
nonprok_tax[1] <- "SingleM non-bac/arc"
singlem_tax <- rbind(singlem_tax, SingleM_nonprok = nonprok_tax)

singlem_ps <- phyloseq(
  otu_table(singlem_counts, taxa_are_rows = TRUE),
  sample_data(singlem_meta),
  tax_table(singlem_tax)
) %>%
  tax_fix(min_length = 1) %>%
  phyloseq_validate()

# Prokaryotes only, all samples and per habitat
singlem_prok <- subset_taxa(singlem_ps, Domain %in% c("Bacteria", "Archaea"))
singlem_sed  <- ps_filter(singlem_prok, Source == "Sediment")
singlem_sea  <- ps_filter(singlem_prok, Source == "Seawater")

message("Samples: ", nsamples(singlem_sed), " sediment, ", nsamples(singlem_sea), " seawater")

saveRDS(singlem_sed, file.path(out_dir, "singlem_prok_SED.rds"))
saveRDS(singlem_sea, file.path(out_dir, "singlem_prok_SEA.rds"))


# ==============================================================================
# 3. SEDIMENT VS SEAWATER
# ==============================================================================

# ---- PCA of family-level CLR profiles (Fig. 1d) ------------------------------
p_habitat_pca <- singlem_prok %>%
  tax_transform(trans = "clr", rank = rank_used) %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    size = 6, fill = habitat_col, shape = trip_col, plot_taxa = 10:1,
    tax_vec_style_all = vec_constraint(linewidth = 0, arrow = NULL),
    tax_lab_style = taxa_label_style()
  ) +
  stat_ellipse(aes(colour = .data[[habitat_col]]), linewidth = 0.2) +
  scale_fill_manual(values = source_cols) +
  scale_colour_manual(values = source_cols) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        legend.position = "none") +
  coord_fixed(xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75))
p_habitat_pca
save_fig(p_habitat_pca, "Fig1d_PCA_habitat.pdf")

# ---- Alpha diversity: family-level Shannon index -----------------------------
alpha_tbl <- singlem_prok %>%
  ps_calc_diversity(rank = rank_used, index = "shannon") %>%
  samdat_tbl()

alpha_summary <- alpha_tbl %>%
  group_by(.data[[habitat_col]]) %>%
  summarise(mean_Shannon = mean(shannon_Family, na.rm = TRUE),
            sd_Shannon   = sd(shannon_Family,   na.rm = TRUE))
alpha_summary
wilcox.test(reformulate(habitat_col, "shannon_Family"), data = alpha_tbl)

p_alpha <- ggplot(alpha_tbl, aes(y = .data[[habitat_col]], x = shannon_Family,
                                 fill = .data[[habitat_col]])) +
  geom_boxplot(width = 0.3) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  scale_fill_manual(values = source_cols) +
  labs(x = "Shannon index (family level)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none")
save_fig(p_alpha, "Shannon_habitat.pdf", height = 3)

# ---- Beta diversity: PERMANOVA and dispersion by habitat ---------------------
habitat_dist <- aitchison_dist(singlem_prok)
habitat_meta <- as(sample_data(singlem_prok), "data.frame")

set.seed(123)
permanova_habitat <- adonis2(reformulate(habitat_col, "habitat_dist"),
                             data = habitat_meta, permutations = n_perm)
permanova_habitat
save_tab(permanova_habitat, "PERMANOVA_habitat.xlsx")

habitat_disp <- betadisper(habitat_dist, habitat_meta[[habitat_col]])
set.seed(123)
permutest(habitat_disp, permutations = n_perm)


# ==============================================================================
# 4. WITHIN-HABITAT STRUCTURE: PCA AND SEQUENTIAL PERMANOVA
# ==============================================================================

# Sequential PERMANOVA with zoning status entered first:
#   Aitchison distance ~ zoning status + sampling trip + GBR sector + reef site
# Zoning status is a reef-level attribute, so its sample-level p-value is
# optimistic; see the reef-level PERMANOVA in section 5.
run_within_habitat <- function(ps, label) {

  # PCA coloured by sampling trip (Fig. 1g-h)
  p_pca <- plot_trip_pca(ps)
  save_fig(p_pca, paste0("Fig1gh_PCA_trip_", label, ".pdf"))

  d    <- aitchison_dist(ps)
  meta <- as(sample_data(ps), "data.frame")

  set.seed(123)
  perm_full <- adonis2(reformulate(c(zone_col, trip_col, sector_col, reef_col), "d"),
                       data = meta, permutations = n_perm, by = "terms")
  save_tab(perm_full, paste0("PERMANOVA_sequential_", label, ".xlsx"))
  save_fig(plot_permanova_r2(perm_full), paste0("PERMANOVA_sequential_", label, ".pdf"),
           height = 4)

  list(pca = p_pca, permanova = perm_full)
}

within_sed <- run_within_habitat(singlem_sed, "SED")
within_sea <- run_within_habitat(singlem_sea, "SEA")

within_sed$permanova
within_sea$permanova
within_sed$pca | within_sea$pca


# ==============================================================================
# 5. WITHIN-HABITAT STRUCTURE: REEF-LEVEL PERMANOVA
# ==============================================================================

# Replicates are summed per reef so that each reef contributes one observation,
# avoiding pseudo-replication for reef-level factors (zoning status, sector, trip).
run_reef_level <- function(ps, label) {
  ps_reef <- merge_samples2(ps, reef_col, fun_otu = sum)
  d_reef  <- aitchison_dist(ps_reef)
  meta    <- as(sample_data(ps_reef), "data.frame")

  stopifnot(identical(rownames(meta), labels(d_reef)))                    # rows match distances
  stopifnot(!anyNA(meta[, c(zone_col, trip_col, sector_col)]))            # one value per reef

  set.seed(123)
  perm_reef <- adonis2(reformulate(c(zone_col, trip_col, sector_col), "d_reef"),
                       data = meta, permutations = n_perm, by = "terms")
  save_tab(perm_reef, paste0("PERMANOVA_reeflevel_", label, ".xlsx"))
  perm_reef
}

reef_sed <- run_reef_level(singlem_sed, "SED")
reef_sea <- run_reef_level(singlem_sea, "SEA")
reef_sed
reef_sea


# ==============================================================================
# 6. QC: SEQUENCING FACILITY CHECK (SEDIMENT)
# ==============================================================================

# Sediment samples were sequenced at three facilities; check that samples do not
# cluster by facility
p_seq_check <- singlem_sed %>%
  tax_transform(trans = "clr", rank = rank_used) %>%
  ord_calc(method = "PCA") %>%
  ord_plot(size = 4, colour = "SEQUENCING", fill = "SEQUENCING", shape = trip_col,
           plot_taxa = FALSE) +
  stat_ellipse(aes(colour = SEQUENCING), linewidth = 0.2) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_colour_manual(values = seq_cols) +
  scale_fill_manual(values = seq_cols) +
  theme_bw()
save_fig(p_seq_check, "QC_PCA_sequencing_facility_SED.pdf")


# ==============================================================================
sessionInfo()