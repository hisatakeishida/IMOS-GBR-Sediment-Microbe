# ==============================================================================
# 6.10_singlem_rda.R
#
# Supporting analysis (not reported): environmental and design drivers of
# within-habitat community composition (SingleM mOTUs, family level).
#   1. Setup and input
#   2. Pairwise PERMANOVA and dispersion tests between sampling trips
#   3. Response matrix and predictors
#   4. RDA and partial RDA of each term given the others
#   5. Constrained ordination plot
#
# Usage : Rscript 3_combine/6_Mixomics/scripts/6.10_singlem_rda.R <HAB>
# Input : profile_<HAB>_SINGLEM.rds (6.01)
# Output: tables (xlsx) and figures (png)
# ==============================================================================


# ==============================================================================
# 1. SETUP AND INPUT
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(vegan)
  library(openxlsx)
  library(tidyverse)
})
select <- dplyr::select

args <- commandArgs(trailingOnly = TRUE)
hab  <- if (length(args) >= 1) args[1] else "SED"

in_dir  <- "<path/to/6_Mixomics>"
profile_dir <- file.path(in_dir, "0_profiles")   # output of 6.01
out_dir <- file.path(in_dir, "10_SingleM_RDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rda_rank  <- "Family"
trip_cols <- c("Trip 1" = "#f2b09e", "Trip 2" = "#5b5b5b", "Trip 3" = "#f3c789", "Trip 4" = "slateblue")

ps <- readRDS(file.path(profile_dir, paste0("profile_", hab, "_SINGLEM.rds")))


# ==============================================================================
# 2. PAIRWISE TRIP TESTS
# ==============================================================================
# Aitchison distance; PERMANOVA and dispersion (PERMDISP) for each pair of trips

pairwise_trip_tests <- function(ps, perms = 999) {
  trips <- ps %>% samdat_tbl() %>% pull(TripName) %>% na.omit() %>% unique()
  map_dfr(combn(trips, 2, simplify = FALSE), function(pair) {
    sub <- ps_filter(ps, TripName %in% pair)
    d   <- sub %>% tax_transform("identity", rank = "unique") %>% dist_calc("aitchison") %>% dist_get()
    md  <- as(sample_data(sub), "data.frame")
    tibble(pair        = paste(pair, collapse = " vs "),
           permanova_p = adonis2(d ~ TripName, data = md, permutations = perms)$`Pr(>F)`[1],
           permdisp_p  = permutest(betadisper(d, md$TripName), permutations = perms)$tab$`Pr(>F)`[1])
  }) %>%
    mutate(permanova_p_adj = p.adjust(permanova_p, method = "BH"))
}

set.seed(123)
trip_tests <- pairwise_trip_tests(ps)
print(trip_tests)


# ==============================================================================
# 3. RESPONSE MATRIX AND PREDICTORS
# ==============================================================================

clr_fam <- ps %>% tax_agg(rda_rank) %>% tax_transform("clr") %>% otu_get() %>% as("matrix")

env_vars <- c(
  "metagenome_size", "Efmean", "Latitude", "Longitude",
  "Seawater_Temperature", "Salinity", "Turbidity", "Fluorescence", "Chl_a", "Phaeophytin_a",
  "Particulate_Nitrogen__PN", "Particulate_Organic_Carbon__POC", "Particulate_Phosphorus__PP",
  "Dissolved_Organic_Carbon__DOC", "Phosphate__PO4", "Ammonium__NH4", "Nitrite__NO2",
  "Nitrate__NO3", "Silica__Si", "Total_Dissolved_Nitrogen__TDN",
  "Total_Dissolved_Phosphorus__TDP", "Total_Suspended_Solids__TSS"
)

# Continuous variables as z-scores ("_z"); categorical variables coded as numbers ("_ps")
meta_rda <- ps %>%
  samdat_tbl() %>%
  mutate(across(all_of(env_vars), ~ as.numeric(scale(.x)), .names = "{.col}_z"),
         TripSeason_ps    = if_else(TripSeason == "Winter", 1, 0),
         FishingStatus_ps = if_else(FishingStatus == "Open", 1, 0),
         TripName_ps      = recode(TripName, "Trip 1" = 1, "Trip 2" = 2, "Trip 3" = 3, "Trip 4" = 4),
         GBRSector_ps     = as.numeric(factor(GBRSector)),
         ReefName_ps      = as.numeric(factor(ReefName)))
stopifnot(identical(rownames(clr_fam), meta_rda$.sample_name))

# TripName_ps enters the RDA as a single numeric trend over trips 1 to 4
rda_terms <- c("Efmean_z", "TripSeason_ps", "TripName_ps", "Latitude_z")


# ==============================================================================
# 4. RDA AND PARTIAL RDA
# ==============================================================================

rda_full <- rda(reformulate(rda_terms, response = "clr_fam"), data = meta_rda, na.action = na.omit)
print(RsquareAdj(rda_full))
set.seed(123)
anova_full  <- anova.cca(rda_full, permutations = 999)
anova_terms <- anova.cca(rda_full, permutations = 999, by = "term")

partial_rda <- function(term) {
  cond <- sprintf("Condition(%s)", paste(setdiff(rda_terms, term), collapse = " + "))
  fit  <- rda(reformulate(c(term, cond), response = "clr_fam"), data = meta_rda, na.action = na.omit)
  r2   <- RsquareAdj(fit)
  set.seed(123)
  aov  <- anova.cca(fit, permutations = 999)
  tibble(variable = term, r2 = r2[["r.squared"]], adj_r2 = r2[["adj.r.squared"]],
         constrained_percent = 100 * fit$CCA$tot.chi / fit$tot.chi,
         F_val = aov$F[1], pval = aov$`Pr(>F)`[1])
}

partial_results <- map_dfr(rda_terms, partial_rda) %>% arrange(desc(adj_r2))
full_r2 <- RsquareAdj(rda_full)
fractions <- tibble(
  fraction   = c(partial_results$variable, "Confounded", "Residual"),
  adj_R2_pct = 100 * c(partial_results$adj_r2,
                       full_r2$adj.r.squared - sum(partial_results$adj_r2),
                       1 - full_r2$adj.r.squared)
)
print(partial_results)
print(fractions)

write.xlsx(list(pairwise_trips = trip_tests, rda_model = as.data.frame(anova_full),
                rda_terms = as.data.frame(anova_terms), partial_rda = partial_results,
                fractions = fractions),
           file.path(out_dir, paste0("SingleM_RDA_", hab, "_", rda_rank, ".xlsx")), rowNames = TRUE)


# ==============================================================================
# 5. CONSTRAINED ORDINATION PLOT
# ==============================================================================

ps_env <- ps
sample_data(ps_env) <- meta_rda %>% as.data.frame() %>% column_to_rownames(".sample_name") %>% sample_data()

p_rda <- ps_env %>%
  tax_agg(rda_rank) %>%
  tax_transform("clr") %>%
  ord_calc(constraints = rda_terms, method = "RDA", scale_cc = FALSE) %>%
  ord_plot(axes = c(1, 2), plot_taxa = 1:5, colour = "TripName", shape = "GBRSector", size = 2) +
  stat_chull(aes(group = TripName, colour = TripName), geom = "polygon", alpha = 0) +
  scale_shape_manual(values = c(21:25, 7, 12)) +
  scale_colour_manual(values = trip_cols) +
  theme_bw()
ggsave(file.path(out_dir, paste0("SingleM_RDA_", hab, "_", rda_rank, ".png")), p_rda, width = 8, height = 7)
