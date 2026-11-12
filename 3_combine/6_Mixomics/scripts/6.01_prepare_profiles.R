# ==============================================================================
# 6.01_prepare_profiles.R
#
# Community profiles for the reef zoning analyses (6.02 to 6.04), for each
# habitat (sediment, seawater) and profile type (MAG, mOTU).
#   1. Setup
#   2. MAG-based profiles (CoverM read counts; prokaryotic MAGs)
#   3. mOTU-based profiles (SingleM; prokaryotes)
#   4. Analysis inputs: samples with complete design variables, counts, CLR
#
# MAG profiles: CoverM read counts per MAG (95% ANI, 75% aligned, 95% identity;
# see 4_MAG_recovery). mOTU profiles: SingleM relative coverage scaled to
# estimated prokaryotic reads per sample (metagenome size x read fraction from
# SingleM microbial_fraction; see 1_SingleM).
#
# Input : CoverM workbook (sheets <HAB>_ENV, <HAB>_95_75_COUNT, <HAB>_TAX);
#         SingleM workbook (sheets seased_env_RS226, seased_cov_RS226,
#         seased_tax_RS226)
# Output: 0_profiles/profile_<HAB>_<METHOD>.rds: phyloseq object of all samples
#         (also used in 7.01); 0_profiles/zoning_input_<HAB>_<METHOD>.rds: list
#         with the phyloseq object, counts (samples x features), CLR values,
#         and sample metadata of samples with complete design variables;
#         HAB = SED or SEA, METHOD = MAG or SINGLEM
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(openxlsx)
  library(tidyverse)
})

mag_xlsx     <- "<path/to/1_coverm_ANI95_ALI15_MAP75_id_95_WORKING.xlsx>"
singlem_xlsx <- "<path/to/IMOS_10Mseqs_SingleM_FINAL.xlsx>"
out_dir      <- "<path/to/6_Mixomics>/0_profiles"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Design variables (column names in the metadata sheets)
zone_col    <- "FishingStatus"     # reef zoning status
zone_levels <- c("Closed", "Open") # Closed = no-take marine reserve (NTMR), Open = fished
sector_col  <- "GBRSector"
reef_col    <- "ReefName"
trip_col    <- "TripDescription"
design_cols <- c(zone_col, sector_col, reef_col, trip_col, "TripSeason")

habitats <- c(SED = "Sediment", SEA = "Seawater")


# ==============================================================================
# 2. MAG-BASED PROFILES
# ==============================================================================

read_mag_ps <- function(hab) {
  meta <- read.xlsx(mag_xlsx, sheet = paste0(hab, "_ENV"),        rowNames = TRUE)
  otu  <- read.xlsx(mag_xlsx, sheet = paste0(hab, "_95_75_COUNT"), rowNames = TRUE)
  tax  <- read.xlsx(mag_xlsx, sheet = paste0(hab, "_TAX"),         rowNames = TRUE) %>% as.matrix()

  ps <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE),
                 sample_data(meta), tax_table(tax))

  subset_taxa(ps, Domain %in% c("Bacteria", "Archaea")) %>%
    tax_fix(min_length = 1) %>%
    phyloseq_validate()
}

mag_ps <- lapply(names(habitats), read_mag_ps) %>% setNames(names(habitats))


# ==============================================================================
# 3. mOTU-BASED PROFILES
# ==============================================================================

sm_meta <- read.xlsx(singlem_xlsx, sheet = "seased_env_RS226", rowNames = TRUE)
sm_otu  <- read.xlsx(singlem_xlsx, sheet = "seased_cov_RS226", rowNames = TRUE)
sm_tax  <- read.xlsx(singlem_xlsx, sheet = "seased_tax_RS226", rowNames = TRUE) %>% as.matrix()

# Relative coverage scaled to estimated prokaryotic reads per sample
sm <- sm_meta[colnames(sm_otu), ]
stopifnot(!anyNA(sm$metagenome_size), !anyNA(sm$read_fraction))
prok_reads <- sm$metagenome_size * sm$read_fraction / 100

sm_counts <- as.matrix(sm_otu) %>%
  sweep(2, colSums(.), `/`) %>%
  sweep(2, prok_reads, `*`) %>%
  round()

singlem_ps_all <- phyloseq(otu_table(sm_counts, taxa_are_rows = TRUE),
                           sample_data(sm_meta), tax_table(sm_tax)) %>%
  tax_fix(min_length = 1) %>%
  phyloseq_validate() %>%
  subset_taxa(Domain %in% c("Bacteria", "Archaea"))

singlem_ps <- lapply(habitats, function(h) ps_filter(singlem_ps_all, Source == h))


# ==============================================================================
# 4. ANALYSIS INPUTS
# ==============================================================================

clr_counts <- function(counts, pseudocount = 1) {  # samples x features
  lx <- log(counts + pseudocount)
  lx - rowMeans(lx)
}

prepare_input <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  stopifnot(all(design_cols %in% names(meta)))

  # Samples with complete design variables
  ps   <- prune_samples(complete.cases(meta[, design_cols]), ps)
  meta <- as(sample_data(ps), "data.frame")
  meta[[zone_col]] <- factor(meta[[zone_col]], levels = zone_levels)
  stopifnot(!anyNA(meta[[zone_col]]))  # fails if other zone labels are present

  # Zoning is a reef-level attribute: one zone per reef
  stopifnot(all(tapply(as.character(meta[[zone_col]]), meta[[reef_col]],
                       function(z) length(unique(z))) == 1))

  counts <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) counts <- t(counts)
  counts <- counts[, colSums(counts) > 0, drop = FALSE]
  clr    <- clr_counts(counts)
  stopifnot(identical(rownames(clr), rownames(meta)))

  list(ps = ps, counts = counts, clr = clr, meta = meta)
}

profiles <- list(MAG = mag_ps, SINGLEM = singlem_ps)

for (method in names(profiles)) {
  for (hab in names(habitats)) {
    saveRDS(profiles[[method]][[hab]],
            file.path(out_dir, paste0("profile_", hab, "_", method, ".rds")))
    inp <- prepare_input(profiles[[method]][[hab]])
    cat(sprintf("%s %s: %d samples, %d features, %d reefs\n", hab, method,
                nrow(inp$clr), ncol(inp$clr), length(unique(inp$meta[[reef_col]]))))
    print(table(inp$meta[[zone_col]], inp$meta[[sector_col]]))
    saveRDS(inp, file.path(out_dir, paste0("zoning_input_", hab, "_", method, ".rds")))
  }
}
