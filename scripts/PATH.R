
singlem_xlsx <- "~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/1_SingleM_ANALYSIS/IMOS_10Mseqs_SingleM_FINAL.xlsx"
out_dir <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results"
setwd(out_dir)


select <- dplyr::select

in_rds  <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/1_gene_centric/1_SingleM/1_community_profiles/singlem_prok_SEASED.rds"
out_dir <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results"
setwd(out_dir)



suppressPackageStartupMessages({
  library(tidyverse)
})

results_dir <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/1_gene_centric/1_SingleM/4_appraise_novelty/0_results"
out_dir <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results"
setwd(out_dir)

habitat_cols <- c(Sediment = "#33a02c", Seawater = "#1f78b4")
