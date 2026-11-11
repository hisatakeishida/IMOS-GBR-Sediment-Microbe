
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


suppressPackageStartupMessages({
  library(ComplexUpset)
  library(tidyverse)
})
dereplication_95_ANI <- read.xlsx("~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/3_MAG_GTDB/dereplication_95_ANI.xlsx", sheet="Summary_ANI95_AF15_OMD_RMD_GEMs")
summary_path <- "<path/to/4_MAG_recovery/9_novelty>/sed_sea_public_ANI95_AF15.cluster_summary.tsv"
out_file     <- "<path/to/4_MAG_recovery/9_novelty>/Fig2d_novelty_upset.svg"

