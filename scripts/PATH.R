suppressPackageStartupMessages({
  library(openxlsx)
  library(patchwork)
  library(tidyverse)
})

infile_bs <-"~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/imosmanuscript/ARCHIVE/BGCs_GCFs_IMOS_OMD.xlsx"

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


df_binary <- read.xlsx("~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/3_MAG_GTDB/dereplication_95_ANI.xlsx", sheet="Summary_ANI95_AF15_OMD_RMD_GEMs") %>%
  column_to_rownames("Rep") %>%
  mutate(across(everything(), ~ as.integer(.x > 0)))

out_file     <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results/Fig2d_novelty_upset.svg"

magstats <- read.xlsx("~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/2_MAG_ANALYSIS/Sea__Sed_indiass_coass.xlsx",sheet = "derep", rowNames=TRUE)

stats_path  <- "~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/2_MAG_ANALYSIS/Sea__Sed_indiass_coass.xlsx"
stats_sheet <- "derep"
out_dir <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results"
setwd(out_dir)

taxonomy_path  <- "~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/1_IMOS/imosmanuscript/3_Sup_Table/ARCH/table_S5.xlsx"
taxonomy_sheet <- "curatedtaxonomy_Sed1201_Sea876"
out_file       <- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results/Fig2b_phylum_counts.svg"
out_file2<- "/Users/s4573340/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/IMOS-GBR-Sediment-Microbe/4_results/Fig2d_novelty_upset.all.svg"

