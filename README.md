# Metagenomics analysis of sediment microbes across the Great Barrier Reef

## This repository contains codes and data associated with the paper: 


### Overview 
This project implemented gene-centric and genome-centric analyses to investigate taxonomic and metabolic diversity of surface sediment microbiomes of the Great Barrier Reef. This project is part of the Great Barrier Reef Microbial Genomics Database sub-Facility within the Marine Microbiome Initiative Facility of Australia's Integrated Marine Observing System (IMOS). More details at (https://imos.org.au/facilities/marinemicrobiomeinitiative/gbr-microbial-genomic-database). 

### Table of Contents

#### Sample, Sequencing information, and QC
- **[1. Sample, Sequencing information and QC](0_metadata/0_sample_data_qc.md)**

#### Seawater samples analysed along sediment samples
- **[1. Read and MAG information](4_seawater/4_seawater.md)**

#### Gene centric analysis
- **[1. mOTU-based analyses with SingleM](1_gene_centric/1_SingleM/1_SingleM.md)**

- **[2. non-redundant protein analyses](1_gene_centric/2_Protein/2_Protein.md)**

- **[3. k-mer-based analyses with Sylph](1_gene_centric/3_Sylph/3_Sylph.md)**

#### Genome centric analysis
- **[1. MAG recovery, QC, Dereplication, Taxonomy, Novelty](2_genome_centric/1_genome_recovery/1_genome_recovery.md)**

- **[2. MAG annotation](2_genome_centric/2_genome_annot/2_genome_annot.md)**

- **[3. MAG abundance](2_genome_centric/3_genome_abundance/3_genome_abundance.md)**

- **[4. MAG Mixomics analysis](2_genome_centric/4_Mixomics/4_Mixomics.md)**

#### Integration 
- **[1. SingleM vs CoverM abudance comparison](3_SingleM_VS_CoverM/3_SingleM_VS_CoverM.md)**

#### Rscripts
**[MASH visualization](5_Rscripts/1_MASH.R)**
**[SingleM visualization](5_Rscripts/2_SingleM.R)**
**[RPKM_KO visualization](5_Rscripts/3_RPKM_calculation_KO.R)**
**[RPKM_CAZyme visualization](5_Rscripts/4_RPKM_calculation_CAZyme.R)**

### Authors and Affiliations

### Questions or issues 
Please contact Hisatake Ishida at hisatake.ishida@uq.net.au
