# Information for samples, sequencing, and quality control 

## Table of Contents
1. [Sample and Sequencing information](#sample_info)
2. [Sequencing information](#seq_info)
3. [QC information](#qc_info)

## 1. Sample and Sequencing information <a name="sample_info"></a>
- 192 sediment samples (48 sites * 4 technical replicates) were collected across the GBR in 2019 and 2020
- Of 4 replicates, 3 shallow sequencing (2-10 Gbp) and 1 deep sequencing (30-70 Gbp)
- If sponge microbiome, expect 80-100 HQ MAGs from ~40 Gbp sequencing
- All sites associated with LTMP (long-term monitoring program, reef, algae cover, nutrient etc)
- Seidments were collected from the boat using GRAB-type approach
- Sample code
  -   Site_replicate#_SequenceSampleID

- We are not doing co-assembly 
-  27 Nanopore sequencing  
- Illumina Novaseq sequencing
- Sequenced at MICROBA, ACE, QUT
- Focus on 27 samples with hybrid assembly (expecting 10-15 MAGs per sample, thus 250-300 MAGs in total)

## 2. QC information <a name="qc_info"></a>
- porechop for Nanopore 
- trimmomatic for Illumina 
- use seqkit-stat for read information 

## 3. Initial diversity investigation using SingleM pipe <a name="qc_info"></a>
- Search for GTDB OTU data in reads 

## 4. Other reference database 

## Reef microbe? why bother
- Microbe responds to environmental changes/disturbances rapidly
- Can we use microbes as bioindicator of environmental stress
  - Coral and sediment too variable
  - Seawater consistent across time and space
  - Sponge too consistent
- Water MAG study presented functional(gene) more useful than taxonomic information as bioindicator
  - Ammoniamonoxygenase increased when Nitrogen concentration increased
- Database (Ocean Microbiome Database Pooli et al Nature 2022)
- Seawater MAGs (Robbins et al XXXX)
  -   GraftM showed Prochlorococcus, Synechococcus, Pelagibacterales abundant
  -   High strain heterogeneity makes assembly challenging
  -   SNP analysis in the future
  -   Patin and Goodwin mSystems 2022
  -   Illumina shows what is in every samples and Nanopore as a mean to get genomes
  -   Taxonomic bias between platforms (Illumina vs Nanopore)
  -   40% reads mapped to Prok MAGs, 10% to Viral MAGs, 15% to Euk MAGs
  -   

