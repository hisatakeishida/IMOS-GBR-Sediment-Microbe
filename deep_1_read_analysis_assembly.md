# Read-based community profiling, Metagenomic Assembly, and Binning

## Table of Contents
2. [Microbial Community Profiling](#div_classify)
3. [Metagenomic assembly and Binning of Bacterial and Archeal MAGs](#assembly_bining)
4. [Functional annotation](#func_annot)
5. 
6. 

## 2. Microbial Community Profiling using SingleM <a name="div_classify"></a>
- Search for GTDB OTU data in metagenomic reads using SingleM
- Singlem summarise
- Relative abundance 
- Beta diversity 

## 3. Metagenomic assembly and Binning of Bacterial and Archeal MAGs <a name="assembly_bining"></a>
- 27 metagenome individually assembled using hybrid (nanopore + illumina) with Aviary
- Aggregate refined MAGs and dereplicated using dRep
- SingleM to determine genome recovery at genus-level 
- Taxonomic annotation of each MAG using GTDB-Tk
- CheckM to examine completeness
- CoverM to calculate relative abundance 

## 3. Functional annotation <a name="func_annot"></a>
- MetaErg + CD-HIT
- METABOLIC
- EnrichM
