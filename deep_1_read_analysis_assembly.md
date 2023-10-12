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
```
conda activate singlem
singlem pipe -1 {inr1} -2 {inr2} --threads NUM -p ${base}_prof.tsv --taxonomic-profile-krona ${base}_prof_krona.html --otu-table ${base}_otu.tsv
```

## 3. Metagenomic assembly and Binning of Bacterial and Archeal MAGs <a name="assembly_bining"></a>
- 3.1 27 metagenome individually assembled using hybrid (nanopore + illumina) with Aviary
```

```

aviary recover --refinery-max-iterations 10? 

- 3.2 Aggregate refined MAGs and dereplicated using dRep
- 3.3 SingleM to determine genome recovery at genus-level 
- 3.4 Taxonomic annotation of each MAG using GTDB-Tk
- 3.5 CheckM to examine completeness
- 3.6 CoverM to calculate relative abundance 

## 3. Functional annotation <a name="func_annot"></a>
- MetaErg + CD-HIT
- METABOLIC
- EnrichM
