# Combined gene- and genome-centric analyses

Analyses that bring together the [gene-centric](../1_gene_centric/README.md) and [genome-centric](../2_genome_centric/README.md) results: the effect of reef zoning status on community structure, using both mOTU-based and MAG-based profiles, and a direct comparison of what the two approaches reveal about taxonomic composition and metabolic potential.
 
| Section | Page | Scripts | Manuscript content |
|---|---|---|---|
| 6. Mixomics | [Community structure and reef zoning status](6_Mixomics/6_Mixomics.md) | [`6.01` to `6.10`](6_Mixomics/scripts/) | Fig. S3 (`6.05` onwards are supporting) |
| 7. Gene and genome analysis | [Comparison of gene- and genome-centric analyses](7_Gene_Genome_analysis/7_Gene_Genome_analysis.md) | [`7.01` and `7.02`](7_Gene_Genome_analysis/scripts/) | Fig. 3a-c, Fig. S2 |

## [6. Mixomics](6_Mixomics/6_Mixomics.md)
Run on four datasets: MAG-based and mOTU-based profiles of sediment and seawater.
- [Community profiles](6_Mixomics/6_Mixomics.md#profiles): CoverM and SingleM profiles, samples with complete design variables, and CLR values (`6.01`)
- [PERMANOVA](6_Mixomics/6_Mixomics.md#permanova): sequential PERMANOVA on Aitchison distances, zoning status + sampling trip + GBR sector + reef site (Fig. S3a; `6.02`)
- [MINT sPLS-DA](6_Mixomics/6_Mixomics.md#mint): prediction of zoning status with GBR sector as the study variable, and indicator features (Fig. S3b; `6.03`)
- [Permutation tests](6_Mixomics/6_Mixomics.md#permutation): zoning status permuted among reefs within sectors (Fig. S3c-d; `6.04`)
- [Supporting analyses](6_Mixomics/6_Mixomics.md#supporting): ordination and dbRDA, PLS-DA and sPLS-DA, MINT PLS-DA and choice of study variable, leave-one-reef-out validation, random forest, and SingleM RDA (`6.05` to `6.10`)

## [7. Gene and genome analysis](7_Gene_Genome_analysis/7_Gene_Genome_analysis.md)
- [Taxonomic composition](7_Gene_Genome_analysis/7_Gene_Genome_analysis.md#taxonomy): sediment community composition from SingleM and CoverM, as log2 fold changes, rank comparisons, and class-level composition per reef (Fig. S2; `7.01`)
- [Metabolic marker genes](7_Gene_Genome_analysis/7_Gene_Genome_analysis.md#metabolism): % of MAGs per phylum encoding each function, MAG enrichment, and gene-centric differential abundance in one heatmap (Fig. 3a-c; `7.02`)

## Inputs
- [SingleM](../1_gene_centric/1_SingleM/1_SingleM.md) profiles and microbial fractions
- [Protein](../1_gene_centric/2_Protein/2_Protein.md) differential abundance (`2.16`)
- [MAG recovery](../2_genome_centric/4_MAGs/4_MAG_recovery.md) read mapping (`4.12`)
- [MAG annotation](../2_genome_centric/5_MAGsAnnot/5_MAG_annotate.md) marker gene enrichment (`5.08`)
