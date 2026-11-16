# Gene-centric analyses

Read-based analyses of sediment and seawater metagenomes, independent of genome binning: taxonomic profiles from phylogenetic marker genes (mOTUs), functional profiles from a non-redundant protein catalogue, and a k-mer-based comparison against GTDB genomes.
 
| Section | Page | Scripts | Manuscript content |
|---|---|---|---|
| 1. SingleM | [mOTU-based community profiling](1_SingleM/1_SingleM.md) | [`1.01` to `1.06`](1_SingleM/scripts/) | Fig. 1c, 1d, 1g-h; Tables S4 to S6 |
| 2. Protein | [Functional profiling with a non-redundant protein catalogue](2_Protein/2_Protein.md) | [`2.01` to `2.16`](2_Protein/scripts/) | Fig. 1e-f, Fig. 3c; Table S17 |
| 3. Sylph | [k-mer-based taxonomic profiling](3_Sylph/3_Sylph.md) | [`3.01` and `3.02`](3_Sylph/scripts/) | Supporting; not reported |

## [1. SingleM](1_SingleM/1_SingleM.md)
SingleM v1.18 on reads subsampled to 10 million read pairs per sample.
- [Community profiles](1_SingleM/1_SingleM.md#profiles): read subsampling, SingleM pipe against GTDB R10-RS226, and microbial fraction (`1.01`, `1.02`)
- [Profile parsing](1_SingleM/1_SingleM.md#parse): coverage matrices per habitat (`1.03`)
- [Ordination and PERMANOVA](1_SingleM/1_SingleM.md#ordination): PCA, Shannon diversity, and sequential PERMANOVA (Fig. 1d, 1g-h; `1.04`)
- [Abundance, prevalence, and differential abundance](1_SingleM/1_SingleM.md#abundance): MaAsLin2 and LinDA, sediment versus seawater (Fig. 1c; Tables S4 to S6; `1.05`)
- [Marker gene novelty](1_SingleM/1_SingleM.md#appraise): SingleM appraise against assemblies, MAGs, and GTDB (`1.06`)

## [2. Protein](2_Protein/2_Protein.md)
Proteins predicted from all sediment and seawater assemblies, clustered into a non-redundant catalogue.
- [Gene catalogue](2_Protein/2_Protein.md#catalogue): complete ORFs clustered with MMseqs2 easy-linclust (0.9 identity, 0.8 coverage) (`2.01`, `2.02`)
- [Functional annotation](2_Protein/2_Protein.md#annotation): KOfam (DRAM), CAZymes (run_dbcan), hydrogenases (HydDB), metabolic marker genes, and *amoA* (GraftM; Table S17) (`2.03` to `2.08`)
- [Read mapping](2_Protein/2_Protein.md#mapping): DIAMOND blastx of R1 reads (≥ 140 bp) to the catalogue and to SingleM marker genes (`2.09` to `2.12`)
- [RPKM normalisation and ordination](2_Protein/2_Protein.md#rpkm): RPKM normalised by mean SingleM marker gene RPKM, and PCoA (Fig. 1e-f; `2.13` to `2.15`)
- [Differential abundance](2_Protein/2_Protein.md#diffabund): MaAsLin2 and LinDA on metabolic marker genes and hydrogenase subgroups (Fig. 3c; `2.16`)

## [3. Sylph](3_Sylph/3_Sylph.md)
- [Sylph profiles](3_Sylph/3_Sylph.md#sylph-profiles): sylph profile against GTDB R10-RS226 genomes, and merged tables (`3.01`, `3.02`)
- [Summary](3_Sylph/3_Sylph.md#summary): fewer than 1% of sediment reads, compared with about half of seawater reads, were assigned to GTDB genomes

## Outputs used elsewhere
- SingleM profiles and microbial fractions: mOTU-based profiles in the [reef zoning analyses](../3_combine/6_Mixomics/6_Mixomics.md) and the [gene- versus genome-centric comparison](../3_combine/7_Gene_Genome_analysis/7_Gene_Genome_analysis.md)
- Marker gene differential abundance (`2.16`): [Fig. 3 heatmap](../3_combine/7_Gene_Genome_analysis/7_Gene_Genome_analysis.md#metabolism)
