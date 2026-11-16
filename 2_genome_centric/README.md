# Genome-centric analyses

Analyses based on metagenome-assembled genomes (MAGs): recovery, quality control, dereplication, taxonomy, phylogeny, and novelty of sediment MAGs, and annotation of the metabolic and biosynthetic potential of sediment and seawater MAGs.
 
| Section | Page | Scripts | Manuscript content |
|---|---|---|---|
| 4. MAG recovery | [MAG recovery](4_MAGs/4_MAG_recovery.md) | [`4.01` to `4.16`](4_MAGs/scripts/) | Fig. 1b, Fig. 2a-g, Fig. S1; Table 1, Tables S2, S7, and S8 |
| 5. MAG annotation | [MAG annotation](5_MAGsAnnot/5_MAG_annotate.md) | [`5.01` to `5.14`](5_MAGsAnnot/scripts/) | Figs 3 to 5; Tables S10 to S22 |

## [4. MAG recovery](4_MAGs/4_MAG_recovery.md)
191 sediment metagenomes (17 with ONT long reads) yielding 2,698 MAGs with quality score ≥ 50% and 1,201 species-representative MAGs (MGD-sediment).
- [Single-sample assembly](4_MAGs/4_MAG_recovery.md#assembly): Aviary assemble, 17 hybrid and 174 Illumina-only assemblies (Table S2; `4.01`)
- [Selection of datasets for binning](4_MAGs/4_MAG_recovery.md#selection): SingleM-based choice of 20 binning datasets per assembly (`4.02`, `4.03`)
- [Multi-sample binning](4_MAGs/4_MAG_recovery.md#binning): Aviary recover (`4.04`)
- [Co-assembly and co-binning](4_MAGs/4_MAG_recovery.md#coassembly): BinChicken, 71 co-assemblies (Table S2; `4.05`, `4.06`)
- [MAG quality](4_MAGs/4_MAG_recovery.md#qc): CheckM and CheckM2, quality score ≥ 50% by either tool (2,698 MAGs; Table S7; `4.07`)
- [Dereplication](4_MAGs/4_MAG_recovery.md#derep): CoverM cluster at 95% ANI (1,201 MAGs; Table S8; `4.08`)
- [Taxonomy and phylogeny](4_MAGs/4_MAG_recovery.md#taxonomy): GTDB-Tk R10-RS226, PhyloRank curation, IQ-TREE trees, and genome features (Fig. 2a-b, 2e-g, Fig. S1; `4.09`, `4.15`, `4.16`)
- [Novelty](4_MAGs/4_MAG_recovery.md#novelty): clustering with GTDB and public genome collections (Table 1, Fig. 2c-d; `4.10`, `4.11`)
- [Read mapping to MAGs](4_MAGs/4_MAG_recovery.md#mapping): CoverM read mapping, giving MAG-based community profiles (`4.12`)
- [Unbinned contigs](4_MAGs/4_MAG_recovery.md#unbinned): MMseqs2 taxonomy of contigs not placed in MAGs, and recruitment of reads not mapped to MAGs (`4.13`, `4.14`)

## [5. MAG annotation](5_MAGsAnnot/5_MAG_annotate.md)
Sediment (1,201) and seawater (876) MAGs; habitat comparisons use the 796 MAGs with CheckM2 completeness ≥ 85% (347 sediment, 449 seawater).
- [DRAM](5_MAGsAnnot/5_MAG_annotate.md#dram): KOfam, Pfam, and MEROPS annotation (`5.01`)
- [CAZymes and CGCs](5_MAGsAnnot/5_MAG_annotate.md#cazymes): run_dbcan CAZymes, CAZyme gene clusters, and substrates (Fig. 4, Table S20; `5.02`)
- [Sulphatases](5_MAGsAnnot/5_MAG_annotate.md#sulphatases): SulfAtlas (`5.03`)
- [Hydrogenases and metabolic marker genes](5_MAGsAnnot/5_MAG_annotate.md#markers): HydDB and metabolic marker genes (Fig. 3a-b and 3d, Table S16; `5.04`)
- [*amoA*](5_MAGsAnnot/5_MAG_annotate.md#amoa): GraftM detection and clade assignment (Table S18; `5.05`)
- [Enrichment tests](5_MAGsAnnot/5_MAG_annotate.md#enrichment): Mann-Whitney U tests between habitats with BH correction (Fig. 3b, Tables S10 to S15 and S19; `5.06` to `5.08`)
- [Biosynthetic gene clusters](5_MAGsAnnot/5_MAG_annotate.md#bgcs): antiSMASH, clust-o-matic gene cluster families with OMD, RMD, and MIBiG, and BiG-SLiCE distances (Fig. 5, Tables S21 and S22; `5.09` to `5.14`, with the helpers `clustomatic.modified.py` and `perform_l2norm_clustering.py`)

## Outputs used elsewhere
- MAG-based community profiles (`4.12`): [reef zoning analyses](../3_combine/6_Mixomics/6_Mixomics.md) and the [taxonomic comparison](../3_combine/7_Gene_Genome_analysis/7_Gene_Genome_analysis.md#taxonomy)
- Marker gene enrichment (`5.08`): [Fig. 3 heatmap](../3_combine/7_Gene_Genome_analysis/7_Gene_Genome_analysis.md#metabolism)
