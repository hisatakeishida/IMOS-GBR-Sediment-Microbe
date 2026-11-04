# Metagenomics analysis of sediment microbes across the Great Barrier Reef

This repository contains code and data associated with the paper:

> Ishida, H., Robbins, S. J., Bell, S. C., Lim, J. Y., Webster, N., Bourne, D. G., Laffy, P. W., Hugenholtz, P., Yeoh, Y. K., & Chan, C. X. **Comparative metagenomics reveals diversity and metabolic heterogeneity of sediment microbes from the Great Barrier Reef.** *In preparation.*

## Overview
This project implemented gene-centric and genome-centric analyses to investigate the taxonomic and metabolic diversity of surface sediment microbiomes of the Great Barrier Reef (GBR), relative to microbes in the co-located surface seawater ([Robbins et al., 2026](https://doi.org/10.1038/s41586-026-10778-z)). Sediments were collected from 48 reefs (23 no-take marine reserves and 25 fished reefs) across seven GBR sectors over four sampling trips (November 2019 to July 2020), yielding 191 Illumina metagenomes (17 with ONT long reads), 2,698 metagenome-assembled genomes (MAGs), and 1,201 species-level MAGs (MGD-sediment).

This project is part of the Great Barrier Reef Microbial Genomes Database sub-Facility within the Marine Microbiome Initiative Facility of Australia's Integrated Marine Observing System (IMOS). More details are available on the [IMOS GBR Microbial Genomic Database page](https://imos.org.au/facilities/marinemicrobiomeinitiative/gbr-microbial-genomic-database).

More details for the seawater data analysed alongside sediments (MGD-seawater) are available on Dr Marko Terzin's github pages and Zendo https://zenodo.org/records/20503966.
- [Gene-centric analyses](https://github.com/mterzin/IMOS_GBR_MGD_read-centric_analysis)  associated with  [Terzin et al., 2025](https://doi.org/10.1038/s41586-026-10778-z)
- [Genome-centric analyses](https://github.com/mterzin/fishy_microbes) associated with [Terzin et al., 2026](https://doi.org/10.1038/s41586-026-10778-z) and [Robbins et al., 2026](https://doi.org/10.1038/s41586-026-10778-z)

## Data availability
Raw sequencing data will be available upon publication on the IMOS-AODN repository (XXXXXX) and ENA under project number (XXXXXX). Sediment MAGs generated from this study will be archived on Zenodo (XXXXXX) upon publication. 

## Repository structure
### [`0_metadata/`](0_metadata/): samples, sequencing, and QC ([page](0_metadata/0_sample_data_qc.md))
- 192 sediment samples, 48 reefs, zoning status, and sampling trips
- Illumina and ONT library preparation and sequencing yield
- Adapter removal with fastp and porechop

### [`4_seawater/`](4_seawater/): ([page](4_seawater/4_seawater.md))


### [`1_gene_centric/`](1_gene_centric/): read-based analyses

#### [`1_SingleM/`](1_gene_centric/1_SingleM/): mOTU-based analyses ([page](1_gene_centric/1_SingleM/1_SingleM.md))
- [`1_community_profiles/`](1_gene_centric/1_SingleM/1_community_profiles/): SingleM pipe on reads subsampled to 10 M read pairs (Fig. 1c)
- [`2_ordination_PERMANOVA/`](1_gene_centric/1_SingleM/2_ordination_PERMANOVA/): PCA and sequential PERMANOVA (Fig. 1d, 1g-h)
- [`3_differential_abundance/`](1_gene_centric/1_SingleM/3_differential_abundance/): MaAsLin2 and LinDA, sediment vs seawater (Table S4)
- [`4_abundance_prevalence/`](1_gene_centric/1_SingleM/4_abundance_prevalence/): family-level abundance and prevalence (Tables S5 and S6)
- [`5_appraise_novelty/`](1_gene_centric/1_SingleM/5_appraise_novelty/): SingleM appraise against assemblies, MAGs, and GTDB R10-RS226; microbial_fraction

#### [`2_Protein/`](1_gene_centric/2_Protein/): non-redundant protein analyses ([page](1_gene_centric/2_Protein/2_Protein.md))
- [`1_gene_catalogue/`](1_gene_centric/2_Protein/1_gene_catalogue/): Prodigal ORFs clustered with MMseqs2 easy-linclust (39.2 M non-redundant proteins)
- [`2_annotation/`](1_gene_centric/2_Protein/2_annotation/): KOfam (DRAM), CAZymes (run_dbcan), hydrogenases (HydDB), 44 marker genes, and *amoA* (GraftM; Table S17)
- [`3_read_mapping_RPKM/`](1_gene_centric/2_Protein/3_read_mapping_RPKM/): DIAMOND blastx, RPKM normalised to SingleM marker genes, and PCoA (Fig. 1e-f)
- [`4_differential_abundance/`](1_gene_centric/2_Protein/4_differential_abundance/): MaAsLin2 and LinDA on functional profiles (Fig. 3c)

#### [`3_Sylph/`](1_gene_centric/3_Sylph/): k-mer-based analyses ([page](1_gene_centric/3_Sylph/3_Sylph.md))
- Supporting analysis; not reported in the manuscript

### [`2_genome_centric/`](2_genome_centric/): MAG-based analyses

#### [`1_genome_recovery/`](2_genome_centric/1_genome_recovery/): MAG recovery, QC, dereplication, taxonomy, and novelty ([page](2_genome_centric/1_genome_recovery/1_genome_recovery.md))
- [`1_assembly/`](2_genome_centric/1_genome_recovery/1_assembly/): Aviary assemble, 17 hybrid and 174 Illumina-only assemblies (Table S2)
- [`2_coassembly/`](2_genome_centric/1_genome_recovery/2_coassembly/): BinChicken coassemble, 71 co-assemblies (Table S2)
- [`3_binning/`](2_genome_centric/1_genome_recovery/3_binning/): Aviary recover multi-sample binning
- [`4_MAG_QC/`](2_genome_centric/1_genome_recovery/4_MAG_QC/): CheckM and CheckM2, quality score ≥ 50% (2,698 MAGs; Table S7)
- [`5_dereplication/`](2_genome_centric/1_genome_recovery/5_dereplication/): CoverM cluster at 95% ANI (1,201 MAGs; Table S8)
- [`6_taxonomy/`](2_genome_centric/1_genome_recovery/6_taxonomy/): GTDB-Tk R10-RS226 and PhyloRank curation (Fig. S1)
- [`7_phylogeny/`](2_genome_centric/1_genome_recovery/7_phylogeny/): IQ-TREE bacterial and archaeal trees, and genome features (Fig. 2a-b, 2e-g)
- [`8_novelty/`](2_genome_centric/1_genome_recovery/8_novelty/): comparison with GTDB and public genome collections (Table 1, Fig. 2c-d)

#### [`2_genome_annot/`](2_genome_centric/2_genome_annot/): MAG annotation ([page](2_genome_centric/2_genome_annot/2_genome_annot.md))
- [`1_DRAM/`](2_genome_centric/2_genome_annot/1_DRAM/): KOfam, Pfam, and MEROPS annotation
- [`2_CAZymes_CGCs/`](2_genome_centric/2_genome_annot/2_CAZymes_CGCs/): run_dbcan CAZymes, CGCs, and dbCAN-PUL substrates (Fig. 4, Table S20)
- [`3_sulphatases/`](2_genome_centric/2_genome_annot/3_sulphatases/): SulfAtlas
- [`4_hydrogenases_marker_genes/`](2_genome_centric/2_genome_annot/4_hydrogenases_marker_genes/): HydDB and 44 metabolic marker genes (Fig. 3a-b, 3d; Table S16)
- [`5_amoA/`](2_genome_centric/2_genome_annot/5_amoA/): GraftM *amoA* detection and clade assignment (Table S18)
- [`6_enrichment_tests/`](2_genome_centric/2_genome_annot/6_enrichment_tests/): Mann-Whitney U tests on 796 MAGs (Tables S10 to S15, S19)
- [`7_BGCs/`](2_genome_centric/2_genome_annot/7_BGCs/): antiSMASH, clust-o-matic GCFs, and MIBiG and BiG-FAM novelty (Fig. 5, Tables S21 and S22)

#### [`3_genome_abundance/`](2_genome_centric/3_genome_abundance/): MAG abundance ([page](2_genome_centric/3_genome_abundance/3_genome_abundance.md))
- [`1_CoverM_mapping/`](2_genome_centric/3_genome_abundance/1_CoverM_mapping/): read mapping of each habitat to its own MAGs (MAG-based community profiles)
- [`2_unbinned_contigs/`](2_genome_centric/3_genome_abundance/2_unbinned_contigs/): MMseqs2 taxonomy of unbinned contigs

#### [`4_Mixomics/`](2_genome_centric/4_Mixomics/): reef zoning status analyses ([page](2_genome_centric/4_Mixomics/4_Mixomics.md))
- [`1_PERMANOVA/`](2_genome_centric/4_Mixomics/1_PERMANOVA/): sequential PERMANOVA on MAG-based profiles (Fig. S3a)
- [`2_MINT_sPLSDA/`](2_genome_centric/4_Mixomics/2_MINT_sPLSDA/): tuning, final models, and indicator MAGs (Fig. S3b)
- [`3_permutation_tests/`](2_genome_centric/4_Mixomics/3_permutation_tests/): reef-level permutation nulls for MAG and mOTU profiles (Fig. S3c-d)

### [`3_SingleM_VS_CoverM/`](3_SingleM_VS_CoverM/): gene- vs genome-centric abundance ([page](3_SingleM_VS_CoverM/3_SingleM_VS_CoverM.md))
- Relative abundance compared across taxonomic ranks and sampling trips (Fig. S2, Tables S5 and S9)

### [`5_Rscripts/`](5_Rscripts/): statistics and visualisation
- [`1_MASH.R`](5_Rscripts/1_MASH.R): MASH visualisation
- [`2_SingleM.R`](5_Rscripts/2_SingleM.R): SingleM visualisation
- [`3_RPKM_calculation_KO.R`](5_Rscripts/3_RPKM_calculation_KO.R): RPKM calculation and visualisation (KO)
- [`4_RPKM_calculation_CAZyme.R`](5_Rscripts/4_RPKM_calculation_CAZyme.R): RPKM calculation and visualisation (CAZyme)

## Data availability

| Data | Location |
|---|---|
| Raw sequencing reads (Illumina and ONT) | IMOS-AODN (`TODO`) and ENA (`TODO`) |
| Sediment MAGs (MGD-sediment) | Zenodo (`TODO`) |
| Seawater MAGs (MGD-seawater) | <https://zenodo.org/records/20503966> |

## Authors and affiliations

Hisatake Ishida<sup>1</sup>, Steven J. Robbins<sup>1</sup>, Sara C. Bell<sup>2</sup>, Jin Yan Lim<sup>2</sup>, Nicole Webster<sup>1,3</sup>, David G. Bourne<sup>2,4,5</sup>, Patrick W. Laffy<sup>2,5</sup>, Philip Hugenholtz<sup>1</sup>, Yun Kit Yeoh<sup>2,5</sup>, Cheong Xin Chan<sup>1</sup>

1. Australian Centre for Ecogenomics, School of Chemistry and Molecular Biosciences, The University of Queensland, Brisbane, QLD 4072, Australia
2. Australian Institute of Marine Science, Townsville, QLD 4810, Australia
3. Institute for Marine and Antarctic Studies, University of Tasmania, Hobart, TAS 7001, Australia
4. College of Science and Engineering, James Cook University, Townsville, QLD 4811, Australia
5. AIMS@JCU, James Cook University, Townsville, QLD 4811, Australia

## Questions or issues

Please contact Hisatake Ishida at hisatake.ishida@uq.net.au, or open an [issue](../../issues).