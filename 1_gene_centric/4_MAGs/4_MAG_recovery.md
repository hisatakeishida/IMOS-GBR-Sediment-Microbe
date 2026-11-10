# MAG recovery

**Aim:** recover prokaryotic MAGs from sediment metagenomes by single-sample assembly and co-assembly with multi-sample binning, then assess their quality, dereplicate them to species level, assign taxonomy, infer phylogenies, and assess their novelty against GTDB and public genome collections. Reads are mapped back to the MAGs and to unbinned contigs to assess how much of each community the MAGs represent.

Related manuscript content: Fig. 1b, Fig. 2a-d, Fig. S1; Table 1, Tables S2, S7, and S8.

Scripts for this analysis are in the top-level [`scripts/`](../../scripts/) folder (`4.01` to `4.16`).

| Step | Script | Description |
|---|---|---|
| 4.01 | [`scripts/4.01_aviary_assemble_array.sh`](../../scripts/4.01_aviary_assemble_array.sh) | Single-sample assembly with Aviary assemble (hybrid or Illumina-only) |
| 4.02 | [`scripts/4.02_singlem_reads_assemblies_array.sh`](../../scripts/4.02_singlem_reads_assemblies_array.sh) | SingleM pipe (GTDB R09-RS220) on each sample's reads and assembly |
| 4.03 | [`scripts/4.03_select_binning_samples.sh`](../../scripts/4.03_select_binning_samples.sh) | SingleM appraise and selection of 20 binning datasets per assembly |
| 4.04 | [`scripts/4.04_aviary_recover_array.sh`](../../scripts/4.04_aviary_recover_array.sh) | Multi-sample binning with Aviary recover |
| 4.05 | [`scripts/4.05_binchicken_coassemble.sh`](../../scripts/4.05_binchicken_coassemble.sh) | Choice of sample pairs for co-assembly, and of co-binning datasets, with BinChicken coassemble |
| 4.06 | [`scripts/4.06_aviary_coassembly_array.sh`](../../scripts/4.06_aviary_coassembly_array.sh) | Co-assembly and co-binning with Aviary (`MODE=assemble` or `MODE=recover`) |
| 4.07 | [`scripts/4.07_mag_quality.sh`](../../scripts/4.07_mag_quality.sh) | CheckM and CheckM2 quality; MAGs with quality score ≥ 50% by either tool |
| 4.08 | [`scripts/4.08_dereplicate.sh`](../../scripts/4.08_dereplicate.sh) | Dereplication at 95% ANI with CoverM cluster (`STAGE=1` or `STAGE=2`) |
| 4.09 | [`scripts/4.09_gtdbtk_taxonomy_phylogeny.sh`](../../scripts/4.09_gtdbtk_taxonomy_phylogeny.sh) | Taxonomy with GTDB-Tk classify_wf, de_novo_wf trees and PhyloRank for curation, and bacterial and archaeal phylogenies with IQ-TREE |
| 4.10 | [`scripts/4.10_novelty_cluster.sh`](../../scripts/4.10_novelty_cluster.sh) | Clustering with public genome collections at 95% ANI, and genomes per collection in each cluster |
| 4.11 | [`scripts/4.11_novelty_upset.R`](../../scripts/4.11_novelty_upset.R) | UpSet plot of clusters shared between collections (Fig. 2d) |
| 4.12 | [`scripts/4.12_coverm_mags_array.sh`](../../scripts/4.12_coverm_mags_array.sh) | Read mapping to MAGs with CoverM genome |
| 4.13 | [`scripts/4.13_unmapped_reads_array.sh`](../../scripts/4.13_unmapped_reads_array.sh) | Reads not mapped to sediment MAGs |
| 4.14 | [`scripts/4.14_unbinned_contigs_array.sh`](../../scripts/4.14_unbinned_contigs_array.sh) | Unbinned contigs: MMseqs2 taxonomy and reads per domain |
| 4.15 | [`scripts/4.15_genome_features.R`](../../scripts/4.15_genome_features.R) | Genome size, GC content, and predicted genes of sediment and seawater MAGs (Fig. 2e-g) |
| 4.16 | [`scripts/4.16_phylum_counts.R`](../../scripts/4.16_phylum_counts.R) | Number of sediment and seawater MAGs per phylum (Fig. 2b) |

## Table of contents
1. [Single-sample assembly](#assembly)
2. [Selection of datasets for binning](#selection)
3. [Multi-sample binning](#binning)
4. [Co-assembly and co-binning](#coassembly)
5. [MAG quality](#qc)
6. [Dereplication](#derep)
7. [Taxonomy and phylogeny](#taxonomy)
8. [Novelty](#novelty)
9. [Read mapping to MAGs](#mapping)
10. [Unbinned contigs](#unbinned)

## 1. Single-sample assembly <a name="assembly"></a>
Directory: [`1_assembly/`](1_assembly/)

- Aviary assemble v0.9.1: metaFlye v2.9.4 (ONT reads) with metaSPAdes v3.15.4 (Illumina reads) for the 17 samples with both data types; metaSPAdes v3.15.4 for the other 174 samples
- Platform: Bash (SLURM array, one task per sample)
- Input: adapter-removed Illumina reads (fastp v0.23.4) and ONT reads (porechop v0.2.4) (see [Samples, sequencing, and QC](../../0_metadata/0_sample_data_qc.md)), listed in a sample table (`<sample> <R1> <R2> [<ONT>]`)
- Output: 191 single-sample assemblies (17 hybrid, 174 Illumina-only; Table S2, Fig. 1b)
- Script: [`4.01_aviary_assemble_array.sh`](../../scripts/4.01_aviary_assemble_array.sh)

## 2. Selection of datasets for binning <a name="selection"></a>
Directory: [`2_binning_selection/`](2_binning_selection/)

- SingleM v1.18.0 pipe (GTDB R09-RS220 metapackage; 59 single-copy marker genes) on each sample's reads and assembly
- SingleM appraise (`--imperfect`) of each assembly against the reads of all samples; the 20 samples with the highest percentage of their bacterial marker gene OTUs recovered in the assembly are used for its binning
- Platform: Bash (SLURM)
- Output: list of 20 binning datasets per assembly (Table S2)
- Scripts: [`4.02_singlem_reads_assemblies_array.sh`](../../scripts/4.02_singlem_reads_assemblies_array.sh), [`4.03_select_binning_samples.sh`](../../scripts/4.03_select_binning_samples.sh)

## 3. Multi-sample binning <a name="binning"></a>
Directory: [`3_binning/`](3_binning/)

- Aviary recover v0.9.1 (`--binning-only`), with DAS Tool v1.1.2 selecting the best bins from MetaBAT1 v2.15, MetaBAT2 v2.15, SemiBin v2.1.0, VAMB v3.0.2, and Rosella v0.5.2; the sample's ONT reads are also used for hybrid assemblies
- Aviary also runs CheckM and CheckM2 on the final bins (used in section 5)
- Platform: Bash (SLURM array, one task per assembly)
- Script: [`4.04_aviary_recover_array.sh`](../../scripts/4.04_aviary_recover_array.sh)

## 4. Co-assembly and co-binning <a name="coassembly"></a>
Directory: [`4_coassembly/`](4_coassembly/)

- BinChicken v0.12.11 coassemble chooses pairs of samples with similar composition (59 marker genes, as above) and 20 datasets per co-assembly for co-binning, and writes the Aviary commands
- Co-assembly with MEGAHIT v1.2.9 within Aviary assemble; co-binning with Aviary recover
- Platform: Bash (SLURM)
- Output: 71 co-assemblies and their bins (Table S2, Fig. 1b)
- Scripts: [`4.05_binchicken_coassemble.sh`](../../scripts/4.05_binchicken_coassemble.sh), [`4.06_aviary_coassembly_array.sh`](../../scripts/4.06_aviary_coassembly_array.sh)

```bash
sbatch scripts/4.05_binchicken_coassemble.sh
sbatch --export=ALL,MODE=assemble scripts/4.06_aviary_coassembly_array.sh
sbatch --export=ALL,MODE=recover  scripts/4.06_aviary_coassembly_array.sh   # after assembly
```

## 5. MAG quality <a name="qc"></a>
Directory: [`5_MAG_QC/`](5_MAG_QC/)

- CheckM v1.2.2 and CheckM2 v1.0.2 results from Aviary recover are combined; bins are kept if either tool gives a quality score (completeness − 3 × contamination) ≥ 50%
- Platform: Bash
- Input: binning outputs from sections 3 and 4
- Output: 2,698 MAGs, their quality table (Table S7), and a CheckM2 table for dereplication
- Script: [`4.07_mag_quality.sh`](../../scripts/4.07_mag_quality.sh)

## 6. Dereplication <a name="derep"></a>
Directory: [`6_dereplication/`](6_dereplication/)

- CoverM cluster v0.7.0 at 95% ANI (minimum aligned fraction 15%), with representatives chosen by the `Parks2020_reduced` formula (CheckM2 values)
- Two stages: single-sample assembly MAGs first, then their representatives together with co-assembly MAGs; eukaryotic MAGs are then removed
- Platform: Bash (SLURM)
- Output: 1,201 species-representative sediment MAGs (mean completeness 74.17 ± 13.35%, contamination 2.42 ± 2.00%, CheckM2; Table S8)
- Script: [`4.08_dereplicate.sh`](../../scripts/4.08_dereplicate.sh)

```bash
sbatch --export=ALL,STAGE=1 scripts/4.08_dereplicate.sh
sbatch --export=ALL,STAGE=2 scripts/4.08_dereplicate.sh   # after stage 1
```

## 7. Taxonomy and phylogeny <a name="taxonomy"></a>
Directory: [`7_taxonomy/`](7_taxonomy/)

- GTDB-Tk v2.4.1 `classify_wf` with GTDB R10-RS226, for the 1,201 sediment MAGs and 876 seawater MAGs (MGD-seawater; Robbins et al., 2026)
- Taxonomy curated manually by inspecting GTDB-Tk `de_novo_wf` genome trees (outgroups p__Chloroflexota for bacteria and p__Undinarchaeota for archaea) decorated with PhyloRank v0.1.12 (`decorate`, `outliers`), comparing relative evolutionary distance (RED) at each node with median RED per rank
- Phylogenies from concatenated GTDB-Tk bac120 (bacteria) and arc53 (archaea) marker gene alignments with IQ-TREE v2.2.2.3 (1,000 ultrafast bootstraps): LG + G (bacteria), LG + C10 + F + G (archaea)
- Platform: Bash (SLURM)
- Genome size, GC content, and number of predicted genes (CheckM2) compared between sediment and seawater MAGs with Wilcoxon rank-sum tests
- Output: curated taxonomy of 2,077 MAGs (Fig. S1); bacterial and archaeal trees (Fig. 2a-b); MAGs per phylum (Fig. 2b); genome features (Fig. 2e-g)
- Script: [`4.09_gtdbtk_taxonomy_phylogeny.sh`](../../scripts/4.09_gtdbtk_taxonomy_phylogeny.sh) (curation itself is manual), [`4.15_genome_features.R`](../../scripts/4.15_genome_features.R), [`4.16_phylum_counts.R`](../../scripts/4.16_phylum_counts.R)

## 8. Novelty <a name="novelty"></a>
Directory: [`9_novelty/`](9_novelty/)

- CoverM cluster v0.7.0 at 95% ANI (minimum aligned fraction 15%; `Parks2020_reduced`) of sediment and seawater MAGs together with public MAG collections: Ocean Microbiomics Database (OMD), Reef Microbiomics Database (RMD), Genomes from Earth's Microbiomes (GEMs), shallow-water coastal marine sediments (Silva-Solar et al., 2024; Viver et al., 2025), and South China Sea reef sediments (Dong et al., 2023)
- Number of genomes from each collection in each cluster, and UpSet plot of clusters shared between collections
- Representation in GTDB R10-RS226 from the GTDB-Tk species assignments (section 7)
- Platform: Bash (SLURM) and R (ComplexUpset)
- Output: 1,182 of 1,201 sediment species without a GTDB representative, and 1,158 not clustering with any public genome (Table 1, Fig. 2c-d)
- Scripts: [`4.10_novelty_cluster.sh`](../../scripts/4.10_novelty_cluster.sh), [`4.11_novelty_upset.R`](../../scripts/4.11_novelty_upset.R)

## 9. Read mapping to MAGs <a name="mapping"></a>
Directory: [`10_read_mapping/`](10_read_mapping/)

- CoverM genome, mapping each habitat's reads to its own MAGs (sediment reads to 1,201 sediment MAGs; seawater reads to 876 seawater MAGs), excluding alignments with < 95% identity or < 75% aligned length
- Platform: Bash (SLURM array, one task per sample)
- Output: MAG-based community profiles (used in [Mixomics](../6_Mixomics/6_Mixomics.md)); 0.48% to 18.35% of sediment reads mapped per sample (mean 4.43%)
- Script: [`4.12_coverm_mags_array.sh`](../../scripts/4.12_coverm_mags_array.sh)

```bash
sbatch --export=ALL,HABITAT=sed,MODE=count --array=1-191 scripts/4.12_coverm_mags_array.sh
```

## 10. Unbinned contigs <a name="unbinned"></a>
Directory: [`11_unbinned_contigs/`](11_unbinned_contigs/)

- Reads not mapped to sediment MAGs (≥ 95% identity over ≥ 75% of the read) are extracted for each sample
- Unbinned contigs of each single-sample assembly are classified with MMseqs2 v15.6f452 `easy-taxonomy` against a custom database of RefSeq, EukProt, PLSDB, IMG/VR, and proteins predicted from GTDB R08-RS214 genomes (Robbins et al., 2026)
- MAG-unmapped reads from all samples of the same reef are mapped to the unbinned contigs (CoverM contig), and mapped reads are summed per domain
- Platform: Bash (SLURM array)
- Output: reads per domain for unbinned contigs (mean 91.54% Bacteria, 0.53% Archaea, 1.06% viruses, 0.26% eukaryotes)
- Scripts: [`4.13_unmapped_reads_array.sh`](../../scripts/4.13_unmapped_reads_array.sh), [`4.14_unbinned_contigs_array.sh`](../../scripts/4.14_unbinned_contigs_array.sh)