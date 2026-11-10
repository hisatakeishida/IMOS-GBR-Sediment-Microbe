# Community-level functional profiling with a non-redundant protein catalogue

**Aim:** investigate the functional potential of sediment and seawater communities by building a non-redundant protein catalogue from all assemblies, annotating it (KOs, CAZymes, hydrogenases, metabolic marker genes, and *amoA*), and estimating protein abundance from reads.

Related manuscript content: Fig. 1e-f, Fig. 3c; Table S17.

Scripts for this analysis are in the top-level [`scripts/`](../../scripts/) folder (`2.01` to `2.16`). Assemblies are described in [Genome recovery](../../2_genome_centric/1_genome_recovery/1_genome_recovery.md).

| Step | Script | Description |
|---|---|---|
| 2.01 | [`scripts/2.01_rename_complete_orfs.sh`](../../scripts/2.01_rename_complete_orfs.sh) | Prefix Prodigal protein headers with habitat (and sample), concatenate per habitat, and keep complete ORFs (`partial=00`) with SeqKit v2.5.1 |
| 2.02 | [`scripts/2.02_mmseqs_cluster.sh`](../../scripts/2.02_mmseqs_cluster.sh) | Cluster complete ORFs from both habitats with MMseqs2 v15.6f452 `easy-linclust` (90% identity, 80% coverage) |
| 2.03 | [`scripts/2.03_dram_annotate_array.sh`](../../scripts/2.03_dram_annotate_array.sh) | KOfam annotation with DRAM v1.5 `annotate_genes`, in chunks of 100,000 proteins |
| 2.04 | [`scripts/2.04_dbcan_array.sh`](../../scripts/2.04_dbcan_array.sh) | CAZyme families with run_dbcan v5.0.1 `CAZyme_annotation`, one task per DRAM chunk |
| 2.05 | [`scripts/2.05_hyddb.sh`](../../scripts/2.05_hyddb.sh) | Hydrogenases with DIAMOND blastp v2.1.0 against HydDB, with class-specific identity thresholds |
| 2.06 | [`scripts/2.06_marker_genes_blastp_array.sh`](../../scripts/2.06_marker_genes_blastp_array.sh) | Metabolic marker genes with DIAMOND blastp v2.1.0, one task per marker gene database, with gene-specific identity thresholds |
| 2.07 | [`scripts/2.07_amoa_graftm.sh`](../../scripts/2.07_amoa_graftm.sh) | *amoA* with GraftM v0.15.1 (package 7.22) |
| 2.08 | [`scripts/2.08_combine_annotations.sh`](../../scripts/2.08_combine_annotations.sh) | Combine results into one table per annotation (KO, CAZyme, HydDB, marker genes) |
| 2.09 | [`scripts/2.09_diamond_makedb.sh`](../../scripts/2.09_diamond_makedb.sh) | DIAMOND database of the protein catalogue |
| 2.10 | [`scripts/2.10_protein_length.sh`](../../scripts/2.10_protein_length.sh) | Protein lengths in nucleotides (catalogue and SingleM marker genes), for RPKM, and DIAMOND database of the SingleM marker genes |
| 2.11 | [`scripts/2.11_diamond_blastx_array.sh`](../../scripts/2.11_diamond_blastx_array.sh) | Subsample R1 reads (≥ 140 bp) to 10 M and align them with DIAMOND blastx v2.1.0 to the catalogue (`DB=catalogue`) or the SingleM marker genes (`DB=singlem`) |
| 2.12 | [`scripts/2.12_diamond_counts.sh`](../../scripts/2.12_diamond_counts.sh) | Count reads per protein and per SingleM marker gene (long format) |
| 2.13 | [`scripts/2.13_protein_RPKM.R`](../../scripts/2.13_protein_RPKM.R) | RPKM of catalogue proteins normalised by mean SingleM marker gene RPKM (phyloseq object with KO annotation) |
| 2.14 | [`scripts/2.14_protein_annotation_objects.R`](../../scripts/2.14_protein_annotation_objects.R) | Phyloseq objects of normalised RPKM for KOs, marker genes with HydDB subgroups, and CAZyme families |
| 2.15 | [`scripts/2.15_protein_ordination.R`](../../scripts/2.15_protein_ordination.R) | PCoA on Aitchison distances of normalised RPKM for KOs (Fig. 1e) and CAZymes (Fig. 1f) |
| 2.16 | [`scripts/2.16_protein_differential_abundance.R`](../../scripts/2.16_protein_differential_abundance.R) | Differential abundance of metabolic marker genes (with HydDB subgroups) between sediment and seawater with MaAsLin2 and LinDA (Fig. 3c) |

## Table of contents
1. [Gene catalogue](#catalogue)
2. [Functional annotation](#annotation)
3. [Read mapping](#mapping)
4. [RPKM normalisation and ordination](#rpkm)
5. [Differential abundance](#diffabund)

## 1. Gene catalogue <a name="catalogue"></a>
Directory: [`1_gene_catalogue/`](1_gene_catalogue/)

- Prodigal v2.6.3 predicts open reading frames (ORFs) on the contigs of each assembly
- Protein headers are prefixed with habitat (and sample), proteins are concatenated per habitat, and only complete ORFs (`partial=00`) are kept
- MMseqs2 `easy-linclust` clusters the complete proteins from both habitats at 90% identity across 80% of protein length
- Platform: Bash (SLURM)
- Input: assembled contigs (`<sample>.final_contigs.fasta`)
- Output: non-redundant protein catalogue (`imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta`) and its list of protein IDs
- Scripts: [`2.01_rename_complete_orfs.sh`](../../scripts/2.01_rename_complete_orfs.sh), [`2.02_mmseqs_cluster.sh`](../../scripts/2.02_mmseqs_cluster.sh)

```bash
# Gene prediction, per assembly
prodigal -i ${sample_name}.final_contigs.fasta -p meta -f gff -q \
    -o ${sample_name}.prodigal.gff \
    -a ${sample_name}.protein.translations.faa \
    -d ${sample_name}.nucleotide.fna

sbatch scripts/2.01_rename_complete_orfs.sh
sbatch scripts/2.02_mmseqs_cluster.sh
```

| Data | Num seqs |
|---|---|
| Sediment, all predicted proteins | 1,853,838,210 |
| Sediment, complete ORFs | 22,803,724 |
| Seawater, all predicted proteins | 182,544,749 |
| Seawater, complete ORFs | 38,136,288 |
| Both habitats, complete ORFs | 60,940,012 |
| Non-redundant catalogue | 39,225,676 |

## 2. Functional annotation <a name="annotation"></a>
Directory: [`2_annotation/`](2_annotation/)

- KOfam: DRAM v1.5 `annotate_genes`, run on 393 chunks of 100,000 proteins
- CAZyme families: run_dbcan v5.0.1 `CAZyme_annotation` on the same chunks (recommended family from `overview.tsv`)
- Hydrogenases: DIAMOND blastp v2.1.0 against HydDB, keeping hits above class-specific identity ([NiFe] Group 4 > 50%, other [NiFe] > 30%, [FeFe] > 45%, [Fe] > 50%)
- Metabolic marker genes: DIAMOND blastp v2.1.0 against each marker gene database (query cover ≥ 80%, gene-specific identity of 30% to 80%), keeping the best hit per protein across databases
- *amoA*: GraftM v0.15.1 (package 7.22), separating bacterial and archaeal *amoA* (Table S17)
- Platform: Bash (SLURM)
- Input: non-redundant protein catalogue from section 1
- Output: `ko_annot.tsv`, `ko_def.tsv`, `dbcan_annot.tsv`, `hyddb_annot.tsv`, `marker_genes_annot.tsv`, and GraftM output
- Scripts: [`2.03`](../../scripts/2.03_dram_annotate_array.sh) to [`2.08`](../../scripts/2.08_combine_annotations.sh)

```bash
sbatch scripts/2.03_dram_annotate_array.sh
sbatch scripts/2.04_dbcan_array.sh                  # after 2.03
sbatch scripts/2.05_hyddb.sh
sbatch scripts/2.06_marker_genes_blastp_array.sh
sbatch scripts/2.07_amoa_graftm.sh
sbatch scripts/2.08_combine_annotations.sh          # after 2.03 to 2.06
```

## 3. Read mapping <a name="mapping"></a>
Directory: [`3_read_mapping_RPKM/`](3_read_mapping_RPKM/)

- R1 reads ≥ 140 bp are subsampled to 10 M reads per sample (seqtk, seed 100); samples with ≤ 10 M reads are excluded
- Reads are aligned to the whole protein catalogue with DIAMOND blastx v2.1.0 (E-value 1e-5, min-score 40, query cover 80%, identity 70%, max-hsps 1, max-target-seqs 1)
- Reads are also aligned to the single-copy marker genes of the SingleM GTDB R10-RS226 metapackage with the same settings
- Platform: Bash (SLURM, one task per sample)
- Input: protein catalogue from section 1; adapter-removed R1 reads
- Hits are counted per protein and per marker gene for each sample
- Output: DIAMOND databases, protein lengths (`comp_protein_length.fixed.tsv`, `singlem_marker_length.tsv`), per-sample DIAMOND hit tables, and long-format counts (`diamond_long.catalogue.tsv`, `diamond_long.singlem.tsv`)
- Scripts: [`2.09_diamond_makedb.sh`](../../scripts/2.09_diamond_makedb.sh), [`2.10_protein_length.sh`](../../scripts/2.10_protein_length.sh), [`2.11_diamond_blastx_array.sh`](../../scripts/2.11_diamond_blastx_array.sh), [`2.12_diamond_counts.sh`](../../scripts/2.12_diamond_counts.sh)

```bash
sbatch scripts/2.09_diamond_makedb.sh
sbatch scripts/2.10_protein_length.sh
# One array per habitat and database (catalogue or SingleM markers)
sbatch --export=ALL,HABITAT=sed,DB=catalogue --array=1-191 scripts/2.11_diamond_blastx_array.sh
sbatch --export=ALL,HABITAT=sea,DB=catalogue --array=1-182 scripts/2.11_diamond_blastx_array.sh
sbatch --export=ALL,HABITAT=sed,DB=singlem   --array=1-191 scripts/2.11_diamond_blastx_array.sh
sbatch --export=ALL,HABITAT=sea,DB=singlem   --array=1-182 scripts/2.11_diamond_blastx_array.sh

sbatch scripts/2.12_diamond_counts.sh              # after 2.11
```

## 4. RPKM normalisation and ordination <a name="rpkm"></a>
Directory: [`3_read_mapping_RPKM/`](3_read_mapping_RPKM/)

- Removes protein counts representing < 0.0001% of a sample's mapped reads
- Converts counts to reads per kilobase per million (RPKM) for proteins and SingleM marker genes, then divides protein RPKM by the mean RPKM of all SingleM marker genes in each sample
- Builds phyloseq objects of normalised RPKM (proteins as taxa) for KOs, metabolic marker genes with HydDB subgroups (HydDB replaces the generic [NiFe] and [FeFe] hydrogenase markers), and CAZyme families
- PCoA on Aitchison distances for KOs (each KO-annotated protein as a feature; Fig. 1e) and CAZymes (proteins summed by family; Fig. 1f)
- Platform: R (phyloseq, microViz, data.table)
- Input: long-format counts and lengths from section 3, annotation tables from section 2, sample metadata (xlsx)
- Output: `ko.COUNT.rds`, `ko.RPKM.rds`, `ps_ko.RDS`, `ps_greening_hyddb.RDS`, `ps_dbcan.RDS`, mapped reads per sample, and Fig. 1e-f
- Scripts: [`2.13_protein_RPKM.R`](../../scripts/2.13_protein_RPKM.R), [`2.14_protein_annotation_objects.R`](../../scripts/2.14_protein_annotation_objects.R), [`2.15_protein_ordination.R`](../../scripts/2.15_protein_ordination.R)

```bash
Rscript scripts/2.13_protein_RPKM.R
Rscript scripts/2.14_protein_annotation_objects.R
Rscript scripts/2.15_protein_ordination.R
```

## 5. Differential abundance <a name="diffabund"></a>
Directory: [`4_differential_abundance/`](4_differential_abundance/)

- Metabolic marker genes (with HydDB subgroups), as relative abundance within the gene set
- MaAsLin2 v1.16.0: TSS normalisation, AST transformation, habitat and season as fixed effects, reef as a random effect, features in ≥ 5% of samples
- LinDA v0.2.0: proportions, `~ habitat + season + (1 | reef)`, features in ≥ 5% of samples
- Sediment is the reference level; a gene is differentially abundant if both methods give adjusted p ≤ 0.05 for habitat with the same direction (Fig. 3c)
- Platform: R (phyloseq, microViz, Maaslin2, MicrobiomeStat)
- Input: `ps_greening_hyddb.RDS` from section 4
- Output: `MAASLIN_LINDA_greeningDB.xlsx`
- Script: [`2.16_protein_differential_abundance.R`](../../scripts/2.16_protein_differential_abundance.R)