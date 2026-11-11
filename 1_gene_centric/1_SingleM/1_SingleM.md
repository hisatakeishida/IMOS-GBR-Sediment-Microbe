# mOTU-based community profiling with SingleM

**Aim:** characterise prokaryotic diversity in sediment and seawater communities using operational taxonomic units derived from phylogenetic marker genes (mOTUs), with [SingleM](https://github.com/wwood/singlem) v1.18.0 and the GTDB R10-RS226 metapackage.

Related manuscript content: Fig. 1c, 1d, 1g-h; Tables S4 to S6.

Scripts for this analysis are in the [`scripts/`](scripts/) folder (`1.01` to `1.06`).
| Step | Script | Description |
|---|---|---|
| 1.01 | [`scripts/1.01_subsample_10M_array.sh`](scripts/1.01_subsample_10M_array.sh) | Count read pairs (SeqKit v2.5.1) and subsample R1 and R2 to 10 M read pairs with seqtk (seed 100, so pairs stay matched). Samples with ≤ 10 M read pairs are excluded |
| 1.02 | [`scripts/1.02_singlem_pipe_array.sh`](scripts/1.02_singlem_pipe_array.sh) | `singlem pipe` and `singlem microbial_fraction` on the subsampled reads, with the GTDB R10-RS226 metapackage (`S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb`), run separately for sediment and seawater |
| 1.03 | [`scripts/1.03_singlem_prof_coverage.py`](scripts/1.03_singlem_prof_coverage.py) | Combine per-sample SingleM profiles into a single coverage matrix (rows: samples; columns: taxa; absent taxa = 0), one run per habitat |
| 1.04 | [`scripts/1.04_singlem_community_analysis.R`](scripts/1.04_singlem_community_analysis.R) | Build phyloseq objects scaled by microbial fraction; PCA, Shannon diversity, and PERMANOVA (habitat, then zoning status, sampling trip, GBR sector, and reef site within each habitat); reef-level PERMANOVA; sequencing facility QC |
| 1.05 | [`scripts/1.05_singleM_differential_abundance.R`](scripts/1.05_singleM_differential_abundance.R) | Family-level abundance and prevalence per habitat (Tables S5 and S6); top 15 class composition (Fig. 1c); differential abundance between habitats with MaAsLin2 and LinDA, consensus of both (Table S4) |
| 1.06 | [`scripts/1.06_singlem_appraise.sh`](scripts/1.06_singlem_appraise.sh) | `singlem pipe` on assemblies and MAGs, then `singlem appraise` of read mOTUs against assemblies, MAGs, and GTDB R10-RS226 genomes at species and genus level (89% identity), run separately for sediment and seawater |

## Table of contents
1. [Community profiles](#profiles)
2. [Profile parsing](#parse)
3. [Ordination and PERMANOVA](#ordination)
4. [Abundance, prevalence, and differential abundance](#abundance)
5. [Marker gene novelty (appraise)](#appraise)

## 1. Community profiles <a name="profiles"></a>
Directory: [`1_community_profiles/`](1_community_profiles/)

- Reads were subsampled to 10 million read pairs per sample to account for uneven sequencing depth
- `singlem pipe` generates mOTU tables and taxonomic profiles
- `singlem microbial_fraction` estimates the prokaryotic fraction of each metagenome
- Platform: Bash (SLURM arrays, one task per sample)
- Input: adapter-removed paired-end reads (see [Samples, sequencing, and QC](../../0_metadata/0_sample_data_qc.md))
- Output: subsampled reads, then per-sample taxonomic profiles (`prof/`), OTU tables (`otu/`), and microbial fractions (`mfrac/`) for each habitat

```bash
# 1.01: subsample (set LIST and READSOUT for each habitat)
sbatch 1_gene_centric/1_SingleM/scripts/1.01_subsample_10M_array.sh

# 1.02: SingleM profiling, once per habitat
sbatch --export=ALL,HABITAT=sed 1_gene_centric/1_SingleM/scripts/1.02_singlem_pipe_array.sh
sbatch --export=ALL,HABITAT=sea 1_gene_centric/1_SingleM/scripts/1.02_singlem_pipe_array.sh
```

## 2. Profile parsing <a name="parse"></a>
Directory: [`1_community_profiles/`](1_community_profiles/)

- Combines per-sample SingleM profiles into a single coverage matrix (rows: samples; columns: taxa; absent taxa = 0)
- Platform: Python (pandas)
- Input: directory of `*_prof.tsv` files from section 1 (one run per habitat)
- Output: coverage matrix (CSV), e.g. `singlem_10Mseqs.profsum.RS226.sed.csv`
- Script: [`scripts/1.03_singlem_prof_coverage.py`](scripts/1.03_singlem_prof_coverage.py)

```bash
# 1.03: profile parsing, once per habitat
python 1_gene_centric/1_SingleM/scripts/1.03_singlem_prof_coverage.py \
  --indir <path/to/singlem_output>/sed/R226/prof \
  --out singlem_10Mseqs.profsum.RS226.sed.csv
```

## 3. Ordination and PERMANOVA <a name="ordination"></a>
Directory: [`2_ordination_PERMANOVA/`](2_ordination_PERMANOVA/)

- Builds a phyloseq object (SingleM coverage, taxonomy, and sample metadata); coverage is scaled to estimated prokaryotic read counts using `microbial_fraction`, and non-prokaryotic reads are removed
- Sediment vs seawater: PCA of family-level CLR profiles (Fig. 1d), Shannon diversity (Wilcoxon rank-sum test), and PERMANOVA with dispersion test on Aitchison distances
- Within each habitat: PCA coloured by sampling trip (Fig. 1g-h) and sequential PERMANOVA on Aitchison distances (9,999 permutations): `distance ~ zoning status + sampling trip + GBR sector + reef site` (Fig. 1g-h tables)
- Reef-level PERMANOVA (replicates summed per reef) to test reef-level factors without pseudoreplication
- QC: PCA of sediment samples by sequencing facility
- Platform: R (phyloseq, microViz, vegan)
- Input: `IMOS_10Mseqs_SingleM_FINAL.xlsx` (coverage matrix from section 2, taxonomy, and metadata)
- Output: figures, PERMANOVA tables, and phyloseq objects (`singlem_prok_SEASED.rds`, `singlem_prok_SED.rds`, `singlem_prok_SEA.rds`) used by later analyses
- Script: [`scripts/1.04_singlem_community_analysis.R`](scripts/1.04_singlem_community_analysis.R)

## 4. Abundance, prevalence, and differential abundance <a name="abundance"></a>
Directory: [`3_abundance_DA/`](3_abundance_DA/)

- Mean relative abundance and prevalence of each family in each habitat, identifying families prevalent in > 80% of samples
- Composition of the top 15 prokaryotic classes, with samples grouped by habitat and sampling trip (Fig. 1c)
- Family-level differential abundance between sediment and seawater (seawater as reference), with reef site as a random effect, using MaAsLin2 v1.16.0 and LinDA v0.2.0 on families detected in ≥ 10% of samples
- A family is differentially abundant if both methods give adjusted p ≤ 0.05 in the same direction
- Platform: R (phyloseq, microViz, Maaslin2, MicrobiomeStat)
- Input: `singlem_prok_SEASED.rds` from section 3
- Output: Tables S4 (differential abundance), S5 (sediment), and S6 (seawater); Fig. 1c
- Script: [`scripts/1.05_singleM_differential_abundance.R`](scripts/1.05_singleM_differential_abundance.R)

## 5. Marker gene novelty (appraise) <a name="appraise"></a>
Directory: [`4_appraise_novelty/`](4_appraise_novelty/)

- `singlem pipe` was run on assembled contigs, species-representative MAGs, and 143,614 GTDB R10-RS226 representative genomes to produce mOTU tables
- `singlem appraise` compares read mOTUs with these tables to show how much of the community is lost during assembly and binning, and how novel the mOTUs are compared with GTDB
- Recovery is assessed at species level (exact match) and genus level (`--imperfect --sequence-identity 0.89`), separately for sediment and seawater
- Platform: Bash (SLURM)
- Input: read OTU tables from section 1, assemblies and species-representative MAGs (see [Genome recovery](../../2_genome_centric/1_genome_recovery/)), and GTDB R10-RS226 representative genomes
- Output: proportion of read mOTUs that match GTDB genomes at species level (sediment: 1.4% bacterial, 0.1% archaeal; seawater: 35.7% bacterial, 26.7% archaeal)
- Script: [`scripts/1.06_singlem_appraise.sh`](scripts/1.06_singlem_appraise.sh)

