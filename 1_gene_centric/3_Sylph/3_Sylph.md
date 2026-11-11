# k-mer-based taxonomic profiling with sylph

**Aim:** profile sediment and seawater communities against GTDB R10-RS226 genomes with [sylph](https://github.com/bluenote-1577/sylph), as a k-mer-based comparison with the mOTU profiles from [SingleM](../1_SingleM/1_SingleM.md).

This is a supporting analysis and is not reported in the manuscript.

Scripts for this analysis are in the [`scripts/`](scripts/) folder (`3.01` and `3.02`).

| Step | Script | Description |
|---|---|---|
| 3.01 | [`scripts/3.01_sylph_profile_array.sh`](scripts/3.01_sylph_profile_array.sh) | `sylph profile` (with `--estimate-read-counts`) against the GTDB R10-RS226 database (`gtdb-r226-c200-dbv1.syldb`), then `sylph-tax taxprof` for GTDB taxonomy, one task per sample |
| 3.02 | [`scripts/3.02_sylph_merge.sh`](scripts/3.02_sylph_merge.sh) | `sylph-tax merge` of per-sample profiles into sequence abundance and relative abundance tables |

## Sylph profiles
Directory: [`3_Sylph/`](./)

- Platform: Bash (SLURM arrays, one task per sample)
- Input: paired-end reads subsampled to 10 million read pairs per sample (see [SingleM](../1_SingleM/1_SingleM.md), step 1.01)
- Output: per-sample profiles (`1_count/`), taxonomic profiles (`2_sylphmpa/`), and merged tables (`3_summary/`) for each habitat

```bash
# One-off setup
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb
sylph-tax download --download-to <path/to/sylph-tax_db>

# 3.01: profile, once per habitat
sbatch --export=ALL,HABITAT=sed --array=1-<n_sed> 1_gene_centric/3_Sylph/scripts/3.01_sylph_profile_array.sh
sbatch --export=ALL,HABITAT=sea --array=1-<n_sea> 1_gene_centric/3_Sylph/scripts/3.01_sylph_profile_array.sh

# 3.02: merge, once per habitat (after 3.01)
sbatch --export=ALL,HABITAT=sed 1_gene_centric/3_Sylph/scripts/3.02_sylph_merge.sh
sbatch --export=ALL,HABITAT=sea 1_gene_centric/3_Sylph/scripts/3.02_sylph_merge.sh
```

## Summary
Percentage of reads detected in the GTDB R10-RS226 database by `sylph profile` (from the "detected in database by profile" line of each job's log; mean ± SD across samples), for two read sets: paired-end reads subsampled to 10 M read pairs, and R1 reads ≥ 140 bp (the read set used for [protein read mapping](../2_Protein/2_Protein.md#mapping)).

| Reads | Habitat | Samples | Reads detected in database (%) |
|---|---|---|---|
| Paired-end, 10 M read pairs | Sediment | 189 | 0.85 ± 0.65 |
| Paired-end, 10 M read pairs | Seawater | 189 | 49.86 ± 10.15 |
| R1 ≥ 140 bp | Sediment | 191 | 0.88 ± 0.56 |
| R1 ≥ 140 bp | Seawater | 182 | 48.34 ± 9.51 |

With either read set, less than 1% of sediment reads, compared with about half of seawater reads, could be assigned to GTDB genomes, consistent with the much lower representation of sediment mOTUs in GTDB found by [SingleM appraise](../1_SingleM/1_SingleM.md#appraise).

Because sylph only detects genomes in its reference database, whereas SingleM recovers marker gene OTUs directly from reads (including lineages without reference genomes), SingleM is better suited to describing the overall diversity of poorly characterised habitats such as reef sediments.

```bash
# Mean and SD across samples, from the 3.01 logs (run per habitat)
for f in sylph_<habitat>*.e; do
    grep "detected in database by profile" ${f}
done | awk '{x = $6; sum += x; sumsq += x * x; n++}
    END {mean = sum / n; printf "Count: %d\nMean: %.6f\nSD: %.6f\n", n, mean, sqrt(sumsq / n - mean * mean)}'
```
