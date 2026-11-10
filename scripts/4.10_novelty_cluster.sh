#!/bin/bash --login
#SBATCH --job-name=coverm_novelty
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500G
#SBATCH --time=72:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o coverm_novelty_%j.o
#SBATCH -e coverm_novelty_%j.e

# MAG recovery, step 10: novelty of sediment and seawater MAGs against public
# genome collections, by clustering all genomes together at 95% ANI with
# CoverM cluster v0.7.0 (minimum aligned fraction 15%; Parks2020_reduced with
# CheckM2 values), then counting the genomes from each collection in each
# cluster (input for the UpSet plot in 4.11).
#   Public collections: Ocean Microbiomics Database (OMD), Reef Microbiomics
#   Database (RMD), Genomes from Earth's Microbiomes (GEMs), shallow-water
#   coastal marine sediments (Silva-Solar et al., 2024; Viver et al., 2025),
#   South China Sea reef sediments (Dong et al., 2023)
#
# QUALITY: one CheckM2 quality_report.tsv covering all genomes below
# (genome name in column 1, without extension).

set -euo pipefail

# Collection name -> genome directory (names become the columns of the summary)
declare -A SOURCES=(
    [gbr_sediment]="<path/to/4_MAG_recovery/6_dereplication/representatives>"   # 1,201 sediment MAGs
    [gbr_seawater]="<path/to/seawater_MAGs>"                                     # 876 seawater MAGs
    [omd]="<path/to/OMD_genomes>"
    [rmd]="<path/to/RMD_genomes>"
    [gems]="<path/to/GEMs_genomes>"
    [coastal_sediment]="<path/to/coastal_marine_sediment_MAGs>"
    [china_sea_reef_sediment]="<path/to/South_China_Sea_reef_sediment_MAGs>"
)
QUALITY="<path/to/all_genomes.checkm2_quality_report.tsv>"
OUTDIR="<path/to/4_MAG_recovery/9_novelty>"
PREFIX=sed_sea_public_ANI95_AF15

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${TMP}/genomes

# ---- Collect genomes (as .fna) and record their collection ---------------------------
genome_source=${OUTDIR}/genome_source.tsv
> ${genome_source}
for src in "${!SOURCES[@]}"; do
    find ${SOURCES[${src}]} -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fna.gz" -o -name "*.fa.gz" \) | \
    while read -r f; do
        base=$(basename ${f}); base=${base%.gz}; base=${base%.fna}; base=${base%.fa}
        if [[ ${f} == *.gz ]]; then zcat ${f} > ${TMP}/genomes/${base}.fna; else cp ${f} ${TMP}/genomes/${base}.fna; fi
        echo -e "${base}\t${src}" >> ${genome_source}
    done
done
cut -f 2 ${genome_source} | sort | uniq -c

# ---- Cluster --------------------------------------------------------------------------
conda activate coverm   # CoverM v0.7.0

coverm cluster \
    --genome-fasta-directory ${TMP}/genomes -x fna \
    --ani 95 \
    --min-aligned-fraction 15 \
    --quality-formula Parks2020_reduced \
    --checkm2-quality-report ${QUALITY} \
    -t ${SLURM_CPUS_PER_TASK} \
    --output-cluster-definition ${OUTDIR}/${PREFIX}.cluster.def.txt \
    --output-representative-list ${OUTDIR}/${PREFIX}.rep.list.txt

# ---- Genomes per collection in each cluster (rows: clusters; columns: collections) -------
sources=$(printf "%s\n" "${!SOURCES[@]}" | sort | paste -sd " ")
awk -F'\t' -v OFS='\t' -v sources="${sources}" '
    function strip(p) {sub(/.*\//, "", p); sub(/\.fna$/, "", p); return p}
    NR == FNR {src[$1] = $2; next}
    {rep = strip($1); n[rep, src[strip($2)]]++; reps[rep]}
    END {
        ns = split(sources, s, " ")
        header = "Rep"; for (i = 1; i <= ns; i++) header = header OFS s[i]; print header
        for (r in reps) {line = r; for (i = 1; i <= ns; i++) line = line OFS ((r, s[i]) in n ? n[r, s[i]] : 0); print line}
    }' ${genome_source} ${OUTDIR}/${PREFIX}.cluster.def.txt > ${OUTDIR}/${PREFIX}.cluster_summary.tsv

# Sediment clusters without any genome from another collection
awk -F'\t' 'NR == 1 {for (i = 2; i <= NF; i++) col[$i] = i; next}
    $col["gbr_sediment"] > 0 {other = 0; for (c in col) if (c != "gbr_sediment") other += $col[c]; if (other == 0) novel++; sed++}
    END {print "Sediment clusters: " sed "; sediment only: " novel}' ${OUTDIR}/${PREFIX}.cluster_summary.tsv