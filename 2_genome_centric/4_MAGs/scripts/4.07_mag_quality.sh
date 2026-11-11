#!/bin/bash --login
#SBATCH --job-name=mag_quality
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o mag_quality_%j.o
#SBATCH -e mag_quality_%j.e

# MAG recovery, step 7: combine CheckM v1.2.2 and CheckM2 v1.0.2 results from
# Aviary recover (4.4 and 4.6), and keep bins with quality score
# (completeness - 3 x contamination) >= 50% by either tool.
#
# Bin IDs are prefixed with their source (<assembly or coassembly>___<bin>).
# Output: merged quality table, table of retained MAGs (2,698), and their
# FASTA files.

set -euo pipefail

BIN_DIRS=("<path/to/4_MAG_recovery/3_binning>" "<path/to/4_MAG_recovery/4_coassembly/bins>")   # outputs of 4.4 and 4.6
OUTDIR="<path/to/4_MAG_recovery/5_MAG_QC>"
MIN_QUALITY=50

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}/mags

c1_all=${TMP}/checkm1.tsv
c2_all=${TMP}/checkm2.tsv
> ${c1_all}
> ${c2_all}

# ---- Collect CheckM and CheckM2 results with source-prefixed bin IDs ----------
for dir in "${BIN_DIRS[@]}"; do
    for source_dir in ${dir}/*/; do
        source=$(basename ${source_dir})
        c1=${source_dir}/bins/checkm1_output.detailed.out
        c2=${source_dir}/bins/checkm2_output/quality_report.tsv
        [ -f ${c1} ] && [ -f ${c2} ] || { echo "${source}: no CheckM results (no bins?)"; continue; }
        tail -n +2 ${c1} | awk -F'\t' -v OFS='\t' -v p="${source}___" '{$1 = p $1; print}' >> ${c1_all}
        tail -n +2 ${c2} | awk -F'\t' -v OFS='\t' -v p="${source}___" '{$1 = p $1; print}' >> ${c2_all}
    done
done

# ---- Merge and calculate quality ------------------------------------------------
# CheckM (tab table, -o 2): $6 completeness, $7 contamination
# CheckM2 (quality_report.tsv): $2 completeness, $3 contamination
merged=${OUTDIR}/bins_checkm1and2.tsv
awk -F'\t' -v OFS='\t' '
    NR == FNR {c2_comp[$1] = $2; c2_cont[$1] = $3; next}
    {
        q1 = $6 - 3 * $7
        q2 = ($1 in c2_comp) ? c2_comp[$1] - 3 * c2_cont[$1] : "NA"
        split($1, s, "___")
        print $1, s[1], $6, $7, c2_comp[$1], c2_cont[$1], q1, q2
    }' ${c2_all} ${c1_all} | \
    sed '1i Bin_ID\tBin_Source\tCompleteness_C1\tContamination_C1\tCompleteness_C2\tContamination_C2\tQuality_C1\tQuality_C2' \
    > ${merged}

awk -F'\t' -v q=${MIN_QUALITY} 'NR == 1 || $7 >= q || ($8 != "NA" && $8 >= q)' ${merged} > ${merged%.tsv}.Q${MIN_QUALITY}.tsv
echo "Retained MAGs: $(( $(wc -l < ${merged%.tsv}.Q${MIN_QUALITY}.tsv) - 1 ))"

# CheckM2 table for the retained MAGs (used for dereplication in 4.8)
head -n 1 $(ls ${BIN_DIRS[0]}/*/bins/checkm2_output/quality_report.tsv | head -1) > ${OUTDIR}/mags_Q${MIN_QUALITY}.checkm2_quality_report.tsv
awk -F'\t' 'NR == FNR {if (FNR > 1) keep[$1]; next} $1 in keep' \
    ${merged%.tsv}.Q${MIN_QUALITY}.tsv ${c2_all} >> ${OUTDIR}/mags_Q${MIN_QUALITY}.checkm2_quality_report.tsv

# ---- Copy retained MAGs -------------------------------------------------------------
# Co-binning outputs archive data/ (4.6); final bins are read from bins/final_bins
tail -n +2 ${merged%.tsv}.Q${MIN_QUALITY}.tsv | cut -f 1 | while read -r bin_id; do
    source=${bin_id%%___*}
    bin=${bin_id#*___}
    for dir in "${BIN_DIRS[@]}"; do
        f=${dir}/${source}/bins/final_bins/${bin}.fna
        [ -f ${f} ] && cp ${f} ${OUTDIR}/mags/${bin_id}.fna
    done
done
ls ${OUTDIR}/mags | wc -l
