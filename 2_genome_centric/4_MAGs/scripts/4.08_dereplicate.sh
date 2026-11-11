#!/bin/bash --login
#SBATCH --job-name=coverm_dereplicate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o coverm_dereplicate_%j.o
#SBATCH -e coverm_dereplicate_%j.e

# MAG recovery, step 8: dereplicate MAGs at 95% ANI (approximately species
# level) with CoverM cluster v0.7.0, choosing representatives with the
# Parks2020_reduced quality formula (CheckM2 values).
#
# Dereplication was done in two stages:
#   Stage 1: single-sample assembly MAGs
#   Stage 2: stage 1 representatives together with co-assembly MAGs
# "Eukaryotic MAGs" were then removed, giving 1,201 representatives.
#
#   sbatch --export=ALL,STAGE=1 scripts/4.8_dereplicate.sh
#   sbatch --export=ALL,STAGE=2 scripts/4.8_dereplicate.sh

set -euo pipefail

STAGE=${STAGE:?Set STAGE=1 or STAGE=2}
MAG_DIR="<path/to/4_MAG_recovery/5_MAG_QC/mags>"                                   # output of 4.7
QUALITY="<path/to/4_MAG_recovery/5_MAG_QC/mags_Q50.checkm2_quality_report.tsv>"     # output of 4.7
OUTDIR="<path/to/4_MAG_recovery/6_dereplication>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${TMP}/genomes

case ${STAGE} in
    1)  # Single-sample assembly MAGs (IDs not starting with "coassembly_")
        find ${MAG_DIR} -name "*.fna" ! -name "coassembly_*" -exec cp {} ${TMP}/genomes/ \;
        prefix=stage1_single_assembly_ANI95_AF15 ;;
    2)  # Stage 1 representatives and co-assembly MAGs
        while read -r rep; do cp ${rep} ${TMP}/genomes/; done < ${OUTDIR}/stage1_single_assembly_ANI95_AF15.rep.list.txt
        cp ${MAG_DIR}/coassembly_*.fna ${TMP}/genomes/
        prefix=stage2_single_coassembly_ANI95_AF15 ;;
    *)  echo "Unknown STAGE: ${STAGE}" >&2; exit 1 ;;
esac
ls ${TMP}/genomes | wc -l

conda activate coverm   # CoverM v0.7.0

coverm cluster \
    --genome-fasta-directory ${TMP}/genomes -x fna \
    --ani 95 \
    --min-aligned-fraction 15 \
    --quality-formula Parks2020_reduced \
    --checkm2-quality-report ${QUALITY} \
    -t ${SLURM_CPUS_PER_TASK} \
    --output-cluster-definition ${OUTDIR}/${prefix}.cluster.def.txt \
    --output-representative-list ${OUTDIR}/${prefix}.rep.list.txt

wc -l ${OUTDIR}/${prefix}.rep.list.txt

# Stage 2 only: copy representatives, then remove eukaryotic MAGs (TODO: method)
if [ "${STAGE}" = "2" ]; then
    mkdir -p ${OUTDIR}/representatives
    while read -r rep; do cp ${rep} ${OUTDIR}/representatives/; done < ${OUTDIR}/${prefix}.rep.list.txt
fi
