#!/bin/bash --login
#SBATCH --job-name=coverm_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o coverm_mags_%A_%a.o
#SBATCH -e coverm_mags_%A_%a.e
#SBATCH --array=1-<number_of_samples>

# MAG recovery, step 12: map each habitat's reads to its own MAGs with CoverM
# genome (sediment reads to 1,201 sediment MAGs; seawater reads to 876
# seawater MAGs), excluding alignments with < 95% identity or < 75% aligned
# length. Gives MAG-based community profiles and the proportion of reads
# mapped to MAGs (reported in the log).
#
#   sbatch --export=ALL,HABITAT=sed,MODE=count --array=1-191 2_genome_centric/4_MAGs/scripts/4.12_coverm_mags_array.sh
#   MODE: count or relative_abundance

set -euo pipefail

HABITAT=${HABITAT:?Set HABITAT=sed or HABITAT=sea}
MODE=${MODE:-count}
READS_DIR="<path/to/adapter_removed_reads>/${HABITAT}"   # <sample>_R{1,2}_allRunsCombined.fastq.gz
MAG_DIR="<path/to/MAGs>/${HABITAT}"                       # 1,201 sediment or 876 seawater MAGs
OUTDIR="<path/to/4_MAG_recovery/10_read_mapping>/${HABITAT}/${MODE}"

mkdir -p ${OUTDIR}

R1=$(ls ${READS_DIR}/*_R1_allRunsCombined*.fastq.gz | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename "${R1%%_R1_*}")
R2=$(ls ${READS_DIR}/${name}_R2_*.fastq.gz)
echo "Habitat: ${HABITAT}; sample: ${name}"

conda activate coverm   # CoverM v0.7.0

coverm genome \
    --genome-fasta-directory ${MAG_DIR} -x fna \
    -1 ${R1} -2 ${R2} \
    -m ${MODE} \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --min-covered-fraction 0 \
    -t ${SLURM_CPUS_PER_TASK} \
    -o ${OUTDIR}/${name}.readperid.95.readali.75.${MODE}.txt
