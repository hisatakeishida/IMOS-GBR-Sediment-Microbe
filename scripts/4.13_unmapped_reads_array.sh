#!/bin/bash --login
#SBATCH --job-name=mag_unmapped_reads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o mag_unmapped_reads_%A_%a.o
#SBATCH -e mag_unmapped_reads_%A_%a.e
#SBATCH --array=1-191

# MAG recovery, step 13: reads that do not map to the sediment MAGs
# (alignments with >= 95% identity over >= 75% of the read are counted as
# mapped), for mapping to unbinned contigs in 4.14.

set -euo pipefail

READS_DIR="<path/to/adapter_removed_reads>/sed"
MAG_CONCAT="<path/to/sediment_MAGs_concatenated.fna>"   # all sediment MAGs in one FASTA, unique contig IDs
OUTDIR="<path/to/4_MAG_recovery/11_unbinned_contigs/unmapped_reads>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

R1=$(ls ${READS_DIR}/*_R1_allRunsCombined.fastq.gz | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${R1} _R1_allRunsCombined.fastq.gz)
R2=${READS_DIR}/${name}_R2_allRunsCombined.fastq.gz
echo "Sample: ${name}"

conda activate coverm   # CoverM v0.7.0 (with samtools)

coverm make --reference ${MAG_CONCAT} -1 ${R1} -2 ${R2} \
    -t ${SLURM_CPUS_PER_TASK} -o ${TMP}/${name}

coverm filter --bam-files ${TMP}/${name}/*.bam \
    --inverse \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    -t ${SLURM_CPUS_PER_TASK} \
    -o ${TMP}/${name}_unmapped.bam

samtools sort -n -@ ${SLURM_CPUS_PER_TASK} ${TMP}/${name}_unmapped.bam -o ${TMP}/${name}_unmapped.sorted.bam
samtools fastq -@ ${SLURM_CPUS_PER_TASK} ${TMP}/${name}_unmapped.sorted.bam \
    -1 ${OUTDIR}/${name}_unmapped_R1_allRunsCombined.fastq.gz \
    -2 ${OUTDIR}/${name}_unmapped_R2_allRunsCombined.fastq.gz \
    -0 /dev/null -s /dev/null -n
