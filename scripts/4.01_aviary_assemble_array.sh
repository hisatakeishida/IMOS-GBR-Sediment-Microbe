#!/bin/bash --login
#SBATCH --job-name=aviary_assemble
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=<memory, e.g. 600G>
#SBATCH --time=144:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o aviary_assemble_%A_%a.o
#SBATCH -e aviary_assemble_%A_%a.e
#SBATCH --array=1-191

# MAG recovery, step 1: single-sample assembly with Aviary assemble v0.9.1.
#   Hybrid (17 samples): metaFlye v2.9.4 (ONT) with metaSPAdes v3.15.4 (Illumina)
#   Illumina-only (174 samples): metaSPAdes v3.15.4
# Reads were adapter-trimmed beforehand (fastp v0.23.4, porechop v0.2.4), so
# Aviary QC is skipped.
#
# SAMPLES: tab-separated, one sample per line:
#   <sample>  <R1.fastq.gz>  <R2.fastq.gz>  [<ONT.fastq.gz>]
# The fourth column is given only for the 17 hybrid samples.
# Memory needed ranged from 600 to 1,700 GB for the deepest samples.

set -euo pipefail

SAMPLES="<path/to/assembly_samples.tsv>"
OUTDIR="<path/to/4_MAG_recovery/1_assembly>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

IFS=$'\t' read -r sample R1 R2 ONT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES})
echo "Sample: ${sample}"

LONGREADS=""
if [ -n "${ONT:-}" ]; then
    LONGREADS="--longreads ${ONT} --long_read_type ont_hq"
    echo "Hybrid assembly with ${ONT}"
fi

conda activate aviary   # Aviary v0.9.1

aviary assemble \
    -1 ${R1} -2 ${R2} ${LONGREADS} \
    --skip-qc True \
    --max-threads ${SLURM_CPUS_PER_TASK} --n-cores ${SLURM_CPUS_PER_TASK} \
    --max_memory $(( SLURM_MEM_PER_NODE / 1024 )) \
    --conda-frontend conda \
    --tmpdir ${TMP} \
    --output ${TMP}/${sample}

cp -r ${TMP}/${sample} ${OUTDIR}/
