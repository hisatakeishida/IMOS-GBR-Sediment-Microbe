#!/bin/bash --login
#SBATCH --job-name=dbcan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o dbcan_%A_%a.o
#SBATCH -e dbcan_%A_%a.e
#SBATCH --array=1-393

# Annotation: CAZyme families with run_dbcan v5.0.1 CAZyme_annotation, one task
# per DRAM chunk from 2.3 (same protein sequences, same <start>_<end>_ prefix).

set -euo pipefail

DRAM_DIR="<path/to/2_annotation/dram>"
DBCAN_DB="<path/to/run_dbcan/db>"
OUTDIR="<path/to/2_annotation/dbcan>"

chunk=$(ls -d ${DRAM_DIR}/*__annotation | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${chunk} __annotation)
echo "Chunk: ${name}"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

conda activate run_dbcan_env   # run_dbcan v5.0.1

run_dbcan CAZyme_annotation \
    --input_raw_data ${chunk}/genes.faa \
    --mode protein \
    --output_dir ${TMP}/${name}_dbcan \
    --db_dir ${DBCAN_DB} \
    --threads ${SLURM_CPUS_PER_TASK}

cp ${TMP}/${name}_dbcan/overview.tsv ${OUTDIR}/${name}_dbcan.tsv
