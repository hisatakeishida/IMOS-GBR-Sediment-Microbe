#!/bin/bash --login
#SBATCH --job-name=dbcan_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o dbcan_mags_%A_%a.o
#SBATCH -e dbcan_mags_%A_%a.e
#SBATCH --array=1-2077

# MAG annotation, step 2: CAZyme families, CAZyme gene clusters (CGCs), and
# CGC substrates (similarity to PULs in dbCAN-PUL) with run_dbcan v5.0.1
# easy_substrate, using the proteins and gene coordinates predicted by DRAM.

set -euo pipefail

DRAM_DIR="<path/to/5_MAG_annotation/1_DRAM>"   # output of 5.01
DBCAN_DB="<path/to/run_dbcan/db>"
OUTDIR="<path/to/5_MAG_annotation/2_CAZymes_CGCs>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

dir=$(ls -d ${DRAM_DIR}/*__annotation | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${dir} __annotation)
echo "MAG: ${name}"

conda activate run_dbcan_env   # run_dbcan v5.0.1

run_dbcan easy_substrate \
    --input_raw_data ${dir}/genes.faa \
    --mode protein \
    --input_gff ${dir}/genes.gff \
    --gff_type prodigal \
    --db_dir ${DBCAN_DB} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output_dir ${TMP}/${name}

cp -r ${TMP}/${name} ${OUTDIR}/
