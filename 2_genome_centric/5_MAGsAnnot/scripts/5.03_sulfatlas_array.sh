#!/bin/bash --login
#SBATCH --job-name=sulfatlas_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o sulfatlas_mags_%A_%a.o
#SBATCH -e sulfatlas_mags_%A_%a.e
#SBATCH --array=1-2077

# MAG annotation, step 3: sulphatases with DIAMOND blastp v2.1.0 against
# SulfAtlas v2.3.1 (E-value 1e-20, best hit per protein).

set -euo pipefail

DRAM_DIR="<path/to/5_MAG_annotation/1_DRAM>"          # output of 5.01
SULFATLAS="<path/to/sulfatlas_v2.3.1.dmnd>"
OUTDIR="<path/to/5_MAG_annotation/3_sulphatases>"

mkdir -p ${OUTDIR}

dir=$(ls -d ${DRAM_DIR}/*__annotation | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${dir} __annotation)
echo "MAG: ${name}"

module load diamond/2.1.0-gcc-11.3.0

diamond blastp --threads ${SLURM_CPUS_PER_TASK} \
    -d ${SULFATLAS} -q ${dir}/genes.faa \
    --evalue 1e-20 --max-target-seqs 1 \
    -o ${OUTDIR}/${name}.tsv
