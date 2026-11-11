#!/bin/bash --login
#SBATCH --job-name=amoa_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o amoa_mags_%A_%a.o
#SBATCH -e amoa_mags_%A_%a.e
#SBATCH --array=1-2077

# MAG annotation, step 5: ammonia monooxygenase (amoA) with GraftM v0.15.1
# (package 7.22, pmoA/amoA), which also places sequences into functional clades
# (bacterial or archaeal amoA).

set -euo pipefail

MAG_LIST="<path/to/mag_list.txt>"   # one MAG FASTA path per line (1,201 sediment and 876 seawater MAGs)
GPKG="<path/to/7.22.particulate_methane_monoocygenase_pmoA_31082015.gpkg>"
OUTDIR="<path/to/5_MAG_annotation/5_amoA>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

mag=$(sort ${MAG_LIST} | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${mag} .fna)
echo "MAG: ${name}"

conda activate singlem   # environment with GraftM v0.15.1

graftM graft --forward ${mag} \
    --graftm_package ${GPKG} \
    --input_sequence_type nucleotide \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output_directory ${TMP}/${name}

cp -r ${TMP}/${name} ${OUTDIR}/
