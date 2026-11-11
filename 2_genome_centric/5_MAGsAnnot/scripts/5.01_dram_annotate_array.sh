#!/bin/bash --login
#SBATCH --job-name=dram_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o dram_mags_%A_%a.o
#SBATCH -e dram_mags_%A_%a.e
#SBATCH --array=1-2077

# MAG annotation, step 1: annotate each MAG with DRAM v1.5 (KOfam, Pfam, and
# MEROPS databases). DRAM also predicts genes (genes.faa, genes.gff), which
# are used by 5.02 to 5.04.

set -euo pipefail

MAG_LIST="<path/to/mag_list.txt>"   # one MAG FASTA path per line (1,201 sediment and 876 seawater MAGs)
OUTDIR="<path/to/5_MAG_annotation/1_DRAM>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

mag=$(sort ${MAG_LIST} | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${mag} .fna)
echo "MAG: ${name}"

module load dram/1.5   # DRAM v1.5

DRAM.py annotate -i ${mag} -o ${TMP}/${name}__annotation --threads ${SLURM_CPUS_PER_TASK}

cp -r ${TMP}/${name}__annotation ${OUTDIR}/
