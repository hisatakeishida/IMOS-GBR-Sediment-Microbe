#!/bin/bash --login
#SBATCH --job-name=diamond_makedb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o diamond_makedb_%j.o
#SBATCH -e diamond_makedb_%j.e

# Read mapping, step 1: DIAMOND database of the non-redundant protein catalogue.

set -euo pipefail

CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
OUTDIR="<path/to/3_read_mapping_RPKM/dmnd>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

module load diamond/2.1.0-gcc-11.3.0

name=$(basename ${CATALOGUE} .fasta)
diamond makedb --in ${CATALOGUE} -d ${TMP}/${name} --threads ${SLURM_CPUS_PER_TASK}

cp ${TMP}/${name}.dmnd ${OUTDIR}/
