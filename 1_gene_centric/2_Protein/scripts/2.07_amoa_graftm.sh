#!/bin/bash --login
#SBATCH --job-name=amoa_graftm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o amoa_graftm_%j.o
#SBATCH -e amoa_graftm_%j.e

# Annotation: ammonia monooxygenase (amoA) with GraftM v0.15.1 (package 7.22,
# pmoA/amoA), which also places sequences into functional clades
# (bacterial or archaeal amoA).

set -euo pipefail

CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
GPKG="<path/to/7.22.particulate_methane_monoocygenase_pmoA_31082015.gpkg>"
OUTDIR="<path/to/2_annotation/amoa>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

conda activate singlem   # environment with GraftM v0.15.1

graftM graft --forward ${CATALOGUE} \
    --graftm_package ${GPKG} \
    --input_sequence_type aminoacid \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output_directory ${TMP}/imos_catalogue_amoa

cp -r ${TMP}/imos_catalogue_amoa ${OUTDIR}/
