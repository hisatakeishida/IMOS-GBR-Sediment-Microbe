#!/bin/bash --login
#SBATCH --job-name=sylph_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o sylph_merge_%j.o
#SBATCH -e sylph_merge_%j.e

# Merge per-sample sylph taxonomic profiles (3.1) into one table per habitat,
# for sequence abundance and relative (taxonomic) abundance.
#
#   sbatch --export=ALL,HABITAT=sed scripts/3.2_sylph_merge.sh
#   sbatch --export=ALL,HABITAT=sea scripts/3.2_sylph_merge.sh

set -euo pipefail

HABITAT=${HABITAT:?Set HABITAT=sed or HABITAT=sea}
OUTDIR="<path/to/3_Sylph>/${HABITAT}"   # output of 3.1

mkdir -p ${OUTDIR}/3_summary

conda activate sylph_env

for column in sequence_abundance relative_abundance; do
    sylph-tax merge ${OUTDIR}/2_sylphmpa/*.sylphmpa \
        --column ${column} \
        -o ${OUTDIR}/3_summary/merged_${column}_${HABITAT}.tsv
done