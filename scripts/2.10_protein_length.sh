#!/bin/bash --login
#SBATCH --job-name=protein_length
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o protein_length_%j.o
#SBATCH -e protein_length_%j.e

# Read mapping, step 2: length of each catalogue protein in nucleotides
# (amino acids x 3), for RPKM. The protein ID is the first word of the header,
# as reported by DIAMOND.

set -euo pipefail

CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
OUTDIR="<path/to/3_read_mapping_RPKM>"

mkdir -p ${OUTDIR}

awk '/^>/ {if (seqlen) print gene "\t" seqlen * 3; split(substr($0, 2), a, " "); gene = a[1]; seqlen = 0; next}
     {seqlen += length($0)}
     END {print gene "\t" seqlen * 3}' ${CATALOGUE} > ${OUTDIR}/comp_protein_length.fixed.tsv

wc -l ${OUTDIR}/comp_protein_length.fixed.tsv
