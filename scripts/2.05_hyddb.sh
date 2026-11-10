#!/bin/bash --login
#SBATCH --job-name=hyddb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o hyddb_%j.o
#SBATCH -e hyddb_%j.e

# Annotation: hydrogenases with DIAMOND blastp v2.1.0 against HydDB, keeping
# hits above class-specific identity thresholds:
#   [NiFe] Group 4 > 50%; other [NiFe] > 30%; [FeFe] > 45%; [Fe] > 50%

set -euo pipefail

CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
HYDDB="<path/to/HydDB/hyddb.dmnd>"
OUTDIR="<path/to/2_annotation/hyddb>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

module load diamond/2.1.0-gcc-11.3.0

diamond blastp --threads ${SLURM_CPUS_PER_TASK} \
    -q ${CATALOGUE} -d ${HYDDB} \
    -o ${TMP}/hyddb_hits.tsv \
    --max-target-seqs 1 \
    --header simple \
    --outfmt 6 qseqid sseqid pident evalue bitscore full_sseq

awk 'BEGIN {IGNORECASE = 1; print "query_id\tclosest_hydrogenase_group"}
NR > 1 {
    split($2, arr, /\|/)
    class = arr[3]
    id = $3 + 0

    # Group 4 first, with a stricter threshold; skip other checks
    if (class ~ /\[NiFe\]_Group_4([a-i])?$/) {
        if (id > 50) print $1 "\t" class
        next
    }
    if      (class ~ /\[NiFe\]/ && id > 30) print $1 "\t" class
    else if (class ~ /\[FeFe\]/ && id > 45) print $1 "\t" class
    else if (class ~ /\[Fe\]/   && id > 50) print $1 "\t" class
}' ${TMP}/hyddb_hits.tsv > ${TMP}/hyddb_classification_filtered.tsv

cp ${TMP}/hyddb_hits.tsv ${TMP}/hyddb_classification_filtered.tsv ${OUTDIR}/
