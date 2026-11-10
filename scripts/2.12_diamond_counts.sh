#!/bin/bash --login
#SBATCH --job-name=diamond_counts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o diamond_counts_%j.o
#SBATCH -e diamond_counts_%j.e

# Read mapping, step 4: count reads per protein (catalogue) and per marker gene
# (SingleM) from the DIAMOND hit tables of 2.11, in long format
# (Sample, GeneID, Count), for both habitats combined.
# Sample names must match the row names of the sample metadata.

set -euo pipefail

DIAMOND_DIR="<path/to/3_read_mapping_RPKM/diamond>"   # {sed,sea}/{catalogue,singlem}/ from 2.11
OUTDIR="<path/to/3_read_mapping_RPKM>"
SUFFIX=".adprem.140bp"   # removed from file names to give sample names

mkdir -p ${OUTDIR}

count_hits() {   # $1: catalogue or singlem; $2: file extension
    local db=$1 ext=$2
    local out=${OUTDIR}/diamond_long.${db}.tsv
    echo -e "Sample\tGeneID\tCount" > ${out}
    for f in ${DIAMOND_DIR}/{sed,sea}/${db}/*${ext}; do
        sample=$(basename ${f} ${ext})
        sample=${sample%${SUFFIX}}
        cut -f 2 ${f} | sort | uniq -c | awk -v s="${sample}" 'BEGIN {OFS = "\t"} {print s, $2, $1}' >> ${out}
    done
}

count_hits catalogue .tsv
count_hits singlem .singlem.tsv
