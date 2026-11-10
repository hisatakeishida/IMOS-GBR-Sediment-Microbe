#!/bin/bash --login
#SBATCH --job-name=combine_annot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=1:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o combine_annot_%j.o
#SBATCH -e combine_annot_%j.e

# Annotation: combine per-chunk and per-database results into one table per
# annotation (gene ID -> annotation), using catalogue protein IDs.
#   1. KOs (DRAM)
#   2. CAZyme families (run_dbcan)
#   3. Hydrogenases (HydDB), with KO added
#   4. Metabolic marker genes (best hit per protein across databases)

set -euo pipefail

DRAM_DIR="<path/to/2_annotation/dram>"            # output of 2.3
DBCAN_DIR="<path/to/2_annotation/dbcan>"          # output of 2.4
HYDDB_DIR="<path/to/2_annotation/hyddb>"          # output of 2.5
MARKER_DIR="<path/to/2_annotation/marker_genes>"  # output of 2.6
OUTDIR="<path/to/2_annotation/summary>"

mkdir -p ${OUTDIR}

# DRAM (and run_dbcan, which used DRAM's genes.faa) prefix gene IDs with the
# chunk name (<start>_<end>_); strip it to recover catalogue protein IDs.
strip_prefix='BEGIN {FS = OFS = "\t"} {sub(/^[^_]*_[^_]*_/, "", $1); print}'

# ---- 1. KOs -------------------------------------------------------------------
echo -e "GeneID\tko_id" > ${OUTDIR}/ko_annot.tsv
echo -e "ko_id\tko_description" > ${OUTDIR}/ko_def.tsv

for d in ${DRAM_DIR}/*__annotation; do
    tail -n +2 ${d}/annotations.tsv | cut -f 1,4 | awk "${strip_prefix}" >> ${OUTDIR}/ko_annot.tsv
    tail -n +2 ${d}/annotations.tsv | cut -f 4,5
done | sort -u >> ${OUTDIR}/ko_def.tsv

# ---- 2. CAZyme families -------------------------------------------------------
# overview.tsv column 7: recommended CAZyme family ("-" if none). Families within
# a protein are deduplicated and sorted.
echo -e "GeneID\tcazy_id" > ${OUTDIR}/dbcan_annot.tsv

for f in ${DBCAN_DIR}/*_dbcan.tsv; do
    tail -n +2 ${f} | cut -f 1,7
done | awk "${strip_prefix}" | awk -F'\t' '$2 != "-"' | \
gawk -F'\t' -v OFS='\t' '{
    n = split($2, a, "|"); delete u; k = 0
    for (i = 1; i <= n; i++) if (!(a[i] in u)) { u[a[i]]; keys[++k] = a[i] }
    asort(keys)
    $2 = keys[1]; for (i = 2; i <= k; i++) $2 = $2 "|" keys[i]
    delete keys; print
}' >> ${OUTDIR}/dbcan_annot.tsv

# ---- 3. Hydrogenases, with KO -------------------------------------------------
echo -e "GeneID\thyddb_id\tko_id" > ${OUTDIR}/hyddb_annot.tsv
awk -F'\t' -v OFS='\t' 'NR == FNR {ko[$1] = $2; next} FNR > 1 {print $1, $2, ko[$1]}' \
    ${OUTDIR}/ko_annot.tsv ${HYDDB_DIR}/hyddb_classification_filtered.tsv >> ${OUTDIR}/hyddb_annot.tsv

# ---- 4. Metabolic marker genes ------------------------------------------------
# Best hit per protein across all marker databases, by E-value (column 7)
echo -e "GeneID\tmarker_gene" > ${OUTDIR}/marker_genes_annot.tsv

for f in ${MARKER_DIR}/*.blastp.tsv; do
    gene=$(basename ${f} .blastp.tsv)
    awk -v gene="${gene}" 'BEGIN {OFS = "\t"} {print $1, gene, $7}' ${f}
done | awk 'BEGIN {OFS = "\t"} {
    if (!($1 in best_eval) || $3 < best_eval[$1]) {
        best_eval[$1] = $3
        best_hit[$1]  = $1 "\t" $2
    }
} END {for (q in best_hit) print best_hit[q]}' >> ${OUTDIR}/marker_genes_annot.tsv
