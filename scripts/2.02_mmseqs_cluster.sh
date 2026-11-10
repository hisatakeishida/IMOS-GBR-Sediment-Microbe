#!/bin/bash --login
#SBATCH --job-name=mmseqs_cluster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500G
#SBATCH --time=72:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o mmseqs_cluster_%j.o
#SBATCH -e mmseqs_cluster_%j.e

# Gene catalogue, step 2: cluster complete ORFs from both habitats into a
# non-redundant protein catalogue (MMseqs2 easy-linclust, 90% identity, 80% coverage).
#   60,940,012 complete ORFs -> 39,225,676 non-redundant proteins

set -euo pipefail

INDIR="<path/to/1_gene_catalogue>"   # output of 2.1
OUTDIR="<path/to/1_gene_catalogue>"
PREFIX=imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
cd ${TMP}

cat ${INDIR}/imos_sediment_concatenated.protein.translations.renamed.onlycomp.faa \
    ${INDIR}/imos_seawater_concatenated.protein.translations.renamed.onlycomp.faa \
    > ${TMP}/sediment_seawater.onlycomp.faa

module load mmseqs2/15-6f452   # MMseqs2 v15.6f452

mmseqs easy-linclust ${TMP}/sediment_seawater.onlycomp.faa ${PREFIX} ${TMP} \
    --min-seq-id 0.9 -c 0.8 \
    --threads ${SLURM_CPUS_PER_TASK} --split-memory-limit 500G

# Protein IDs of the catalogue (used to split it for DRAM in 2.3)
grep '^>' ${PREFIX}_rep_seq.fasta | sed 's/^>//; s/ .*//' > ${PREFIX}_rep_seq.ids.txt

cp ${PREFIX}_rep_seq.fasta ${PREFIX}_rep_seq.ids.txt ${PREFIX}_cluster.tsv ${OUTDIR}/
