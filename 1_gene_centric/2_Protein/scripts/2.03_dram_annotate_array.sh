#!/bin/bash --login
#SBATCH --job-name=dram_annotate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o dram_annotate_%A_%a.o
#SBATCH -e dram_annotate_%A_%a.e
#SBATCH --array=1-393

# Annotation: KOfam with DRAM v1.5 annotate_genes. The catalogue is annotated in
# chunks of 100,000 proteins (39,225,676 proteins = 393 array tasks).
# DRAM prefixes gene IDs with the input file name (<start>_<end>_); this prefix
# is removed in 2.8.

set -euo pipefail

CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
IDS=${CATALOGUE%.fasta}.ids.txt
OUTDIR="<path/to/2_annotation/dram>"
DRAM_ENV="<path/to/DRAM_conda_env>"

TOTAL=$(wc -l < ${IDS})
CHUNK=100000
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))
[ ${END} -gt ${TOTAL} ] && END=${TOTAL}
echo "Proteins ${START} to ${END}"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

sed -n "${START},${END}p" ${IDS} > ${TMP}/ids_${START}_${END}.txt

conda activate seqkit_2.5.1
seqkit grep -j ${SLURM_CPUS_PER_TASK} -f ${TMP}/ids_${START}_${END}.txt ${CATALOGUE} -o ${TMP}/${START}_${END}.faa

conda activate ${DRAM_ENV}   # DRAM v1.5
DRAM.py annotate_genes -i ${TMP}/${START}_${END}.faa -o ${TMP}/${START}_${END}__annotation \
    --threads ${SLURM_CPUS_PER_TASK}

cp -r ${TMP}/${START}_${END}__annotation ${OUTDIR}/
