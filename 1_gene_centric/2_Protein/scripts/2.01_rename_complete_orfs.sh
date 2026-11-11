#!/bin/bash --login
#SBATCH --job-name=complete_orfs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=16:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o complete_orfs_%j.o
#SBATCH -e complete_orfs_%j.e

# Gene catalogue, step 1: rename Prodigal proteins, concatenate per habitat, and
# keep complete ORFs only (Prodigal partial=00).
#
# Headers are prefixed with the habitat (and, for sediment, the sample), e.g.
#   >sediment_11049_1_NODE_1_length_13274_cov_7.949013_1 # 237 # 380 # -1 # ID=1_1;partial=00;...
#   >seawater_11049_1__NODE_1_length_601848_cov_323.609474_1 # ...
# (seawater Prodigal headers already contain the sample name)
#
#                 all proteins     complete ORFs
#   Sediment      1,853,838,210    22,803,724
#   Seawater        182,544,749    38,136,288

set -euo pipefail

SED_PRODIGAL="<path/to/prodigal/sediment>"   # <sample>.protein.translations.faa.gz
SEA_PRODIGAL="<path/to/prodigal/seawater>"   # <sample>.protein.translations.faa
OUTDIR="<path/to/1_gene_catalogue>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

# ---- Rename and concatenate ---------------------------------------------------
sed_all=${TMP}/imos_sediment_concatenated.protein.translations.renamed.faa
sea_all=${TMP}/imos_seawater_concatenated.protein.translations.renamed.faa
> ${sed_all}
> ${sea_all}

for f in ${SED_PRODIGAL}/*.protein.translations.faa.gz; do
    sample=$(basename ${f} .protein.translations.faa.gz)
    zcat ${f} | awk -v p="sediment_${sample}_" '/^>/ {print ">" p substr($0, 2); next} {print}' >> ${sed_all}
done

for f in ${SEA_PRODIGAL}/*.protein.translations.faa; do
    awk '/^>/ {print ">seawater_" substr($0, 2); next} {print}' ${f} >> ${sea_all}
done

# ---- Keep complete ORFs -------------------------------------------------------
conda activate seqkit_2.5.1   # SeqKit v2.5.1
# 00 = complete ORFs

for prots in ${sed_all} ${sea_all}; do
    ids=${prots%.faa}.onlycomp.ids.txt
    awk '/^>/ && /partial=00/ {sub(/^>/, "", $1); print $1}' ${prots} > ${ids}
    seqkit grep -j ${SLURM_CPUS_PER_TASK} -f ${ids} ${prots} > ${OUTDIR}/$(basename ${prots%.faa}).onlycomp.faa
done

seqkit stats -T ${OUTDIR}/*.onlycomp.faa
