#!/bin/bash --login
#SBATCH --job-name=bigslice
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o bigslice_%j.o
#SBATCH -e bigslice_%j.e

# MAG annotation, step 13: BiG-SLiCE features and GCF distances for the BGCs
# of one habitat (input of 5.14).
#   1. BiG-SLiCE v1.1.1 on the antiSMASH region GenBank files (5.09)
#   2. perform_l2norm_clustering.py (from the OMD scripts; Paoli et al., 2022):
#      GCF models from L2-normalised BGC features (BIRCH, THRESHOLD), and the
#      distance of each BGC to its nearest GCF model
#
# TODO: confirm how distances to BiG-FAM were obtained. As run here, the GCF
# models come from this habitat's own BGCs, so the distances are not to
# BiG-FAM unless the BiG-SLiCE run includes the BiG-FAM BGCs or BiG-SLiCE
# --query is used against the BiG-FAM database.
#
# Run once per habitat, e.g.
#   sbatch --export=ALL,DATASET=IMOS_SED 2_genome_centric/5_MAGsAnnot/scripts/5.13_bigslice_gcf_distance.sh
#   sbatch --export=ALL,DATASET=IMOS_SEA 2_genome_centric/5_MAGsAnnot/scripts/5.13_bigslice_gcf_distance.sh

set -euo pipefail

DATASET=${DATASET:?Set DATASET=IMOS_SED or IMOS_SEA}
INPUT_TEMPLATE="<path/to/bigslice_input>/${DATASET}"            # BiG-SLiCE input folder (datasets.tsv, taxonomy/, dataset_1/)
GBK_DIR="<path/to/5_MAG_annotation/7_BGCs/antismash>/${DATASET}" # output of 5.09 (one folder per genome)
BIGSLICE="<path/to/bigslice-1.1.1/myenv/bin/bigslice>"
L2NORM="<path/to/repo>/2_genome_centric/5_MAGsAnnot/scripts/perform_l2norm_clustering.py"
OUTDIR="<path/to/5_MAG_annotation/7_BGCs/bigslice>"
THRESHOLD=0.4   # TODO: confirm (your bigslice.sh used 0.1; the input sheets of 5.14 are named gcf_dist_0.4)

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}
cd ${TMP}

# ---- 1. BiG-SLiCE ------------------------------------------------------------------------
cp -r ${INPUT_TEMPLATE} ${TMP}/${DATASET}
mkdir -p ${TMP}/${DATASET}/dataset_1
cp -r ${GBK_DIR}/* ${TMP}/${DATASET}/dataset_1/

module load hmmer/3.4-gompi-2023a
${BIGSLICE} -i ${TMP}/${DATASET} -t ${SLURM_CPUS_PER_TASK} ${TMP}/${DATASET}_result
cp -r ${TMP}/${DATASET}_result ${OUTDIR}/

# ---- 2. GCF models and BGC distances -----------------------------------------------------
conda activate clustomatic   # Python with pandas, NumPy, scikit-learn

python ${L2NORM} ${OUTDIR}/${DATASET}_result ${THRESHOLD} \
    ${OUTDIR}/${DATASET}_${THRESHOLD}_output_gcf_features.tsv \
    ${OUTDIR}/${DATASET}_${THRESHOLD}_output_bgc_membership.tsv
