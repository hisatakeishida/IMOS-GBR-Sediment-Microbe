#!/bin/bash --login
#SBATCH --job-name=sylph
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o sylph_%A_%a.o
#SBATCH -e sylph_%A_%a.e
#SBATCH --array=1-<number_of_samples>

# k-mer-based taxonomic profiling with sylph against GTDB R10-RS226
# (gtdb-r226-c200-dbv1.syldb), on reads subsampled to 10 M read pairs.
#   1. sylph profile, with estimated read counts
#   2. sylph-tax taxprof: assign GTDB R10-RS226 taxonomy (.sylphmpa)
#
# Run once per habitat, e.g.
#   sbatch --export=ALL,HABITAT=sed --array=1-<n_sed> scripts/3.1_sylph_profile_array.sh
#   sbatch --export=ALL,HABITAT=sea --array=1-<n_sea> scripts/3.1_sylph_profile_array.sh
#
# One-off setup:
#   wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb
#   sylph-tax download --download-to <path/to/sylph-tax_db>

set -euo pipefail

HABITAT=${HABITAT:?Set HABITAT=sed or HABITAT=sea}
INDIR="<path/to/subsampled_reads>/${HABITAT}"        # <sample>.R1.10Mseqs.fq.gz and .R2.10Mseqs.fq.gz
DB="<path/to/gtdb-r226-c200-dbv1.syldb>"
OUTDIR="<path/to/3_Sylph>/${HABITAT}"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}/{1_count,2_sylphmpa}

R1=$(ls ${INDIR}/*.R1.10Mseqs.fq.gz | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${R1} .R1.10Mseqs.fq.gz)
R2=${INDIR}/${name}.R2.10Mseqs.fq.gz
echo "Habitat: ${HABITAT}; sample: ${name}"

conda activate sylph_env   # sylph and sylph-tax
sylph --version

# 1. Profile
sylph profile -d ${DB} -1 ${R1} -2 ${R2} \
    -t ${SLURM_CPUS_PER_TASK} \
    --estimate-read-counts \
    > ${TMP}/${name}_counts.tsv

# 2. Taxonomy (writes <name>_counts.sylphmpa in the working directory)
cd ${TMP}
sylph-tax taxprof ${TMP}/${name}_counts.tsv -t GTDB_r226

cp ${TMP}/${name}_counts.tsv ${OUTDIR}/1_count/
cp ${TMP}/${name}*.sylphmpa ${OUTDIR}/2_sylphmpa/