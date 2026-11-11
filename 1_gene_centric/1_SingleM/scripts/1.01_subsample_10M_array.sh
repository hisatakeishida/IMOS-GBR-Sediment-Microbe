#!/bin/bash --login
#SBATCH --job-name=subsample_10M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o subsample_%A_%a.o
#SBATCH -e subsample_%A_%a.e
#SBATCH --array=1-<number_of_samples>

# Step 1 of the mOTU workflow: subsample adapter-removed reads to 10 M read pairs.
# Samples with <= 10 M read pairs are not subsampled and are excluded from
# downstream SingleM profiling.
#
# LIST     : text file with one R1 path per line (<sample>_R1_allRunsCombined.fastq.gz)
# READSOUT : output directory for subsampled reads (one per habitat)

set -euo pipefail

LIST=<path/to/R1_list.txt>
READSOUT=<path/to/subsampled_reads/{sed,sea}>
THRESHOLD=10000000
SEED=100

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${TMP} ${READSOUT}

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST})
name=$(basename ${R1} _R1_allRunsCombined.fastq.gz)
R2=$(dirname ${R1})/${name}_R2_allRunsCombined.fastq.gz
echo "Sample: ${name}"

# Count read pairs
conda activate seqkit_2.5.1
nreads=$(seqkit stats --tabular --threads ${SLURM_CPUS_PER_TASK} ${R1} | awk 'NR==2 {print $4}')
echo "Read pairs: ${nreads}"

if [ "${nreads}" -le "${THRESHOLD}" ]; then
    echo "${name}: ${nreads} read pairs (<= ${THRESHOLD}); not subsampled, excluded."
    exit 0
fi

# Subsample R1 and R2 with the same seed so read pairs stay matched
module load seqtk
seqtk sample -s${SEED} ${R1} ${THRESHOLD} | pigz -p ${SLURM_CPUS_PER_TASK} > ${TMP}/${name}_R1_allRunsCombined.fastq.gz
seqtk sample -s${SEED} ${R2} ${THRESHOLD} | pigz -p ${SLURM_CPUS_PER_TASK} > ${TMP}/${name}_R2_allRunsCombined.fastq.gz
cp ${TMP}/${name}_R{1,2}_allRunsCombined.fastq.gz ${READSOUT}/
