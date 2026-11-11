#!/bin/bash --login
#SBATCH --job-name=aviary_recover
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=<memory: 100G Illumina-only, 400G hybrid>
#SBATCH --time=<time: 24:00:00 Illumina-only, 96:00:00 hybrid>
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o aviary_recover_%A_%a.o
#SBATCH -e aviary_recover_%A_%a.e
#SBATCH --array=1-191

# MAG recovery, step 4: multi-sample binning of each single-sample assembly
# with Aviary recover v0.9.1 (--binning-only), using the 20 datasets selected
# in 4.3. Aviary runs MetaBAT1 v2.15, MetaBAT2 v2.15, SemiBin v2.1.0,
# VAMB v3.0.2, and Rosella v0.5.2, selects bins with DAS Tool v1.1.2, and
# assesses them with CheckM and CheckM2 (used in 4.7).
# For hybrid assemblies the sample's ONT reads are also used.
#
# SAMPLES: same file as 4.1 (<sample> <R1> <R2> [<ONT>])

set -euo pipefail

SAMPLES="<path/to/assembly_samples.tsv>"
ASSEMBLY_DIR="<path/to/4_MAG_recovery/1_assembly>"                         # output of 4.1
SELECTION_DIR="<path/to/4_MAG_recovery/2_binning_selection/binning_samples>"  # output of 4.3
READS_DIR="<path/to/adapter_removed_reads>"                                # <sample>__<library>_R{1,2}_allRunsCombined.fastq.gz
OUTDIR="<path/to/4_MAG_recovery/3_binning>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${TMP}/reads

IFS=$'\t' read -r sample R1 R2 ONT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES})
echo "Assembly: ${sample}"

cp ${ASSEMBLY_DIR}/${sample}/assembly/final_contigs.fasta ${TMP}/${sample}_final_contigs.fasta

# Reads of the 20 selected datasets
while read -r binsample; do
    cp ${READS_DIR}/${binsample}*_R1_allRunsCombined.fastq.gz ${READS_DIR}/${binsample}*_R2_allRunsCombined.fastq.gz ${TMP}/reads/
done < ${SELECTION_DIR}/${sample}.binning_samples.txt
ls ${TMP}/reads/*_R1_allRunsCombined.fastq.gz | wc -l

LONGREADS=""
if [ -n "${ONT:-}" ]; then
    LONGREADS="--longreads ${ONT} --long_read_type ont_hq"
fi

conda activate aviary   # Aviary v0.9.1

aviary recover --binning-only \
    --assembly ${TMP}/${sample}_final_contigs.fasta \
    -1 ${TMP}/reads/*_R1_allRunsCombined.fastq.gz \
    -2 ${TMP}/reads/*_R2_allRunsCombined.fastq.gz \
    ${LONGREADS} \
    --max-threads ${SLURM_CPUS_PER_TASK} --n-cores ${SLURM_CPUS_PER_TASK} \
    --max_memory $(( SLURM_MEM_PER_NODE / 1024 )) \
    --conda-frontend conda \
    --tmpdir ${TMP} \
    --output ${TMP}/${sample}

cp -r ${TMP}/${sample} ${OUTDIR}/
