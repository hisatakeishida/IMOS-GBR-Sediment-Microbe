#!/bin/bash --login
#SBATCH --job-name=singlem_pipe_r220
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o singlem_pipe_r220_%A_%a.o
#SBATCH -e singlem_pipe_r220_%A_%a.e
#SBATCH --array=1-191

# MAG recovery, step 2: SingleM v1.18.0 pipe (GTDB R09-RS220 metapackage) on
# each sample's reads and on its single-sample assembly, to select the
# datasets used for binning (4.3).
#
# SAMPLES: same file as 4.1 (<sample> <R1> <R2> [<ONT>])

set -euo pipefail

SAMPLES="<path/to/assembly_samples.tsv>"
ASSEMBLY_DIR="<path/to/4_MAG_recovery/1_assembly>"   # output of 4.1
METAP="<path/to/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb>"
OUTDIR="<path/to/4_MAG_recovery/2_binning_selection>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}/{reads_otu,reads_prof,assembly_otu,assembly_prof}

IFS=$'\t' read -r sample R1 R2 ONT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES})
echo "Sample: ${sample}"

conda activate singlem   # SingleM v1.18.0
singlem --version

# Reads
singlem pipe --metapackage ${METAP} \
    -1 ${R1} -2 ${R2} \
    --threads ${SLURM_CPUS_PER_TASK} \
    -p ${OUTDIR}/reads_prof/${sample}_prof.tsv \
    --otu-table ${OUTDIR}/reads_otu/${sample}_otu.tsv

# Assembly
cp ${ASSEMBLY_DIR}/${sample}/assembly/final_contigs.fasta ${TMP}/${sample}_final_contigs.fasta
singlem pipe --metapackage ${METAP} \
    --genome-fasta-files ${TMP}/${sample}_final_contigs.fasta \
    --threads ${SLURM_CPUS_PER_TASK} \
    -p ${OUTDIR}/assembly_prof/${sample}_assembly.prof.tsv \
    --otu-table ${OUTDIR}/assembly_otu/${sample}_assembly.otu_table.csv
