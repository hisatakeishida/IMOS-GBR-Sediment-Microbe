#!/bin/bash --login
#SBATCH --job-name=binchicken_coassemble
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o binchicken_coassemble_%j.o
#SBATCH -e binchicken_coassemble_%j.e

# MAG recovery, step 5: choose pairs of samples for co-assembly, and the 20
# datasets for co-binning each pair, with BinChicken v0.12.11 coassemble,
# based on SingleM marker gene profiles (GTDB R09-RS220 metapackage).
# BinChicken writes Aviary assemble and recover commands, which are run in 4.6.

set -euo pipefail

READS_DIR="<path/to/adapter_removed_reads>"
METAP="<path/to/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb>"
OUTDIR="<path/to/4_MAG_recovery/4_coassembly/binchicken>"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

conda activate binchicken   # BinChicken v0.12.11

# TODO: confirm --max-coassembly-size and any other options used
binchicken coassemble \
    --cores ${SLURM_CPUS_PER_TASK} \
    --singlem-metapackage ${METAP} \
    --forward ${READS_DIR}/*_R1_allRunsCombined.fastq.gz \
    --reverse ${READS_DIR}/*_R2_allRunsCombined.fastq.gz \
    --num-coassembly-samples 2 \
    --max-coassembly-samples 2 \
    --max-recovery-samples 20 \
    --max-coassembly-size <max_coassembly_size_Gbp> \
    --output ${OUTDIR}
