#!/bin/bash --login
#SBATCH --job-name=aviary_coassembly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=120:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o aviary_coassembly_%A_%a.o
#SBATCH -e aviary_coassembly_%A_%a.e
#SBATCH --array=1-71

# MAG recovery, step 6: co-assembly (MEGAHIT v1.2.9 within Aviary assemble)
# and co-binning (Aviary recover) of the 71 sample pairs chosen in 4.5, one
# array task per co-assembly. Each line of the command files from BinChicken
# is one Aviary command (output in the working directory as coassembly_<i>).
#
# Run assembly first, then binning:
#   sbatch --export=ALL,MODE=assemble scripts/4.6_aviary_coassembly_array.sh
#   sbatch --export=ALL,MODE=recover  scripts/4.6_aviary_coassembly_array.sh

set -euo pipefail

MODE=${MODE:?Set MODE=assemble or MODE=recover}
BINCHICKEN_DIR="<path/to/4_MAG_recovery/4_coassembly/binchicken>"   # output of 4.5

case ${MODE} in
    assemble) COMMANDS="${BINCHICKEN_DIR}/coassemble/commands/coassemble_commands.sh"
              OUTDIR="<path/to/4_MAG_recovery/4_coassembly/assemblies>" ;;
    recover)  COMMANDS="${BINCHICKEN_DIR}/coassemble/commands/recover_commands.sh"
              OUTDIR="<path/to/4_MAG_recovery/4_coassembly/bins>" ;;
    *)        echo "Unknown MODE: ${MODE}" >&2; exit 1 ;;
esac

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}
cd ${TMP}

CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${COMMANDS})
index=$(( SLURM_ARRAY_TASK_ID - 1 ))
echo "Co-assembly ${index}: ${CMD}"

conda activate aviary   # Aviary v0.9.1
eval ${CMD}

# Archive the large intermediate data folder before copying (binning only)
if [ "${MODE}" = "recover" ]; then
    cd ${TMP}/coassembly_${index}
    tar -cf data.tar data && rm -r data
    pigz data.tar
fi

cp -r ${TMP}/coassembly_${index} ${OUTDIR}/
