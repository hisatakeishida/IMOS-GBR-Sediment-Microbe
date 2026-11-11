#!/bin/bash --login
#SBATCH --job-name=singlem_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o singlem_%A_%a.o
#SBATCH -e singlem_%A_%a.e
#SBATCH --array=1-<number_of_samples>

# Step 2 of the mOTU workflow: SingleM pipe and microbial_fraction on the
# 10 M read-pair subsampled reads (output of subsample_10M_array.sh).

HABITAT=sed #HABITAT=sea
INDIR=<path/to/subsampled_reads>/${HABITAT}
OUTDIR=<path/to/singlem_output>/${HABITAT}/R226
METAP=<path/to/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb>

mkdir -p ${OUTDIR}/{prof,otu,mfrac}

R1=$(ls ${INDIR}/*_R1_allRunsCombined.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${R1} _R1_allRunsCombined.fastq.gz)
R2=${INDIR}/${name}_R2_allRunsCombined.fastq.gz
echo "Habitat: ${HABITAT}; sample: ${name}"

conda activate singlem   # SingleM v1.18.0 (environment from BinChicken v0.12.11)

singlem pipe --metapackage ${METAP} \
    -1 ${R1} -2 ${R2} \
    --threads ${SLURM_CPUS_PER_TASK} \
    -p ${OUTDIR}/prof/${name}_prof.tsv \
    --otu-table ${OUTDIR}/otu/${name}_otu.tsv

singlem microbial_fraction --metapackage ${METAP} \
    -1 ${R1} -2 ${R2} \
    -p ${OUTDIR}/prof/${name}_prof.tsv \
    --output-tsv ${OUTDIR}/mfrac/${name}_microbial_fraction.tsv
