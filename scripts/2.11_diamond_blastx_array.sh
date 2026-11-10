#!/bin/bash --login
#SBATCH --job-name=diamond_blastx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o diamond_blastx_%A_%a.o
#SBATCH -e diamond_blastx_%A_%a.e
#SBATCH --array=1-<number_of_samples>

# Read mapping, step 3: subsample R1 reads (length >= 140 bp) to 10 M reads and
# align them with DIAMOND blastx v2.1.0 to one database:
#   DB=catalogue : non-redundant protein catalogue (2.9)
#   DB=singlem   : SingleM marker genes (2.10), used to normalise RPKM
# (E-value 1e-5, min-score 40, query cover 80%, identity 70%, max-hsps 1,
# max-target-seqs 1). Samples with <= 10 M reads are excluded. Subsampling
# uses a fixed seed, so both databases are searched with the same reads.
#
# Run once per habitat and database, e.g.
#   sbatch --export=ALL,HABITAT=sed,DB=catalogue --array=1-191 scripts/2.11_diamond_blastx_array.sh
#   sbatch --export=ALL,HABITAT=sea,DB=catalogue --array=1-182 scripts/2.11_diamond_blastx_array.sh
#   sbatch --export=ALL,HABITAT=sed,DB=singlem   --array=1-191 scripts/2.11_diamond_blastx_array.sh
#   sbatch --export=ALL,HABITAT=sea,DB=singlem   --array=1-182 scripts/2.11_diamond_blastx_array.sh

set -euo pipefail

HABITAT=${HABITAT:?Set HABITAT=sed or HABITAT=sea}
DB=${DB:?Set DB=catalogue or DB=singlem}

case ${DB} in
    catalogue) DMND="<path/to/3_read_mapping_RPKM/dmnd>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.dmnd"
               EXT=".tsv" ;;
    singlem)   DMND="<path/to/1_SingleM_Allmarkers.dmnd>"
               EXT=".singlem.tsv" ;;
    *)         echo "Unknown DB: ${DB}" >&2; exit 1 ;;
esac

LIST="<path/to/read_lists>/${HABITAT}_R1_140bp.txt"   # one R1 file (reads >= 140 bp) per line
OUTDIR="<path/to/3_read_mapping_RPKM/diamond>/${HABITAT}/${DB}"
THRESHOLD=10000000
SEED=100

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}

R1=$(sort ${LIST} | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${R1} .fastq.gz)
echo "Habitat: ${HABITAT}; database: ${DB}; sample: ${name}"

conda activate seqkit_2.5.1
nreads=$(seqkit stats --tabular --threads ${SLURM_CPUS_PER_TASK} ${R1} | awk 'NR == 2 {print $4}')
echo "Reads: ${nreads}"

if [ "${nreads}" -le "${THRESHOLD}" ]; then
    echo "${name}: ${nreads} reads (<= ${THRESHOLD}); excluded."
    exit 0
fi

module load seqtk/1.3-gcc-11.3.0
seqtk sample -s${SEED} ${R1} ${THRESHOLD} | pigz -p ${SLURM_CPUS_PER_TASK} > ${TMP}/${name}.R1.140bp.10Mseqs.fq.gz

module load diamond/2.1.0-gcc-11.3.0
diamond blastx --threads ${SLURM_CPUS_PER_TASK} \
    -d ${DMND} \
    -q ${TMP}/${name}.R1.140bp.10Mseqs.fq.gz \
    -o ${TMP}/${name}${EXT} \
    --evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1

cp ${TMP}/${name}${EXT} ${OUTDIR}/