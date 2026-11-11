#!/bin/bash --login
#SBATCH --job-name=antismash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=16:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o antismash_%A_%a.o
#SBATCH -e antismash_%A_%a.e

# MAG annotation, step 9: biosynthetic gene clusters (BGCs) with antiSMASH
# v6.1.1 (Prodigal gene finding, minimum contig length 5 kbp), 50 genomes per
# array task.
#
# Run once per dataset, e.g.
#   sbatch --export=ALL,DATASET=IMOS_SED --array=1-25  2_genome_centric/5_MAGsAnnot/scripts/5.09_antismash_array.sh   # 1,201 MAGs
#   sbatch --export=ALL,DATASET=IMOS_SEA --array=1-18  2_genome_centric/5_MAGsAnnot/scripts/5.09_antismash_array.sh   # 876 MAGs
#   sbatch --export=ALL,DATASET=OMD_DB   --array=1-697 2_genome_centric/5_MAGsAnnot/scripts/5.09_antismash_array.sh   # 34,815 OMD genomes
# (OMD BGCs are restricted to the 8,422 dereplicated OMD genomes in 5.10.)

set -euo pipefail

DATASET=${DATASET:?Set DATASET=IMOS_SED, IMOS_SEA, or OMD_DB}
case ${DATASET} in
    IMOS_SED) INDIR="<path/to/sediment_MAGs>";  EXT=.fna ;;
    IMOS_SEA) INDIR="<path/to/seawater_MAGs>";  EXT=.fna ;;
    OMD_DB)   INDIR="<path/to/OMD_genomes>";    EXT=.fa  ;;
    *)        echo "Unknown DATASET: ${DATASET}" >&2; exit 1 ;;
esac
OUTDIR="<path/to/5_MAG_annotation/7_BGCs/antismash>/${DATASET}"   # region GenBank files per genome
OUTTAR="<path/to/5_MAG_annotation/7_BGCs/antismash_tar>/${DATASET}" # full antiSMASH output
ANTISMASH_DB="<path/to/antismash/database_v6.1.1>"
MINLEN=5000
PER_TASK=50

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${OUTTAR}
cd ${TMP}

START=$(( (SLURM_ARRAY_TASK_ID - 1) * PER_TASK + 1 ))
END=$(( START + PER_TASK - 1 ))

conda activate antismash   # antiSMASH v6.1.1

find ${INDIR} -maxdepth 1 -type f -name "*${EXT}" | sort | sed -n "${START},${END}p" | while read -r genome; do
    name=$(basename ${genome} ${EXT})
    out=${name}_${MINLEN}
    echo "Genome: ${name}"

    antismash ${genome} \
        -c ${SLURM_CPUS_PER_TASK} \
        --databases ${ANTISMASH_DB} \
        --genefinding-tool prodigal \
        --minlength ${MINLEN} \
        --output-dir ${TMP}/${out} \
        --output-basename ${name}_minlen_${MINLEN}

    mkdir -p ${OUTDIR}/${out}
    cp ${TMP}/${out}/*.gbk ${OUTDIR}/${out}/ 2>/dev/null || true
    tar -czf ${OUTTAR}/${out}.tar.gz -C ${TMP} ${out} && rm -r ${TMP}/${out}
done
