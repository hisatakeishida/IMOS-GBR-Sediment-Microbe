#!/bin/bash --login
#SBATCH --job-name=hyddb_markers_mags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o hyddb_markers_mags_%A_%a.o
#SBATCH -e hyddb_markers_mags_%A_%a.e
#SBATCH --array=1-2077

# MAG annotation, step 4: hydrogenases (HydDB) and 44 metabolic marker genes
# with DIAMOND blastp v2.1.0, on the proteins predicted by DRAM.
#   HydDB: best hit, class-specific identity thresholds
#     [NiFe] Group 4 > 50%; other [NiFe] > 30%; [FeFe] > 45%; [Fe] > 50%
#   Marker genes: one search per marker gene database, gene-specific identity
#     thresholds (as in 2.06); hits kept if alignment length >= 40 aa or query
#     or subject coverage >= 80%; best hit per protein across databases by
#     E-value

set -euo pipefail

DRAM_DIR="<path/to/5_MAG_annotation/1_DRAM>"                       # output of 5.01
HYDDB="<path/to/HydDB/hyddb.dmnd>"
GENE_DB_DIR="<path/to/GreeningLab-database/Compiled>"             # one .fasta/.faa per gene
DM_DIR="<path/to/GreeningLab-database/Compiled_diamond>"          # DIAMOND databases (built in 2.06)
OUTDIR="<path/to/5_MAG_annotation/4_hydrogenases_marker_genes>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}/{hyddb_raw,hyddb,marker_raw,marker_best_hits}

dir=$(ls -d ${DRAM_DIR}/*__annotation | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(basename ${dir} __annotation)
faa=${dir}/genes.faa
echo "MAG: ${name}"

module load diamond/2.1.0-gcc-11.3.0

# ---- HydDB --------------------------------------------------------------------------------
diamond blastp --threads ${SLURM_CPUS_PER_TASK} \
    -q ${faa} -d ${HYDDB} \
    --max-target-seqs 1 --header simple \
    --outfmt 6 qseqid sseqid pident evalue bitscore full_sseq \
    -o ${OUTDIR}/hyddb_raw/${name}.tsv

awk 'BEGIN {IGNORECASE = 1; print "query_id\tclosest_hydrogenase_group"}
NR > 1 {
    split($2, arr, /\|/); class = arr[3]; id = $3 + 0
    if (class ~ /\[NiFe\]_Group_4([a-i])?$/) {if (id > 50) print $1 "\t" class; next}
    if      (class ~ /\[NiFe\]/ && id > 30) print $1 "\t" class
    else if (class ~ /\[FeFe\]/ && id > 45) print $1 "\t" class
    else if (class ~ /\[Fe\]/   && id > 50) print $1 "\t" class
}' ${OUTDIR}/hyddb_raw/${name}.tsv > ${OUTDIR}/hyddb/${name}.tsv

# ---- Metabolic marker genes -------------------------------------------------------------
declare -A PIDENT
PIDENT["Photosystem_I_reaction_centre_PsaA_sequences"]=80
PIDENT["Isoprene_monooxygenase_subunit_A"]=70
PIDENT["Selenate_reductase_YgfK"]=70
PIDENT["Arsenite_oxidase_ARO"]=70
PIDENT["Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences"]=65
PIDENT["F-type_ATP_synthase_AtpA"]=60
PIDENT["Photosystem_II_reaction_centre_PsbA_sequences"]=60
PIDENT["Ammonia_monooxygenase_AmoA_Feb2024"]=60
PIDENT["Soluble_methane_monooxygenase_MmoA"]=60
PIDENT["Carbon_monoxide_dehydrogenase_CoxL"]=60
PIDENT["FeFe_hydrogenase"]=60
PIDENT["Nitrite_oxidoreductase_NxrA"]=60
PIDENT["Rubisco_RbcL_95_Filtered"]=60
PIDENT["NADH-ubiquinone_oxidoreductase_NuoF"]=60
PIDENT["Reductive_dehalogenase_RdhA"]=45
PIDENT["Iron_oxidising_cytochrome_Cyc2"]=35
PIDENT["Microbial_rhodopsin_RHO"]=30
DEFAULT_PIDENT=50

mkdir -p ${TMP}/${name}
mapfile -t FASTA_LIST < <(find ${GENE_DB_DIR} -maxdepth 1 \( -name "*.fasta" -o -name "*.faa" \) | sort)
for FASTA in "${FASTA_LIST[@]}"; do
    GENE=$(basename "${FASTA}"); GENE=${GENE%.fasta}; GENE=${GENE%.faa}; GENE=$(echo "${GENE}" | tr ' ' '_')
    diamond blastp --threads ${SLURM_CPUS_PER_TASK} \
        --query ${faa} --db ${DM_DIR}/${GENE}.dmnd \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp \
        --id ${PIDENT[${GENE}]:-${DEFAULT_PIDENT}} \
        --max-target-seqs 1 \
        --out ${TMP}/${name}/${GENE}.blastp.tsv
done
tar -czf ${OUTDIR}/marker_raw/${name}.tar.gz -C ${TMP} ${name}

# Filter hits and keep the best hit per protein (lowest E-value, column 7)
for tsv in ${TMP}/${name}/*.blastp.tsv; do
    gene=$(basename ${tsv} .blastp.tsv)
    awk -v gene="${gene}" 'BEGIN {OFS = "\t"} ($4 >= 40 || $9 >= 80 || $10 >= 80) {print $0, gene}' ${tsv}
done | awk 'BEGIN {OFS = "\t"} {
    if (!($1 in best) || $7 < best[$1]) {best[$1] = $7; line[$1] = $0}
} END {for (q in line) print line[q]}' > ${OUTDIR}/marker_best_hits/${name}.tsv
