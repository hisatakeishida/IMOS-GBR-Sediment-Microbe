#!/bin/bash --login
#SBATCH --job-name=marker_blastp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o marker_blastp_%A_%a.o
#SBATCH -e marker_blastp_%A_%a.e
#SBATCH --array=1-<number_of_marker_databases>

# Annotation: metabolic marker genes (Greening lab databases) with DIAMOND
# blastp v2.1.0, one array task per marker gene database. 
# Minimum query coverage 80%, with gene-specific identity thresholds (Accroding to Ricci et al., 2026):
#   80%: psaA
#   70%: isoA, ygfK, aro
#   65%: hbsT
#   60%: atpA, psbA, amoA, mmoA, coxL, [FeFe]-hydrogenase, nxrA, rbcL, nuoF
#   45%: rdhA;  35%: cyc2;  30%: rho
#   50%: all remaining genes
# Best hits across all databases are combined in 2.8.

set -euo pipefail

GENE_DB_DIR="<path/to/GreeningLab-database/Compiled>"           # one .fasta/.faa per gene
DM_DIR="<path/to/GreeningLab-database/Compiled_diamond>"
CATALOGUE="<path/to/1_gene_catalogue>/imos_sediment_seawater_ass_clustered.prot.componly_0.9id_0.8len_rep_seq.fasta"
OUTDIR="<path/to/2_annotation/marker_genes>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${DM_DIR}

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

# This task's database (bash arrays start at 0, SLURM array at 1)
mapfile -t FASTA_LIST < <(find ${GENE_DB_DIR} -maxdepth 1 \( -name "*.fasta" -o -name "*.faa" \) | sort)
FASTA=${FASTA_LIST[$(( SLURM_ARRAY_TASK_ID - 1 ))]}
GENE=$(basename "${FASTA}")
GENE=${GENE%.fasta}; GENE=${GENE%.faa}; GENE=$(echo "${GENE}" | tr ' ' '_')
PIDENT_VAL=${PIDENT[${GENE}]:-${DEFAULT_PIDENT}}
DM_DB=${DM_DIR}/${GENE}.dmnd
echo "Task ${SLURM_ARRAY_TASK_ID} of ${#FASTA_LIST[@]}: ${GENE}, identity >= ${PIDENT_VAL}%"

module load diamond/2.1.0-gcc-11.3.0

if [[ ! -f ${DM_DB} ]]; then
    diamond makedb --in "${FASTA}" --db ${DM_DB} --threads ${SLURM_CPUS_PER_TASK}
fi

diamond blastp \
    --query ${CATALOGUE} \
    --db ${DM_DB} \
    --out ${TMP}/${GENE}.blastp.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp \
    --id ${PIDENT_VAL} \
    --query-cover 80 \
    --max-target-seqs 1 \
    --threads ${SLURM_CPUS_PER_TASK}

cp ${TMP}/${GENE}.blastp.tsv ${OUTDIR}/
