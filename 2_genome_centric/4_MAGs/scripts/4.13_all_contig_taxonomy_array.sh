#!/bin/bash
#SBATCH --job-name=mmseqs_contigs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=mmseqs_contigs_%A_%a.o
#SBATCH --error=mmseqs_contigs_%A_%a.e
# =============================================================================
# 4.13_all_contig_taxonomy_array.sh
#
# Taxonomy of ALL assembled contigs (binned and unbinned) with MMseqs2
# easy-taxonomy, one SLURM array task per assembly (single-sample assemblies
# and co-assemblies).
#
# Database: custom MMseqs2 taxonomy database of RefSeq, EukProt, PLSDB, IMG/VR,
# and proteins predicted from GTDB R08-RS214 genomes (Robbins et al., 2026).
#
# Assembly list (tab or space separated, one assembly per line, no header):
#   <assembly_id>  <path/to/contigs.fasta>
#
# Usage:
#   n=$(wc -l < assembly_list.tsv)
#   sbatch --array=1-${n} 2_genome_centric/4_MAGs/scripts/4.13_all_contig_taxonomy_array.sh
#
# Output, per assembly, in ${OUT_DIR}/<assembly_id>/:
#   <assembly_id>_lca.tsv            MMseqs2 LCA assignment per contig
#   <assembly_id>_report             Kraken-style taxonomy report
#   <assembly_id>_tophit_report      top-hit report
#   <assembly_id>_domain_summary.tsv contigs and bp per domain
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration (edit, or override with sbatch --export=ALL,VAR=value)
# ---------------------------------------------------------------------------
ASSEMBLY_LIST=${ASSEMBLY_LIST:-assembly_list.tsv}
MMSEQS_DB=${MMSEQS_DB:-/path/to/mmseqs_custom_taxonomy_db/seqTaxDB}
OUT_DIR=${OUT_DIR:-2_genome_centric/4_MAGs/5_contig_taxonomy}
MMSEQS=${MMSEQS:-mmseqs}                     # MMseqs2 v15.6f452
THREADS=${SLURM_CPUS_PER_TASK:-32}
SPLIT_MEM=${SPLIT_MEM:-200G}                 # keep below --mem

# Activate the MMseqs2 environment here if needed, e.g.
# conda activate mmseqs2_v15

# ---------------------------------------------------------------------------
# Select this task's assembly
# ---------------------------------------------------------------------------
TASK_ID=${SLURM_ARRAY_TASK_ID:?Run as a SLURM array job}
read -r ASSEMBLY_ID CONTIGS < <(awk -v n="${TASK_ID}" 'NR == n {print $1, $2}' "${ASSEMBLY_LIST}")

if [[ -z "${ASSEMBLY_ID:-}" || ! -s "${CONTIGS:-}" ]]; then
    echo "ERROR: no assembly or missing contigs file for task ${TASK_ID} in ${ASSEMBLY_LIST}" >&2
    exit 1
fi

SAMPLE_OUT="${OUT_DIR}/${ASSEMBLY_ID}"
PREFIX="${SAMPLE_OUT}/${ASSEMBLY_ID}"
TMP="${TMPDIR:-/tmp}/mmseqs_${ASSEMBLY_ID}_${SLURM_JOB_ID:-$$}"
mkdir -p "${SAMPLE_OUT}" "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

echo "Task ${TASK_ID}: ${ASSEMBLY_ID}"
echo "Contigs: ${CONTIGS}"
echo "Database: ${MMSEQS_DB}"
date

# ---------------------------------------------------------------------------
# 1. MMseqs2 easy-taxonomy on all contigs
# ---------------------------------------------------------------------------
"${MMSEQS}" easy-taxonomy \
        "${CONTIGS}" \
        "${MMSEQS_DB}" \
        "${PREFIX}" \
        "${TMP}" \
        --tax-lineage 1 \
        --threads "${THREADS}" \
        --split-memory-limit "${SPLIT_MEM}"
