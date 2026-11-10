#!/bin/bash --login
#SBATCH --job-name=gtdbtk_taxonomy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o gtdbtk_taxonomy_%j.o
#SBATCH -e gtdbtk_taxonomy_%j.e

# MAG recovery, step 9: taxonomy and phylogeny of the 1,201 sediment and 876
# seawater MAGs (GTDB R10-RS226).
#   1. GTDB-Tk v2.4.1 classify_wf
#   2. GTDB-Tk de_novo_wf trees (bacteria and archaea)
#   3. PhyloRank v0.1.12: decorate the de novo trees and calculate relative
#      evolutionary distance (RED) outliers per rank
#   4. Manual curation of the classify_wf taxonomy using the decorated trees
#      and RED values (not scripted)
#   5. IQ-TREE v2.2.2.3 trees from the bac120 and arc53 marker gene alignments
#      (1,000 ultrafast bootstraps)

set -euo pipefail

SED_MAGS="<path/to/4_MAG_recovery/6_dereplication/representatives>"   # output of 4.8
SEA_MAGS="<path/to/seawater_MAGs>"                                      # 876 MGD-seawater MAGs (Robbins et al., 2026)
OUTDIR="<path/to/4_MAG_recovery/7_taxonomy>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR} ${TMP}/genomes

cp ${SED_MAGS}/*.fna ${SEA_MAGS}/*.fna ${TMP}/genomes/
ls ${TMP}/genomes | wc -l

conda activate gtdbtk-2.4.1   # GTDB-Tk v2.4.1, GTDB R10-RS226 reference data


# ==============================================================================
# 1. CLASSIFICATION
# ==============================================================================

gtdbtk classify_wf \
    --genome_dir ${TMP}/genomes -x fna \
    --skip_ani_screen \
    --cpus ${SLURM_CPUS_PER_TASK} --pplacer_cpus ${SLURM_CPUS_PER_TASK} \
    --out_dir ${OUTDIR}/classify_wf \
    --tmpdir ${TMP}


# ==============================================================================
# 2. DE NOVO TREES
# ==============================================================================
# Outgroups: p__Chloroflexota (bacteria), p__Undinarchaeota (archaea)

gtdbtk de_novo_wf --bacteria \
    --genome_dir ${TMP}/genomes -x fna \
    --outgroup_taxon p__Chloroflexota \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --out_dir ${OUTDIR}/de_novo_wf_bac \
    --tmpdir ${TMP}

gtdbtk de_novo_wf --archaea \
    --genome_dir ${TMP}/genomes -x fna \
    --outgroup_taxon p__Undinarchaeota \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --out_dir ${OUTDIR}/de_novo_wf_ar \
    --tmpdir ${TMP}


# ==============================================================================
# 3. PHYLORANK
# ==============================================================================
# tax_file: GTDB R10-RS226 reference taxonomy plus classify_wf taxonomy of the
# MAGs (genome ID <tab> taxonomy string)

conda activate phylorank   # PhyloRank v0.1.12

for domain in bac ar; do
    marker=$([ ${domain} = bac ] && echo bac120 || echo ar53)
    intree=${OUTDIR}/de_novo_wf_${domain}/gtdbtk.${marker}.decorated.tree
    tax_file="<path/to/${marker}_taxonomy.tsv>"
    out_tree=${OUTDIR}/phylorank_${domain}/${marker}.phylorank_decorated.tree
    mkdir -p ${OUTDIR}/phylorank_${domain}

    phylorank decorate ${intree} ${tax_file} ${out_tree} --skip_rd_refine
    phylorank outliers ${out_tree} ${tax_file} ${OUTDIR}/phylorank_${domain}
done


# ==============================================================================
# 5. PHYLOGENY
# ==============================================================================
# Input: bac120 and ar53 user MSAs from classify_wf (step 1)

conda activate iqtree   # IQ-TREE v2.2.2.3

zcat ${OUTDIR}/classify_wf/align/gtdbtk.bac120.user_msa.fasta.gz > ${TMP}/bac120.msa.fasta
iqtree -s ${TMP}/bac120.msa.fasta -m LG+G -B 1000 -T ${SLURM_CPUS_PER_TASK} \
    --prefix ${OUTDIR}/iqtree_bac120

zcat ${OUTDIR}/classify_wf/align/gtdbtk.ar53.user_msa.fasta.gz > ${TMP}/ar53.msa.fasta
iqtree -s ${TMP}/ar53.msa.fasta -m LG+C10+F+G -B 1000 -T ${SLURM_CPUS_PER_TASK} \
    --prefix ${OUTDIR}/iqtree_ar53