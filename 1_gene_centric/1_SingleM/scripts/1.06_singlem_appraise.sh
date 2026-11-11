#!/bin/bash --login
#SBATCH --job-name=singlem_appraise
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o singlem_appraise_%j.o
#SBATCH -e singlem_appraise_%j.e

# Marker gene (mOTU) representation of metagenomes in assemblies, MAGs, and
# GTDB R10-RS226 representative genomes, using SingleM appraise.
#   1. OTU tables for assemblies, MAGs, and GTDB genomes (singlem pipe)
#   2. Read mOTUs recovered in assemblies (singlem appraise, assemblies)
#   3. Read mOTUs recovered in MAGs and GTDB genomes (singlem appraise, genomes)
#
# Recovery is assessed at species level (exact match) and genus level
# (--imperfect, 89% sequence identity).
#
# Run once per habitat, e.g.
#   sbatch --export=ALL,HABITAT=sed scripts/1.6_singlem_appraise.sh
#   sbatch --export=ALL,HABITAT=sea scripts/1.6_singlem_appraise.sh

set -euo pipefail

HABITAT=${HABITAT:?Set HABITAT=sed or HABITAT=sea}

# ---- Inputs ------------------------------------------------------------------
METAP="<path/to/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb>"

METAGENOME_OTU="<path/to/singlem_output>/${HABITAT}/R226/otu"    # from 1.2_singlem_pipe_array.sh
ASSEMBLY_DIR="<path/to/assemblies>/${HABITAT}"                   # *.fasta (single-sample and co-assemblies)
MAG_DIR="<path/to/species_MAGs>/${HABITAT}"                      # *.fna (1,201 sediment or 876 seawater)
GTDB_OTU="<path/to/gtdb_r226_reps_otu_table.csv>"                # singlem pipe on GTDB R10-RS226 representatives

OUTDIR="<path/to/appraise_output>/${HABITAT}"
GENUS_ID=0.89                                                  # genus-level identity for --imperfect

mkdir -p ${OUTDIR}/{otu_tables,appraise,unaccounted}

conda activate singlem   # SingleM v1.18.0
singlem --version


# ==============================================================================
# 1. OTU TABLES FOR ASSEMBLIES AND MAGs
# ==============================================================================

# Assemblies: one OTU table per assembly
for fa in ${ASSEMBLY_DIR}/*.fasta; do
    name=$(basename ${fa} .fasta)
    singlem pipe --metapackage ${METAP} \
        --genome-fasta-files ${fa} \
        --threads ${SLURM_CPUS_PER_TASK} \
        --otu-table ${OUTDIR}/otu_tables/assembly_${name}.otu.csv
done

# MAGs: one OTU table for all species-representative MAGs
singlem pipe --metapackage ${METAP} \
    --genome-fasta-files ${MAG_DIR}/*.fna \
    --threads ${SLURM_CPUS_PER_TASK} \
    --otu-table ${OUTDIR}/otu_tables/mags.otu.csv

# GTDB R10-RS226 representative genomes (run once, shared by both habitats):
# singlem pipe --metapackage ${METAP} --genome-fasta-files <gtdb_r226_reps>/*.fna.gz \
#     --threads ${SLURM_CPUS_PER_TASK} --otu-table ${GTDB_OTU}


# ==============================================================================
# 2. RECOVERY IN ASSEMBLIES
# ==============================================================================

# Species level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --assembly-otu-tables ${OUTDIR}/otu_tables/assembly_*.otu.csv \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_assemblies.species.tsv

# Genus level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --assembly-otu-tables ${OUTDIR}/otu_tables/assembly_*.otu.csv \
    --imperfect --sequence-identity ${GENUS_ID} \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_assemblies.genus.tsv


# ==============================================================================
# 3. RECOVERY IN MAGs AND GTDB GENOMES
# ==============================================================================

# MAGs, species level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --genome-otu-tables ${OUTDIR}/otu_tables/mags.otu.csv \
    --output-unaccounted-for-otu-table ${OUTDIR}/unaccounted/${HABITAT}_mags_unaccounted.species.csv \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_mags.species.tsv

# MAGs, genus level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --genome-otu-tables ${OUTDIR}/otu_tables/mags.otu.csv \
    --imperfect --sequence-identity ${GENUS_ID} \
    --output-unaccounted-for-otu-table ${OUTDIR}/unaccounted/${HABITAT}_mags_unaccounted.genus.csv \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_mags.genus.tsv

# GTDB R10-RS226, species level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --genome-otu-tables ${GTDB_OTU} \
    --output-unaccounted-for-otu-table ${OUTDIR}/unaccounted/${HABITAT}_gtdb_r226_unaccounted.species.csv \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_gtdb_r226.species.tsv

# GTDB R10-RS226, genus level
singlem appraise --metapackage ${METAP} \
    --metagenome-otu-tables ${METAGENOME_OTU}/* \
    --genome-otu-tables ${GTDB_OTU} \
    --imperfect --sequence-identity ${GENUS_ID} \
    --output-unaccounted-for-otu-table ${OUTDIR}/unaccounted/${HABITAT}_gtdb_r226_unaccounted.genus.csv \
    > ${OUTDIR}/appraise/${HABITAT}_reads_vs_gtdb_r226.genus.tsv
    