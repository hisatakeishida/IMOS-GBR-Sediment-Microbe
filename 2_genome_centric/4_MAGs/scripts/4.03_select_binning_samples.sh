#!/bin/bash --login
#SBATCH --job-name=binning_selection
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o binning_selection_%j.o
#SBATCH -e binning_selection_%j.e

# MAG recovery, step 3: for each assembly, select the 20 Illumina datasets
# with the highest percentage of bacterial marker gene OTUs recovered in the
# assembly (singlem appraise --imperfect, assembly vs reads of all samples).
#
# Output: one file per assembly listing its 20 binning datasets (sample names
# as in the SingleM OTU tables), used by 4.4.

set -euo pipefail

INDIR="<path/to/4_MAG_recovery/2_binning_selection>"   # output of 4.2
METAP="<path/to/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb>"
N_SAMPLES=20
SORT_COL=5   # appraise output column used for ranking: 5 = percent of OTUs assembled

mkdir -p ${INDIR}/{appraise,binning_samples}

conda activate singlem   # SingleM v1.18.0

for otu in ${INDIR}/assembly_otu/*_assembly.otu_table.csv; do
    assembly=$(basename ${otu} _assembly.otu_table.csv)
    echo "Assembly: ${assembly}"

    singlem appraise --metapackage ${METAP} \
        --imperfect \
        --assembly-otu-tables ${otu} \
        --metagenome-otu-tables ${INDIR}/reads_otu/* \
        > ${INDIR}/appraise/${assembly}.imperfect.txt

    grep "d__Bacteria" ${INDIR}/appraise/${assembly}.imperfect.txt | grep -v "total\|average" | \
        sort -k ${SORT_COL},${SORT_COL}nr | head -n ${N_SAMPLES} | cut -f 1 | sort \
        > ${INDIR}/binning_samples/${assembly}.binning_samples.txt
done