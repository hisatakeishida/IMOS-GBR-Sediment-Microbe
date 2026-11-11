#!/bin/bash --login
#SBATCH --job-name=gcf_summary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o gcf_summary_%j.o
#SBATCH -e gcf_summary_%j.e

# MAG annotation, step 12: summarise clust-o-matic GCFs (5.11).
#   1. One row per BGC: dataset, BGC, contig, length, antiSMASH region
#      product, and GCF (MIBiG BGCs labelled MIBIG)
#   2. Number of BGCs from each dataset in each GCF (GCF x dataset)
#   3. Number of BGCs and GCFs per dataset

set -euo pipefail

BGC_DIR="<path/to/5_MAG_annotation/7_BGCs/bgc_proteins>"   # output of 5.10 (*_summary_info.csv)
RMD_INFO="<path/to/RMDv1_derep-antismash_summary_info.csv>"
NAME=clustomatic__imos_sea_sed__derep_omd__derep_rmd__mibig_5000
CLUST_DIR="<path/to/5_MAG_annotation/7_BGCs>/${NAME}"      # output of 5.11
CLUSTERS=${CLUST_DIR}/cluster_assignments.txt

summary=${CLUST_DIR}/${NAME}_BGCS_GCF_withmibig.tsv
matrix=${CLUST_DIR}/${NAME}_BGCS_GCF_withmibig.countmatrix.tsv

# ---- 1. BGC table ---------------------------------------------------------------------------
# Metadata keyed by the BGC ID used in the protein headers ("_" replaced by "-")
cat ${BGC_DIR}/IMOS_SED_summary_info.csv ${BGC_DIR}/IMOS_SEA_summary_info.csv \
    ${BGC_DIR}/OMD_DEREP_summary_info.csv ${RMD_INFO} > ${CLUST_DIR}/all_summary_info.csv

awk -F'\t' -v OFS='\t' '
    BEGIN {print "SOURCE", "BGC_NAME", "BGC_NAME_MODIFIED", "CONTIG_NAME", "BGC_length", "BGC_FUNCTION", "GCF_ID"}
    NR == FNR {split($0, a, ","); id = a[2]; gsub("_", "-", id); info[id] = a[1] OFS a[2] OFS id OFS a[4] OFS a[5] OFS a[6]; next}
    FNR == 1 && $2 !~ /^[0-9]+$/ {next}   # header line, if present
    $1 ~ /^MIBIG/ {print "MIBIG", $1, $1, $1, "NA", "NA", $2; next}
    ($1 in info) {print info[$1], $2; next}
    {print "", "", $1, "", "", "", $2}
' ${CLUST_DIR}/all_summary_info.csv ${CLUSTERS} > ${summary}

# ---- 2. GCF x dataset count matrix -----------------------------------------------------------
awk -F'\t' -v OFS='\t' '
    NR > 1 {n[$7, $1]++; gcf[$7]; src[$1]}
    END {
        printf "GCF_ID"; for (s in src) printf "\t%s", s; print ""
        for (g in gcf) {printf "%s", g; for (s in src) printf "\t%d", n[g, s]; print ""}
    }' ${summary} > ${matrix}

# ---- 3. BGCs and GCFs per dataset ------------------------------------------------------------
echo -e "Dataset\tBGCs\tGCFs"
tail -n +2 ${summary} | awk -F'\t' '{bgc[$1]++; if (!(($1, $7) in seen)) {seen[$1, $7]; gcf[$1]++}}
    END {for (s in bgc) print s "\t" bgc[s] "\t" gcf[s]}'
