#!/bin/bash --login
#SBATCH --job-name=bgc_proteins
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o bgc_proteins_%j.o
#SBATCH -e bgc_proteins_%j.e

# MAG annotation, step 10: prepare BGC protein sets for clust-o-matic (5.11).
# For each dataset:
#   1. Copy antiSMASH region GenBank files, renamed <genome>__BGC<n>.gbk
#   2. Record BGC metadata (dataset, BGC, contig, length, antiSMASH products)
#   3. Extract CDS translations, one FASTA per BGC (headers
#      <BGC ID with "_" replaced by "-">_protein<n>), and keep BGCs >= 5 kbp
#   4. Concatenate per dataset
# OMD BGCs are restricted to the 8,422 OMD genomes dereplicated at 95% ANI.
# MIBiG v4.0 protein sequences are renamed to MIBIG<BGC>_<n>.

set -euo pipefail

ANTISMASH_DIR="<path/to/5_MAG_annotation/7_BGCs/antismash>"   # output of 5.09 (one folder per dataset)
OMD_REPS="<path/to/OMD_ANI95_AF15.rep.list.txt>"               # dereplicated OMD genomes (paths)
MIBIG_FAA="<path/to/mibig_prot_seqs_4.0.fasta>"
OUTDIR="<path/to/5_MAG_annotation/7_BGCs/bgc_proteins>"
MINLEN=5000

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}
LOG=${OUTDIR}/bgc_proteins.log
> ${LOG}

# CDS translations of one GenBank file as FASTA
extract_cds() {
    local gbk=$1 clusterid counter=1 in_cds=0 seq="" line
    clusterid=$(basename ${gbk} .gbk | tr '_' '-')
    while IFS= read -r line; do
        [[ ${line} =~ ^[[:space:]]*CDS ]] && { in_cds=1; seq=""; }
        if [ ${in_cds} -eq 1 ]; then
            if [[ ${line} =~ /translation=\"(.*) ]]; then
                seq="${BASH_REMATCH[1]}"
                if [[ ${seq} =~ \"$ ]]; then
                    seq=$(echo "${seq%\"}" | tr -d '[:space:]')
                    printf ">%s_protein%d\n%s\n" ${clusterid} ${counter} ${seq}
                    counter=$((counter + 1)); in_cds=0
                fi
            elif [[ -n ${seq} ]]; then
                seq="${seq}$(echo "${line}" | sed 's/"$//')"
                if [[ ${line} =~ \"$ ]]; then
                    seq=$(echo "${seq}" | tr -d '[:space:]')
                    printf ">%s_protein%d\n%s\n" ${clusterid} ${counter} ${seq}
                    counter=$((counter + 1)); in_cds=0; seq=""
                fi
            fi
        fi
    done < ${gbk}
}

# First /product of a feature type in a GenBank file
product() {
    awk -v feat="$2" '$1 == feat {f = 1} f && /\/product="/ {match($0, /\/product="([^"]+)"/, a); print a[1]; exit}' "$1"
}

process_dataset() {
    local dataset=$1 src=$2
    local gbkdir=${TMP}/${dataset}_GBK faadir=${TMP}/${dataset}_FAA
    local summary=${OUTDIR}/${dataset}_summary_info.csv
    mkdir -p ${gbkdir} ${faadir}
    > ${summary}

    echo "${dataset}: $(find ${src} -type f -name "*region*.gbk" -printf '%h\n' | sort -u | wc -l) of $(ls ${src} | wc -l) genomes with BGCs" >> ${LOG}

    for d in ${src}/*/; do
        genome=$(basename ${d})
        n=1
        for g in ${d}/*region*.gbk; do
            [ -f "${g}" ] || continue
            cp "${g}" ${gbkdir}/${genome}__BGC${n}.gbk
            n=$((n + 1))
        done
    done

    for f in ${gbkdir}/*.gbk; do
        name=$(basename ${f} .gbk)
        contig=$(awk '/^LOCUS/ {print $2; exit}' ${f})
        len=$(awk '/^LOCUS/ {print $3; exit}' ${f})
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" ${dataset} ${name} $(echo ${name} | tr '_' '-') ${contig} ${len} \
            "$(product ${f} region)" "$(product ${f} cand_cluster)" "$(product ${f} protocluster)" "$(product ${f} proto_core)" >> ${summary}
        if (( len >= MINLEN )); then
            extract_cds ${f} > ${faadir}/${name}.faa
        fi
    done

    cat ${faadir}/*.faa > ${OUTDIR}/${dataset}_bgcs_min${MINLEN}.faa
    echo "${dataset}: $(ls ${faadir} | wc -l) of $(ls ${gbkdir} | wc -l) BGCs >= ${MINLEN} bp" >> ${LOG}
}

# Sediment and seawater MAGs
process_dataset IMOS_SED ${ANTISMASH_DIR}/IMOS_SED
process_dataset IMOS_SEA ${ANTISMASH_DIR}/IMOS_SEA

# Dereplicated OMD genomes only
mkdir -p ${TMP}/OMD_DEREP
while read -r path; do
    name=$(basename ${path} .fa)
    cp -r ${ANTISMASH_DIR}/OMD_DB/${name}_${MINLEN} ${TMP}/OMD_DEREP/
done < ${OMD_REPS}
process_dataset OMD_DEREP ${TMP}/OMD_DEREP

# Reef Microbiomics Database (RMD): BGCs are taken from the antiSMASH output
# published with the RMD (Wiederkehr et al., 2026, Nature; rmd.microbiomics.io),
# restricted to the dereplicated RMD genomes, and processed as above
# (RMDv1_derep-antismash_min5000.faa). TODO: commands.

# MIBiG v4.0: rename to MIBIG<BGC>_<n>
awk -F'|' '/^>/ {c = substr($1, 2); n[c]++; print ">MIBIG" c "_" n[c]; next} {print}' ${MIBIG_FAA} \
    > ${OUTDIR}/mibig_renamed.fasta

cat ${LOG}
