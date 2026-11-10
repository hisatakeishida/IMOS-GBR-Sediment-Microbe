#!/bin/bash --login
#SBATCH --job-name=unbinned_contigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o unbinned_contigs_%A_%a.o
#SBATCH -e unbinned_contigs_%A_%a.e
#SBATCH --array=1-191

# MAG recovery, step 14: taxonomy of unbinned contigs, weighted by reads.
#   1. Unbinned contigs: contigs of each single-sample assembly not in any of
#      its MAGs
#   2. Taxonomy with MMseqs2 v15.6f452 easy-taxonomy against a custom database
#      (RefSeq, EukProt, PLSDB, IMG/VR, and proteins predicted from GTDB
#      R08-RS214 genomes; Robbins et al., 2026)
#   3. Map MAG-unmapped reads (4.13) from all samples of the same reef to the
#      unbinned contigs (CoverM contig, >= 95% identity, >= 75% aligned)
#   4. Sum mapped reads per domain
#
# Assembly names are <reef>_<replicate>[...], e.g. 11049_1, st_Crispin_2.

set -euo pipefail

ASSEMBLY_LIST="<path/to/single_assembly_contigs_list.txt>"   # one <sample>.fasta.gz per line
MAG_DIR="<path/to/sediment_MAGs>"                            # MAG file names start with <sample>
UNMAPPED_DIR="<path/to/4_MAG_recovery/11_unbinned_contigs/unmapped_reads>"   # output of 4.13
MMSEQS_DB="<path/to/mmseqs_database/DB_LCA>"
OUTDIR="<path/to/4_MAG_recovery/11_unbinned_contigs>"

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
mkdir -p ${OUTDIR}/{taxonomy,mapping,domain_summary}

contigs=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ASSEMBLY_LIST})
sample=$(basename ${contigs} | cut -f 1 -d ".")
reef=$(echo ${sample} | sed -E 's/_[0-9].*$//')
echo "Assembly: ${sample}; reef: ${reef}"

# ---- 1. Unbinned contigs ------------------------------------------------------------
zcat ${contigs} | grep "^>" | sed 's/^>//; s/ .*//' | sort > ${TMP}/all.txt
cat ${MAG_DIR}/${sample}[_.]*.fna 2>/dev/null | grep "^>" | sed 's/^>//; s/ .*//' | sort > ${TMP}/binned.txt || true
comm -23 ${TMP}/all.txt ${TMP}/binned.txt > ${TMP}/unbinned.txt
echo "Contigs: $(wc -l < ${TMP}/all.txt); unbinned: $(wc -l < ${TMP}/unbinned.txt)"

module load seqtk
seqtk subseq ${contigs} ${TMP}/unbinned.txt > ${TMP}/${sample}_unbinned.contigs.fa

# ---- 2. Taxonomy --------------------------------------------------------------------
module load mmseqs2   # MMseqs2 v15.6f452
mmseqs easy-taxonomy ${TMP}/${sample}_unbinned.contigs.fa ${MMSEQS_DB} \
    ${OUTDIR}/taxonomy/${sample}.taxonomy ${TMP}/mmseqs_tmp \
    --tax-lineage 1 \
    --threads ${SLURM_CPUS_PER_TASK} --split-memory-limit 100G

# ---- 3. Map MAG-unmapped reads from the same reef -------------------------------------
conda activate coverm   # CoverM v0.7.0
for R1 in ${UNMAPPED_DIR}/${reef}_[0-9]*_unmapped_R1_allRunsCombined.fastq.gz; do
    readsample=$(basename ${R1} _unmapped_R1_allRunsCombined.fastq.gz)
    coverm contig -m count \
        -r ${TMP}/${sample}_unbinned.contigs.fa \
        -1 ${R1} -2 ${R1/_R1_/_R2_} \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        -t ${SLURM_CPUS_PER_TASK} \
        -o ${OUTDIR}/mapping/${sample}__${readsample}.contigs.mappings.txt
done

# ---- 4. Reads per domain --------------------------------------------------------------
# Reads per contig summed over all mapped samples, joined with contig domain
# (lineage column of the LCA output; "Unclassified" if none)
awk -F'\t' -v OFS='\t' 'FNR > 1 {count[$1] += $2} END {for (c in count) print c, count[c]}' \
    ${OUTDIR}/mapping/${sample}__*.contigs.mappings.txt | sort -k1,1 > ${TMP}/counts.tsv

awk -F'\t' -v OFS='\t' '{
    domain = "Unclassified"
    n = split($9, tax, ";")
    for (i = 1; i <= n; i++) if (tax[i] ~ /^d_/) {domain = tax[i]; break}
    print $1, domain
}' ${OUTDIR}/taxonomy/${sample}.taxonomy_lca.tsv | sort -k1,1 > ${TMP}/domain.tsv

join -t $'\t' -a 1 -e "Unclassified" -o 1.1,1.2,2.2 ${TMP}/counts.tsv ${TMP}/domain.tsv | \
    awk -F'\t' -v OFS='\t' -v s=${sample} '{reads[$3] += $2; total += $2}
        END {for (d in reads) print s, d, reads[d], 100 * reads[d] / total}' | \
    sed '1i Assembly\tDomain\tReads\tPercent' > ${OUTDIR}/domain_summary/${sample}.domain_reads.tsv

