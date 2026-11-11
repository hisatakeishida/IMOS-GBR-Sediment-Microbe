#!/bin/bash --login
#SBATCH --job-name=clustomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --account=<account>
#SBATCH --partition=general
#SBATCH -o clustomatic_%j.o
#SBATCH -e clustomatic_%j.e

# MAG annotation, step 11: group BGCs (>= 5 kbp) from sediment and seawater
# MAGs, dereplicated OMD and RMD genomes, and MIBiG v4.0 into gene cluster
# families (GCFs) with clust-o-matic (threshold 0.5).
#
# clust-o-matic: https://github.com/Helmholtz-HIPS/clustomatic_source
#   git clone https://github.com/Helmholtz-HIPS/clustomatic_source.git
#   pip3 install -r clustomatic_source/requirements.txt   # plus DIAMOND on PATH
# clust-o-matic runs an all-vs-all DIAMOND blastp of the BGC proteins, scores
# each pair of BGCs as 1 - 2 x (identical residues of best non-redundant protein
# hits) / (summed protein lengths of both BGCs), and groups BGCs by average-linkage
# clustering cut at THRESHOLD. Protein headers must be >BGCID_proteinID with a
# single underscore, hence the "_" to "-" conversion of BGC IDs in 5.10.
#
# Run with clustomatic.modified.py (in scripts/): the published clustomatic.py
# with an output folder (cluster_assignments.txt, distance_matrix.txt), DIAMOND
# threads and progress messages; scoring and clustering are unchanged. The
# original modified script was not kept, so the version in scripts/ is
# reconstructed from the published script and the run's log.

set -euo pipefail

BGC_DIR="<path/to/5_MAG_annotation/7_BGCs/bgc_proteins>"    # output of 5.10
RMD_FAA="<path/to/RMDv1_derep-antismash_min5000.faa>"        # dereplicated RMD BGCs (from the RMD antiSMASH output; Wiederkehr et al., 2026)
CLUSTOMATIC="<path/to/repo>/2_genome_centric/5_MAGsAnnot/scripts/clustomatic.modified.py"
OUTDIR="<path/to/5_MAG_annotation/7_BGCs>"
NAME=clustomatic__imos_sea_sed__derep_omd__derep_rmd__mibig_5000
THRESHOLD=0.5

TMP=${TMPDIR:-/scratch/temp/${SLURM_JOB_ID}}
cd ${TMP}

cat ${BGC_DIR}/IMOS_SEA_bgcs_min5000.faa \
    ${BGC_DIR}/IMOS_SED_bgcs_min5000.faa \
    ${BGC_DIR}/OMD_DEREP_bgcs_min5000.faa \
    ${RMD_FAA} \
    ${BGC_DIR}/mibig_renamed.fasta > ${TMP}/${NAME}.faa

conda activate clustomatic

python ${CLUSTOMATIC} ${TMP}/${NAME}.faa ${THRESHOLD} ${NAME} ${SLURM_CPUS_PER_TASK}

cp -r ${NAME} ${OUTDIR}/

