# Community profiling with SingleM

## AIM: Investigate prokaryotic diversity using marker-derived operational taxonomic units (mOTU) with SingleM v1.18.0 

## Table of Contents
1. [SingleM](#singlem_run)
2. [SingleM output parsing](#singlem_parse)
3. [SingleM output visualization](#singlem_visu)

## 1. SingleM <a name="singlem_run"></a>
- Using SingleM v1.18.0 pipe (https://github.com/wwood/singlem) to find mOTUs and summarise into profiles 
- Using SingleM v1.18.0 microbial fraction (https://github.com/wwood/singlem) to find prokaryotic fraction of the data
- Platform: Bash
- IN: Reads in fastq file (INDIR)
- OUT: Coverage matrix 

```
# READS
INDIR=""
OUTDIR=""

mkdir -p ${OUTDIR}

name=$(basename $FILENAME _R1_allRunsCombined.fastq.gz)
echo $name

mfrac_out=${OUTDIR}/mfrac
otu_out=${OUTDIR}/otu
prof_out=${OUTDIR}/prof

mkdir -p $mfrac_out $otu_out $prof_out

conda activate singlem 

#path to the metapackage 
metap=S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb

#pipe 
singlem pipe --metapackage $metap -1 $FILENAME -2 ${FILENAME%_R1_allRunsCombined.fastq.gz}_R2_allRunsCombined.fastq.gz --threads ${SLURM_CPUS_PER_TASK} -p ${prof_out}/${name}_prof.tsv --otu-table ${otu_out}/${name}_otu.tsv 

#microbial fraction
singlem microbial_fraction -1 $FILENAME -2 ${FILENAME%_R1_allRunsCombined.fastq.gz}_R2_allRunsCombined.fastq.gz -p ${prof_out}/${name}_prof.tsv --output-tsv ${mfrac_out}/${name}_microbfraction.tsv --metapackage $metap
```

## 2. SingleM output parsing <a name="singlem_parse"></a>
- Combine SingleM profiles from each sample into one matrix (row:tax and column:coverage)
- Platform: python
- IN: All the *_prof.tsv to combine
- OUT: matrix as csv (coverage)
- **[Scripts](/6_pythonscripts/singlem_prof_coverage.py)**

## 3. SingleM output visualization <a name="singlem_visu"></a>
- Visualise SingleM profile matrix 
- Platform: R
- IN: 
    - mOTU coverage matrix in relative abundance (for otu_table-phyloseq)
    - mOTU taxonomy table (for tax_table-phyloseq)
    - metadata (for sample_data-phyloseq)
- The script below will take the relative abundance and scale it based on library size and prokaryotic fraction estimated from SingleM microbial fraction
- OUT: hellinger-PCA, composition stacked bar plot, and core-family identification
- **[Scripts](/5_Rscripts/2_SingleM.R)**
