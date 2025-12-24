# Community-level fucntional profiling based on reconstruction of non-redundunt protein database

## AIM: Investigate functional diversity using protein db (KO and CAZyme)

## Table of Contents
1. [Assembly](#ass)
2. [Gene prediction](#gene_pre)
3. [Dereplication](#derep)
4. [Gene annotation](#gene_annot)
5. [Read mapping](#readmapping)
6. [RPKM preparation](#rpkm_prep)
7. [RPKM and Visualization](#rpkm)

## 1. Assembly <a name="ass"></a>
- Assemble metagenomic reads into contigs 
- See **[Metagenomic assembly](/2_genome_centric/1_genome_recovery/1_genome_recovery.md)**
 
## 2. Gene prediction <a name="gene_pre"></a>
- Using prodigal 2.6.3 for gene prediction 
- Platform: Bash
- IN: Assembled contigs in fasta file (INDIR)
- OUTDIR: includes predicted protein sequences  

```
module load prodigal/2.6.3-gcccore-11.3.0

prodigal -i ${sample_name}.final_contigs.fasta -f gff -p meta -o ${sample_name}.prodigal.gff -a ${sample_name}.protein.translations.faa -d ${sample_name}.nucleotide.faa -q
```

## 3. Dereplication <a name="derep"></a>
- Using mmseqs2 v15.6f452 (https://github.com/soedinglab/MMseqs2) to dereplicate all predicted protein sequences from assembled contigs 
- Platform: Bash
- IN: All protein predicted sequences 
- OUT: Dereplicated protein sequences

```
module load mmseqs2/15-6f452

cat *.protein.translations.faa > all_assembly_prodigal.proteins.faa

mmseqs easy-linclust all_assembly_prodigal.proteins.faa all_assembly_prodigal.proteins.faa.prot /scratch/temp/${SLURM_JOB_ID} --min-seq-id 0.8 -c 0.8 --threads ${SLURM_CPUS_PER_TASK} --split-memory-limit 400G 
```

### Dataset Statistics
| Data              | Type    | Num Seqs      | Sum Length        | Min Len | Avg Len | Max Len |
|-------------------|---------|---------------|-------------------|---------|---------|---------|
| Before clustering | Protein | 2,036,382,959 | 170,506,768,952   | 20      | 83.7    | 20,426  |
| After clustering  | Protein | 1,682,992,275 | 137,421,747,632   | 20      | 81.7    | 20,426  |

## 4. Gene annotation <a name="gene_annot"></a>
- Annotating the dereplicated protein sequences with KO (https://github.com/WrightonLabCSU/DRAM)
- Platform: Bash
- IN: Dereplicated protein sequences
- OUT: KO and CAZyme annotation 

```
module load dram/1.5 

infile=dereplicated_all_assembly_prodigal_proteins.faa

DRAM.py annotate_genes -i "${infile}" -o dereplicated_all_assembly_prodigal_proteins_annotation --threads ${SLURM_CPUS_PER_TASK} 

```

## 5. Read mapping <a name="readmapping"></a>
- Map reads to the proteins with annotation and singlem diamond database
- Platform: Bash

```
#################################
# 1. Find proteins with anotation
#################################

annot_file=annotations.tsv #from dram 
seq_file=genes.faa #from dram 
outfile_ko="ko.faa"
outfile_cazy="cazy.faa"

module load seqtk 

# kegg
awk -F'\t' 'NR > 1 && $4 != ""' $annot_file | cut -f 1  > ko_name.txt 
seqtk subseq $seq_file ko_name.txt > ${outfile_ko}

# cazy
awk -F'\t' 'NR > 1 && $17 != ""' $annot_file | cut -f 1  > cazy_name.txt 
seqtk subseq $seq_file cazy_name.txt  >> ${outfile_cazy}

pigz ${outfile_ko} ${outfile_cazy}

############################
# 2. Diamond db construction
############################

module load diamond/2.1.0-gcc-11.3.0

diamond makedb --in ${outfile_ko} -d ${outfile_ko%.faa.gz} --threads ${SLURM_CPUS_PER_TASK}
diamond makedb --in ${outfile_cazy} -d ${outfile_cazy%.faa.gz} --threads ${SLURM_CPUS_PER_TASK} 
diamond makedb --in SingleM_RS220_allmarkers.faa -d SingleM_RS220_allmarkers --threads ${SLURM_CPUS_PER_TASK}

#######################################
# 3. Run diamond with alignment filters
#######################################

diamond blastx --threads ${SLURM_CPUS_PER_TASK} -d ${outfile_ko%.faa.gz}.dmnd -q ${sample_name}.R1.140bp.10Mseqs.fq.gz -o ${sample_name}.ko.tsv --evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1

diamond blastx --threads ${SLURM_CPUS_PER_TASK} -d ${outfile_cazy%.faa.gz}.dmnd -q ${sample_name}.R1.140bp.10Mseqs.fq.gz -o ${sample_name}.cazy.tsv --evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1

diamond blastx --threads ${SLURM_CPUS_PER_TASK} -d SingleM_RS220_allmarkers.dmnd -q ${sample_name}.R1.140bp.10Mseqs.fq.gz -o ${sample_name}.singlem.tsv --evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1

```

## 6. Matrix and gene-length calculation <a name="rpkm_prep"></a>
- Parse count from diamond into count matrix and calculating gene-lengths
- Platform: Bash

```
#######################################
# 1. Create count matrix in long format 
#######################################

in2=diamond_counts_dir
out=diamond_long_format.tsv

echo -e "Sample\tGeneID\tCount" > "$out"

for file in ${in2}/*.kegg.tsv; do
# for file in ${in2}/*.cazy.tsv; do
    sample=$(basename "$file" .kegg.tsv)
    #sample=$(basename "$file" .cazy.tsv)
    echo $sample 

    # Extract column 2 (GeneID), count occurrences
    awk -v sample="$sample" -F'\t' '{print $2}' "$file" | sort | uniq -c | awk -v sample="$sample" '{print sample "\t" $2 "\t" $1}' >> "$out"
done 

#################################################
# 2. Calculate gene length for RPKM normalization
#################################################

cat ko.faa | awk '/^>/{if(seqlen){print gene "\t" seqlen*3}; gene=substr($0,2); seqlen=0; next}{seqlen+=length($0)} END{print gene "\t" seqlen*3}' > ko.seqlen.tsv

cat cazy.faa | awk '/^>/{if(seqlen){print gene "\t" seqlen*3}; gene=substr($0,2); seqlen=0; next}{seqlen+=length($0)} END{print gene "\t" seqlen*3}' > cazy.seqlen.tsv

cat SingleM_RS220_allmarkers.faa | awk '/^>/{if(seqlen){print gene "\t" seqlen*3}; gene=substr($0,2); seqlen=0; next}{seqlen+=length($0)} END{print gene "\t" seqlen*3}' > SingleM_RS220_allmarkers.seqlen.tsv

######################
# 3. Get gene taxonomy 
######################

genename_path=ko.name.tsv
echo -e "gene_id\tko_id" >  $genename_path
tail -n +2 annotations.tsv | cut -f1,4  >>  $genename_path

genename_path=cazy.name.tsv
echo -e "gene_id\tcazy_id" >  $genename_path
tail -n +2 annotations.tsv | cut -f1,17 >>  $genename_path

```

## 6. RPKM and Visualization <a name="rpkm"></a>
- calculate RPKM for protein and singlem db with NMDS analysis 
- Platform: R
- INFILEs
    - Count matrix in long format (tsv)
    - metadata (xlsx)
    - annotation table (tsv)
    - lengths of protein (tsv)
    - count matrix of SingleM markers (tsv)
    - lengths of singlem protein markers (tsv)
- OUT: NMDS ordination with aitchison distance, library size, rarecurve
- **[KO Scripts: 3_RPKM_calculation_KO.R](/5_Rscripts/3_RPKM_calculation_KO.R)**
- **[CAZyme Scripts: 4_RPKM_calculation_CAZyme.R](/5_Rscripts/4_RPKM_calculation_CAZyme.R)**


