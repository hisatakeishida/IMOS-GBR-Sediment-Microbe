
# Information for samples, sequencing, and quality control 

## Table of Contents
0. [Sample and Sequencing information all](#sample_info)
1. [QC information](#qc_info)
2. [Initial diversity investigation using SingleM] (#div_classify)

## 0. Sample and Sequencing information all samples <a name="sample_info"></a>
-  Illumina Novaseq (deep sequencing:30-70 Gbp)and Oxford Nanopore sequencing on 27 samples across the Great Barrier Reef in 2019 and 2020
- If sponge microbiome, expect 80-100 HQ MAGs from ~40 Gbp sequencing
- All sites associated with LTMP (long-term monitoring program, reef, algae cover, nutrient etc)
- Surface seidments were collected from the boat using Grab-type approach
- Illumina sample code: Site_replicate#__SequenceSampleID_R1/2_allRunsCombined
- Nanopore sample code: Site_replicate#__SequenceSampleID_NanoporeSampleNumber_barcode_allRunsCombined
- We are not doing co-assembly
- All sequenced at ACE
- Expecting 10-15 MAGs per sample, thus 250-300 MAGs in total?

## 1. Quality-Check <a name="qc_info"></a>
Nanopore 
- porechop

Illumina 
- fastqc 
- We may consider downsampling to to 60 Gbp if above 60Gbp

Seqkit stats for both nanopore and illumina

'''
conda activate porechop_0.2.4
porechop -t <threads> -i <input_file> -o <output_file>_porechopped.fastq.gz --discard_middle --min_split_read_size 250

conda activate seqkit_2.5.1
seqkit stats -j 24 --all -T -o /scratch/project/ncmas_d85/ishida/IMOS/deep_samples/nanopore/seqkit_results_preporechop *.gz

conda activate qc 
fastqc *.fastq.gz -o /scratch/project/ncmas_d85/ishida/IMOS/deep_samples/illumina/fastqc/ -t 24 
'''

## 2. Microbial Community Profiling using SingleM pipe <a name="div_classify"></a>
- Search for GTDB OTU data in metagenomic reads using SingleM
