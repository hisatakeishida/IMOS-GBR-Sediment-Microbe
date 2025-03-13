# Information for samples, sequencing, and quality control 

## Table of Contents
0. [Sample and Sequencing information all](#sample_info)
1. [Raw read processing](#Adapter_trimming)

## 0. Sample and Sequencing information all samples <a name="sample_info"></a>
- 192 sediment samples (48 sites * 4 replicates) were collected across the GBR in 2019 and 2020
- All sites are associated with LTMP (long-term monitoring program, reef, algae cover, nutrient etc)
- For 27 sites (out of 48), quadruplicates were split into 3 replicates for shallow sequencing (2-10 Gbp) and 1 replicate for deep sequencing (30-70 Gbp)
- For these 27 sites, same extracted genomic DNA from the deep sequencing, were used for ONT long-read sequencing, with PCR-ligation kit
- Sample code
  -   Site_replicate#_SequenceSampleID
- Illumina Novaseq sequencing
- Sequenced at MICROBA, ACE, QUT

## 1. Adapter trimming <a name="Adapter_trimming"></a>
- Fastp and Porechop were used for Illumina and Nanopore datasets, respectively

- porechop for Nanopore
'''
conda activate porechop_0.2.4
porechop -t <threads> -i <input_file> -o <output_file>_porechopped.fastq.gz --discard_middle --min_split_read_size 250
'''
