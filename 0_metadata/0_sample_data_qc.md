# Information for samples, sequencing, and quality control 

## Table of Contents
1. [Sample information](#sample_info)
2. [Sequencing information](#Seq)
3. [Raw read processing](#Adapter_trimming)

## 1. Sample information <a name="sample_info"></a>
- 192 sediment samples (48 sites * 4 replicates) were collected across the GBR in 2019 and 2020
  - During the austral summer, covering 34 sites
    - the northern GBR (Trip 1 in November-December 2019)
    - the southern GBR (Trip 2 in January 2020)
    - the central GBR (Trip 3 in February 2020). 
  -  During austral winter, covering the remaining 14 sites 
    - the central GBR (Trip 4 in July 2020)
- All sites are associated with AIMS long-term monitoring program (LTMP)
- Sample code: Site_replicate#_SequenceSampleID

## 2. Sequencing information <a name="Seq"></a>
- DNA extraction: Phenol:chloroform:isoamyl alcohol protocol 
- For 27 sites (out of 48), quadruplicates were split into 3 replicates for shallow sequencing (2-10 Gbp) and 1 replicate for deep sequencing (30-70 Gbp)
- For these 27 sites, same extracted genomic DNA from the deep sequencing, were used for ONT long-read sequencing

- Illumina 
  - Library prep: Nextera DNA Flex Library Preparation Kit (Illumina # 20018705)
  - Sequencing platform: Illumina NovaSeq 6000 S4 2*150 bp flow cells at MICROBA and QUT 
- ONT
  - Circulomics SRE XS kit (PacBio, SKU 102-208-200) for removing short fragments
  - Library prep: Ligation Sequencing kit (Oxford Nanopore, LSK-109) with PCR-Barcoding Expansion Kit (EXP-PBC001)
  - Sequencing platform: ONT PromethION 24 using R9.4 flow cells for 72 hours with MinKNOW (v20.06.18) under default settings
  - Guppy v5.0.16 for rebasecalling

## 2. Adapter trimming <a name="Adapter_trimming"></a>
- Fastp v0.23.4 and Porechop v0.2.4 were used for Illumina and Nanopore datasets, respectively

```
#fastp for Illumina
conda activate fastp_0.23.4
fastp -Q --in1 /QRISdata/Q6440/imos_sediment_reads_illumina/lynchs_1__21S000275_R1_allRunsCombined.fastq.gz --in2 

#porechop for ONT data
conda activate porechop_0.2.4
porechop -t <threads> -i <input_file> -o <output_file>_porechopped.fastq.gz --discard_middle --min_split_read_size 250
```
