
# Information for samples, sequencing, and quality control 

## Table of Contents
0. [Sample and Sequencing information all](#sample_info)
1. [QC information](#qc_info)
2. [Read pre-processing](#prepros)
[Taxa of interest](#taxa)

## 0. Sample and Sequencing information all samples <a name="sample_info"></a>
- Illumina Novaseq (deep sequencing:30-110 Gbp)and PCR-coupled Oxford Nanopore sequencing on 27 samples across the Great Barrier Reef in 2019 and 2020
- If sponge microbiome, expect 80-100 HQ MAGs from ~40 Gbp sequencing
- All sites associated with LTMP (long-term monitoring program, reef, algae cover, nutrient etc)
- Surface sediments were collected from the boat using Grab-type approach
- Illumina sample code: Site_replicate#__SequenceSampleID_R1/2_allRunsCombined
- Nanopore sample code: Site_replicate#__SequenceSampleID_NanoporeSampleNumber_barcode_allRunsCombined
- We are not doing co-assembly
- All sequenced at ACE
- Expecting 10-15 MAGs per sample, thus 250-300 MAGs in total?
- PCR information for Nanopore 
  - Long range PCR bias
  - Enzyme or transposn based
  - Need microgram for allegation but only got ~10 nannograms from extraction
  - Removed mid-strand adapters (since we used two flow cell) 

## 1. Quality-Check <a name="qc_info"></a>
Seqkit stats for both nanopore and illumina

```conda activate seqkit_2.5.1
seqkit stats -j 24 --all -T -o <output_file> *.gz

conda activate qc 
fastqc *.fastq.gz -o <output_file.gz> -t 24 
```

## 2. Read pre-processing <a name="prepros"></a>
- Check metadata for before and after quality control
  
### 2a. Nanopore 
- porechop
```
conda activate porechop_0.2.4
porechop -t <threads> -i <input_file> -o <output_file>_porechopped.fastq.gz --discard_middle --min_split_read_size 250
```
- we will be comparing Illumina only vs Illumina+Nanopore with samples with large nanopore reads (11049_4__SE6286_Pr130_02_BC06, mcsweeny_4__SE6278_Pr105_2_BC02, mantis_4__SE6289_Pr109_BC11,22084_4__SE6293_Pr132_01_BC09,grub_4__SE6285_Pr130_01_BC05,11162_4__SE6287_Pr108_BC10,lagoon_4__SE6294_Pr132_02_BC10, broomfield_4__SE6292_Pr109_BC12)

### 2b. Illumina 
- We downsampled reads to to 50 Gbp if above 50Gbp (R1 and R2 base number combined) for more manageable metagenomic assembly
```
module load seqtk
seqtk sample -s 11 <input_file>  <fraction0-1>  | gzip > <output_file.gz>
```

| Samples               | Downsample y_n |
|-----------------------|----------------|
| 11049_4__SE6259       | y              |
| 11162_4__SE6260       | y              |
| 13124_4__SE6261       | y              |
| 21550_4__SE6270       | y              |
| 21580_4__SE6274       | y              |
| 22084_4__SE6266       | y              |
| agincourt1_4__SE6271  | y              |
| arlington_4__SE6272   | y              |
| broomfield_4__SE6265  | y              |
| chicken_4__SE6257     | n              |
| chinaman_4__SE6250    | y              |
| corbett_4__SE6253     | y              |
| davie_4__SE6268       | y              |
| erskine_4__SE6264     | y              |
| grub_4__SE6258        | y              |
| hoskyn_4__SE6269      | y              |
| kelso_4__SE6255       | y              |
| lagoon_4__SE6267      | y              |
| mantis_4__SE6262      | y              |
| masthead_4__SE6263    | y              |
| mcsweeny_4__SE6251    | n              |
| moonsoon_4__SE6252    | y              |
| myrmidon_4__SE6256    | n              |
| north_4__SE6276       | y              |
| sandbank1_4__SE6254   | n              |
| smallLagoon_4__SE6275 | y              |
| taylor_4__SE6273      | y              |


Go to next analysis step [deep_1_read_analysis_assembly.md](deep_1_read_analysis_assembly.md)

## Taxa of interest <a name="taxa"></a>
- 4484-113
- UBP17
- UBA9160
