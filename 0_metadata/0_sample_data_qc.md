# Samples, sequencing, and quality control

## Table of contents
1. [Sample information](#sample_info)
2. [Sequencing information](#seq)
3. [Raw read processing](#read_processing)

## 1. Sample information <a name="sample_info"></a>
- 192 surface sediment samples (48 sites × 4 replicates) were collected across the GBR between November 2019 and July 2020, with the AIMS Long-Term Monitoring Program (LTMP)
  - Austral summer (34 sites)
    - Northern GBR: Trip 1 (November to December 2019)
    - Southern GBR: Trip 2 (January 2020)
    - Central GBR: Trip 3 (February 2020)
  - Austral winter (14 sites)
    - Central GBR: Trip 4 (July 2020)
- Zoning status: 23 no-take marine reserves (NTMRs) and 25 fished reefs
- 191 samples were sequenced; DNA from one sample (reef 11-162) did not meet sequencing requirements
- Sample code: `Site_replicate#_SequenceSampleID`
- Full metadata: Table S1
- [`Field photos`](0_metadata/1_Sampling): photos provided by Sara Bell

## 2. Sequencing information <a name="seq"></a>
- Cell separation: microbial cells were detached from sediment using a cell separation buffer, sonication, and low-speed centrifugation, prefiltered (8 µm and 5 µm), and collected on 0.22 µm Sterivex filters
- DNA extraction: phenol:chloroform:isoamyl alcohol protocol
- Sequencing design:
  - For 27 of the 48 sites, three replicates were sequenced at standard depth (target 5 Gbp; achieved 2 to 10 Gbp) and one replicate at high depth (target 40 Gbp; achieved 30 to 70 Gbp)
  - The same DNA extract used for deep sequencing was also used for ONT long-read sequencing; 17 of 27 samples yielded ONT data of sufficient quality

### Illumina
- Library prep: Nextera DNA Flex Library Preparation Kit (Illumina #20018705)
- Sequencing: Illumina NovaSeq 6000, S4 flow cells, 2 × 150 bp
- Facilities: Microba Life Sciences and Centre for Microbiome Research (QUT) for standard-depth samples; Australian Centre for Ecogenomics (UQ) for deep samples

### ONT
- Size selection: Circulomics SRE XS kit (PacBio, SKU 102-208-200) to remove fragments < 10 kb
- Library prep: Ligation Sequencing Kit (Oxford Nanopore, LSK-109) with PCR-Barcoding Expansion Kit (EXP-PBC001); two samples per flow cell
- Sequencing: ONT PromethION 24, R9.4 flow cells, 72 h, MinKNOW v20.06.18 (default settings), at the Australian Centre for Ecogenomics (UQ)
- Basecalling: Guppy v5.0.16 (super-accuracy model)

## 3. Raw read processing <a name="read_processing"></a>
Adapters were removed with fastp v0.23.4 (Illumina) and Porechop v0.2.4 (ONT). 
Read yields are summarised in Table S3.

```bash

# Raw stats 
conda activate seqkit_2.5.1 

seqkit stats -T <sample>_R1.fastq.gz <sample>_R2.fastq.gz \
--threads <threads> > <sample>_Illumina_raw_seqkit.tsv

seqkit stats -T <sample>.ONT.fastq.gz \
--threads <threads> > <sample>_ONT_raw_seqkit.tsv

# Illumina: adapter removal only (-Q disables quality filtering)
conda activate fastp_0.23.4

fastp -Q \
  --in1 <sample>_R1.fastq.gz --in2 <sample>_R2.fastq.gz \
  --out1 <sample>_R1.adprem.fastq.gz --out2 <sample>_R2.adprem.fastq.gz \
  --detect_adapter_for_pe \
  --thread <threads> \
  --html <sample>.fastp.html --json <sample>.fastp.json

# ONT: adapter removal and splitting of chimeric reads using Porechop 
# https://github.com/rrwick/Porechop.git

conda activate porechop_0.2.4

porechop -t <threads> \
  -i <sample>.ONT.fastq.gz \
  -o <sample>.ONT_porechopped.fastq.gz \
  --discard_middle --min_split_read_size 250

# Post QC stats 
conda activate seqkit_2.5.1 

seqkit stats -T <sample>_R1.adprem.fastq.gz <sample>_R2.adprem.fastq.gz \
--threads <threads> > <sample>_Illumina_adprem_seqkit.tsv

seqkit stats -T <sample>.ONT.fastq.gz \
--threads <threads> > <sample>_ONT_adprem_seqkit.tsv

# Printing
for i in *_ONT_adprem_seqkit.tsv
do
    samp_name=$(basename "$i" _ONT_adprem_seqkit.tsv)
    # sample=$(echo "$samp_name" | awk -F '__' '{print $1}')
    numreads=$(tail -n1 $i | awk '{print $4}')
    sumreads=$(tail -n1 $i | awk '{print $5}')

    echo $samp_name $numreads $sumreads
    # break
done

```