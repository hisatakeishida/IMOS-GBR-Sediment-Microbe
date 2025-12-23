# Alignment-free beta-diversity analyses with MASH (k-mer distance)

## AIM: Investigate diversity of k-mers, proxy for diversity of taxonomy and function without alignment

## Table of Contents
1. [Distance calculation](#dist_calc)
2. [Distance visualization](#dist_visu)


## 1. Calculate MASH-distance <a name="dist_calc"></a>
- Using (https://github.com/marbl/Mash) to calculate MASH distance between all pairs of samples 
- Platform: Bash 
- IN: Reads in fasta or fastq file (INDIR)
- OUT: Distance in text file (mash_dist.txt)

```
module load mash/2.3-gcc-10.3.0

#skeching 
for f in ${INDIR}/* ;do 
    mash sketch $f -p ${SLURM_CPUS_PER_TASK}
done 

#computing lower triangular distance matrix 
mash triangle -r -p ${SLURM_CPUS_PER_TASK} ${INDIR}/*.msh > mash_dist.txt
```

## 2. Visualise MASH-distance <a name="dist_visu"></a>
- Visualise MASH-distance matrix 
- Platform: R
- IN:**[Distance matrix](mash_dist.txt)** (text file) and **[metadata](MASH_ENV.xlsx)**
- OUT:plots in svg (NMDS1-NMDS2 or NMDS1-NMDS2-NMDS3), other stats (e.g. PERMANOVA, ANOSIM)
- **[Script](/5_Rscripts/1_MASH.R)**
