# Zenodo data release

Data files deposited on Zenodo (DOI: `https://doi.org/10.5281/zenodo.24187765`) for:

> Ishida, H., Robbins, S. J., Bell, S. C., Lim, J. Y.,  Terzin, M., Webster, N., Bourne, D. G., Laffy, P. W., Hugenholtz, P., Yeoh, Y. K., & Chan, C. X. **Comparative metagenomics reveals diversity and metabolic heterogeneity of sediment microbes from the Great Barrier Reef.** *In preparation.*

Code that produced each file is in this repository; the step numbers below refer to its scripts. Raw reads are on the [AIMS data repository](https://doi.org/10.25845/kkma-n127) and NCBI BioProject accession [PRJEB112104](https://www.ebi.ac.uk/ena/browser/view/PRJEB112104).

Seawater MAGs (MGD-seawater) are available from Zenodo record [20503966](https://zenodo.org/records/20503966).

## Datasets
- **Gene catalogue**: non-redundant proteins from all sediment and seawater assemblies (complete ORFs; MMseqs2 easy-linclust, 90% identity, 80% coverage)
- **MGD-sediment**: MAGs from 191 GBR surface sediment metagenomes (48 reefs, four sampling trips, November 2019 to July 2020)
- **MGD-seawater**: MAGs from the co-located surface seawater metagenomes (876 species representatives; Terzin et al., 2026; Robbins et al., 2026)

## Files

### Gene catalogue
| File | Contents | Steps |
|---|---|---|
| `01_gene-catalogue_proteins.faa.gz` | Non-redundant protein catalogue (complete ORFs, 90% identity, 80% coverage) | 2.01, 2.02 |
| `02_gene-catalogue_protein-lengths.tsv.gz` | Length of each catalogue protein in nucleotides (`GeneID`, `Length`) | 2.10 |
| `03_gene-catalogue_KO-annotations.tsv.gz` | KOfam annotation (DRAM v1.5, E-value 1e-6) of catalogue proteins (`GeneID`, `ko_id`) | 2.03, 2.08 |
| `04_gene-catalogue_other-annotations.tar.gz` | CAZyme (run_dbcan v5.0.1), HydDB, metabolic marker gene, and *amoA* annotations of catalogue proteins | 2.04 to 2.08 |
| `05_gene-catalogue_read-counts.tsv.gz` | Reads mapped to each catalogue protein per sample, long format (`Sample`, `GeneID`, `Count`); DIAMOND blastx of R1 reads ≥ 140 bp per sample | 2.11, 2.12 |
| `06_singlem-markers_read-counts.tsv.gz` | Reads mapped to each SingleM marker gene per sample, long format (`Sample`, `GeneID`, `Count`); same reads as file 05 | 2.11, 2.12 |
| `07_singlem-markers_lengths.tsv.gz` | Length of each SingleM marker gene in nucleotides (`GeneID`, `Length`) | 2.10 |
| `08_gene-catalogue_phyloseq-objects.tar.gz` | R phyloseq objects of normalised RPKM (catalogue protein RPKM / mean SingleM marker gene RPKM per sample): `dram_ko_e6.RPKM.rds`, `ps_ko.RDS`, `ps_dbcan.RDS`, `ps_greening_hyddb.RDS` | 2.13, 2.14 |

### MAGs
| File | Contents | Steps |
|---|---|---|
| `09_MGD-sediment_2698-MAGs.tar.gz` | All 2,698 sediment MAGs with quality score ≥ 50 (CheckM or CheckM2), from single-sample assemblies and co-assemblies | 4.01 to 4.07 |
| `10_MGD-sediment_1201-MAGs_species-representatives.tar.gz` | 1,201 species-representative sediment MAGs (dereplicated at 95% ANI) | 4.08 |

### MAG annotations
| File | Contents | Steps |
|---|---|---|
| `11_MGD-sediment_1201-MAGs_DRAM-annotations.tar.gz` | DRAM v1.5 annotations of the 1,201 sediment MAGs (`annotations.tsv`, `genes.faa`, `genes.fna`, `genes.gff` per MAG) | 5.01 |
| `12_MGD-sediment-seawater_2077-MAGs_GreeningDB-annotations.tar.gz` | Metabolic marker gene (Greening lab database) hits for the 1,201 sediment and 876 seawater MAGs | 5.04 |
| `13_MGD-sediment_1201-MAGs_3504-BGCs_proteins.faa.gz` | Proteins of 3,504 biosynthetic gene clusters (antiSMASH 6.1.1, BGCs ≥ 5,000 bp) from sediment MAGs | 5.09, 5.10 |
| `14_MGD-seawater_876-MAGs_1472-BGCs_proteins.faa.gz` | Proteins of 1,472 biosynthetic gene clusters (BGCs ≥ 5,000 bp) from seawater MAGs | 5.09, 5.10 |
| `15_FeFe-hydrogenases_MAGs-and-gene-catalogue_min200aa.faa.gz` | [FeFe] hydrogenase catalytic subunits (≥ 200 aa) from sediment and seawater MAGs and the gene catalogue | 2.05, 5.04 |
| `16_NiFe-hydrogenases_MAGs-and-gene-catalogue_min200aa.faa.gz` | [NiFe] hydrogenase large subunits (≥ 200 aa) from sediment and seawater MAGs and the gene catalogue | 2.05, 5.04 |

## Using the files

```bash
# Archives of per-MAG archives: extract the outer archive, then the MAG of interest
tar -xzf 10_MGD-sediment_1201-MAGs_species-representatives.tar.gz
tar -xzf <MAG_ID>.tar.gz

# Compressed tables and FASTA files can be read directly
zcat 05_gene-catalogue_read-counts.tsv.gz | head
```

```r
# Phyloseq objects (R packages phyloseq and microViz)
library(phyloseq)
ps_ko <- readRDS("ps_ko.RDS")
```

Normalised RPKM can be recalculated from files 02, 03, 05, 06, and 07 with [`2.13_protein_RPKM.R`](../1_gene_centric/2_Protein/scripts/2.13_protein_RPKM.R).

## Checksums
MD5 checksums of all files are in `md5sums.txt`:

```bash
md5sum -c md5sums.txt
```


## Contact
Hisatake Ishida, hisatake.ishida@uq.net.au
