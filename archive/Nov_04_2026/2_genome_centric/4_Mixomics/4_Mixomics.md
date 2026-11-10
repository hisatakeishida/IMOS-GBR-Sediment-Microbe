# Can reef prokaryotic MAGs discriminate between no-take and fished reefs?

R code to test whether prokaryotic communities on the Great Barrier Reef (GBR) differ between reefs that are closed to fishing (no-take marine reserves, NTMRs) and reefs that are open to fishing, using metagenome-assembled genomes (MAGs) from the IMOS-MGD metagenomes.

The analysis follows one line of reasoning. We first look at the main clustering patterns without using zoning information, then test whether zoning explains a significant share of community variance once spatiotemporal factors are included, and only then use supervised methods (PLS-DA, sPLS-DA, MINT sPLS-DA) to identify indicator MAGs. Each supervised model is checked against a null built by shuffling zone labels among reefs, and with whole reefs or sectors held out.

The analysis is run separately for each habitat (`SEA` seawater, `SED` sediment).

---

## Contents

- [Study design](#study-design)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Configuration and running](#configuration-and-running)
- [Analysis workflow](#analysis-workflow)
  - [Data preparation and CLR transformation](#data-preparation-and-clr-transformation)
  - [What are the main clustering patterns across our samples?](#what-are-the-main-clustering-patterns-across-our-samples)
  - [Is there a community-level signal of zoning?](#is-there-a-community-level-signal-of-zoning)
  - [Can a supervised approach separate the two zones?](#can-a-supervised-approach-separate-the-two-zones)
  - [Can we refine the clusters by selecting the most influential MAGs?](#can-we-refine-the-clusters-by-selecting-the-most-influential-mags)
  - [Are there stable indicators across GBR sectors?](#are-there-stable-indicators-across-gbr-sectors)
  - [Is the classification better than chance?](#is-the-classification-better-than-chance)
  - [Does the signal hold at the local (reef) scale?](#does-the-signal-hold-at-the-local-reef-scale)
- [Outputs](#outputs)
- [Notes on interpretation](#notes-on-interpretation)
- [Known issues](#known-issues)

---

## Study design

| Variable | Column | Role |
|---|---|---|
| Reef protection status | `FishingStatus` (`Closed`, `Open`) | Categorical outcome (Y); a reef-level attribute |
| GBR sector | `GBRSector` | Covariate; MINT study variable |
| Reef | `ReefName` | Unit of replication for protection status |
| Sampling trip | set by `trip_col` | Covariate; temporal (seasonal) effect |
| Season | `TripSeason` | Used only to filter incomplete samples |

Each reef has a single protection status (checked in 6.1). Samples from the same reef are therefore not independent replicates of zoning, so the tests that matter for zoning are run at the reef level, or with whole reefs or sectors left out of training.

---

## Requirements

R (4.3 or later) with:

| Purpose | Packages |
|---|---|
| Data handling | `phyloseq`, `dplyr`, `tibble`, `tidyr`, `forcats` |
| Community-level tests | `vegan` |
| Supervised models | `mixOmics`, `BiocParallel` |
| Plotting | `ggplot2`, `ggpubr`, `patchwork`, `VennDiagram`, `grid` |
| Parallel cross-validation | `doParallel`, `foreach`, `parallel` |
| Output | `openxlsx`, `knitr` |

The script expects these objects and helpers from earlier scripts in the project:

| Object | Description |
|---|---|
| `trip_col` | Name of the sampling trip column |
| `imos_dir` | Project root |
| `save_outputs` | Logical; write outputs when `TRUE` |
| `fish_cols` | Named colours for `Closed` / `Open` |
| `trip_col_vec`, `trip_col_vec_seasons` | Named colours for the four sampling trips |
| `param` | `BiocParallel` parameters for `perf()` |
| `save_fig()` | Saves a ggplot to `<dir>/<filename>` |
| `plot_fishstat_pca()` | Draws the CLR PCA from a phyloseq object |

---

## Inputs

One phyloseq object per habitat:

```
<OUT_dir>/2_MAG_ANALYSIS/mag_prok_<habitat>_MAG.rds
```

- `otu_table`: CoverM read counts per MAG (prokaryotes only).
- `tax_table`: GTDB ranks (`Class`, `Order`, `Family`, `Genus`, `Species`).
- `sample_data`: the design variables listed above.

---

## Configuration and running

Edit section 6.0:

```r
habitat  <- "SEA"      # or "SED"
OUT_dir  <- "~/Library/CloudStorage/OneDrive-TheUniversityofQueensland/5_mixomics"

run_slow <- TRUE       # tuning, permutation nulls, leave-one-reef-out
n_perm   <- 9999       # PERMANOVA / dbRDA permutations
n_cores  <- max(1, parallel::detectCores() - 2)
logging  <- TRUE       # copy console output to a log file
```

Then run:

```r
source("fishing_status_MAG.R", echo = TRUE)
```

With `logging = TRUE`, printed output is written to `<OUT_dir>/2_MAG_ANALYSIS/log_fishing_status_<habitat>_MAG_<YYYYMMDD_HHMM>.txt`. `echo = TRUE` makes bare object names print, so they are captured too. If the script stops partway, restore the console with `while (sink.number() > 0) sink(); sink(type = "message")`.

Tuning, the permutation test and the leave-one-reef-out validation are slow; a full run takes hours.

---

## Analysis workflow

### Data preparation and CLR transformation

Microbial abundance data are compositional: each sample's counts only carry information relative to one another. Following the mixOmics recommendations for microbial data ([mixMC preprocessing](http://mixomics.org/mixmc/mixmc-preprocessing/)), counts are centred log-ratio (CLR) transformed. The geometric mean cannot be computed with zero counts, so a pseudocount of 1 is added first (`clr_counts()`).

Section 6.1 also:

- keeps samples with complete design variables and drops MAGs with zero counts afterwards;
- checks that each reef has a single protection status and builds a reef lookup table;
- computes `zone_clr_diff`, the mean CLR in NTMRs minus fished reefs for each MAG, which gives the direction of enrichment (positive means NTMR-enriched). This is used instead of (s)PLS loading signs, which are arbitrary;
- defines `shuffle_zone_by_reef()`, which shuffles protection status among reefs within each sector and maps it back to samples. Every permutation null in the script uses it.

A per-sector CLR (`fish_clr_sec`) is also kept for the MINT inputs. Because CLR is computed within each sample, it is numerically identical to `fish_clr`.

### What are the main clustering patterns across our samples?

**6.2 PCA.** Principal Components Analysis is used as an unsupervised view of how reef sites cluster by community profile, before any zoning information is used. Samples are coloured by protection status to check whether zoning is visible on its own, or whether sites cluster mainly by geography (sector) and time (sampling trip).

### Is there a community-level signal of zoning?

PCA is only a visualisation. Community-level tests let us ask explicitly whether protection status explains a significant proportion of variance once sampling trip, sector and reef are included as covariates. Establishing this signal first gives a stronger foundation for the indicator search with supervised methods.

**6.3 PERMANOVA** (`vegan::adonis2`, Aitchison distance, that is Euclidean distance on CLR):

- protection status on its own;
- a sequential model of status, trip, sector and reef (`by = "terms"`), which shows the hierarchy of variance explained;
- a **reef-level model** on reef-mean CLR profiles, with sector entered before status and permutations restricted within sector. Because reefs, not samples, are the unit of replication for zoning, this is the inferential test for protection status.

**6.4 dbRDA.** Distance-based redundancy analysis complements PERMANOVA with a constrained ordination, showing how samples separate along axes that are linear combinations of the design variables. Outputs include a term-wise permutation ANOVA, the ordination by protection status and the share of constrained variance on each axis. Variance partitioning (`varpart`) separates the variance unique to zoning, trip and sector from the variance they share.

### Can a supervised approach separate the two zones?

**6.5 PLS-DA.** PLS-DA uses the class of each sample and finds components that maximise the covariance between the MAG profiles (X) and protection status (Y), rather than the variance of X as in PCA. Better separation than PCA is therefore expected and is not evidence on its own.

- **Choosing the number of components.** `perf()` evaluates how well the model classifies samples it was not trained on, using repeated cross-validation (4-fold, 50 repeats) on a 10-component model. With only two outcome classes, one or two components should be enough; if the error keeps dropping with more components, the later ones are likely capturing trip or sector structure, so the final model is kept small to avoid overfitting.
- **Final model.** Sample plots are drawn by protection status, by sampling trip and by trip season, to see whether zoning separates samples and whether trip effects remain.

### Can we refine the clusters by selecting the most influential MAGs?

**6.6 sPLS-DA.** Many MAGs are likely noisy or uninformative for discriminating zones. Sparse PLS-DA selects a small subset of MAGs on each component.

- **Tuning.** `tune.splsda()` estimates the balanced error rate (BER) over a grid of `keepX` values (1 to 150 MAGs per component), one component at a time, with each component conditional on those already tuned. A one-sided t-test decides whether adding a component improves classification.
- **Final model and performance.** `perf()` reports the overall error rate and the BER, which is the relevant metric because the number of samples per class is unbalanced.
- **Stability of the selection.** During repeated cross-validation, `perf()` records how often each MAG is selected. This answers the question: *how reproducible is the microbial signature when the training set is perturbed?* `selectVar()` output is joined with these stability scores.
- **Visualisation.** Sample plots, the loading plot for component 1 (bars coloured by the class in which each MAG is maximally enriched) and Clustered Image Maps (Euclidean distance, complete linkage) with reef sites coloured by sampling trip or by protection status.

If sPLS-DA mainly selects season-specific taxa and trip effects remain in the sample plots, the signature is not yet a stable indicator of zoning. For monitoring, the aim is indicators that hold across the GBR regardless of season and geography.

### Are there stable indicators across GBR sectors?

**6.6 MINT sPLS-DA.** Multivariate INTegration (MINT; Rohart et al. 2017) treats each GBR sector as a separate study, removing sector-specific variation while selecting MAGs that discriminate zones consistently across sectors. Sector is used as the study variable because AIMS LTMP benthic data are analysed by sector.

- **Number of components.** `perf()` with leave-one-group-out cross-validation (LOGOCV), so each sector is predicted by a model trained only on the others.
- **Trip versus sector as study.** The top component-1 MAGs selected with trip and with sector as the study variable are compared (Venn diagram) to check that the signature does not depend on that choice.
- **Error rate by sector.** Shows whether some sectors are harder to classify than others, overall and for each class.
- **Number of features.** `tune()` with LOGOCV over `keepX` values (10 to 500, BER, centroids distance), followed by the final model.
- **Visualisation.** Global and per-sector sample plots (does trip clustering disappear once sectors are accounted for?) and Clustered Image Maps of the selected MAGs.
- **AUROC.** Training ROC curves, and ROC curves from LOGOCV predictions, which reflect performance on unseen sectors.

### Is the classification better than chance?

**6.6c Permutation test.** Protection status is shuffled among reefs within sectors (999 permutations). For each permutation the MINT sPLS-DA model is refitted with the tuned `ncomp` and `keepX`, and the LOGOCV balanced accuracy is recorded. The observed accuracy is compared with this null distribution to give a permutation p-value.

Labels are shuffled among reefs, not samples, because zoning is fixed per reef; shuffling samples would create label patterns the design cannot produce and make the null too easy to beat. `keepX` and `ncomp` were tuned on the observed labels and held fixed for the permutations, so the null does not include tuning optimism.

### Does the signal hold at the local (reef) scale?

**6.6c Leave-one-reef-out (LORO) validation.** To test reproducibility at the local scale, sPLS-DA is trained on all reefs but one and used to predict the held-out reef, repeated for every reef (`run_loro()`). Predictions are pooled across reefs and summarised as a balanced accuracy and a confusion matrix, benchmarked against a null from the same procedure on reef-shuffled labels.

Plain sPLS-DA is used here rather than MINT, because the held-out reef's sector is still present in the training data, which does not fit MINT's study-wise design.

---

## Outputs

`<h>` is the habitat tag (`SEA` or `SED`).

### `2_MAG_ANALYSIS/0_PCA_PERMANOVA_RDA/`

| File | Content |
|---|---|
| `01_PCA_<h>fish.png` | CLR PCA |
| `02_PERMANOVA_full_<h>_<n_perm>.xlsx` | Sequential PERMANOVA table |
| `02_PERMANOVA_<h>_<n_perm>.png` | Variance explained per term (sample level) |
| `02_PERMANOVA_reeflevel_<h>_<n_perm>.png` | Variance explained per term (reef level) |
| `03_dbRDA_ANOVA_<h>.xlsx` | dbRDA term-wise ANOVA |
| `03_dbRDA_<h>ordination.png` | dbRDA axes 1 and 2 by protection status |
| `03_dbRDA_<h>axis.constrained.variation.png` | Constrained variance per axis |

### `2_MAG_ANALYSIS/1_PLSDA/`

| File | Content |
|---|---|
| `04_PLS_DA_<h>.error.rate.png` | PLS-DA error rate by component |
| `04_PLS_DA_<h>.initial.png` | 10-component model, components 1 and 2 |
| `04_PLS_DA_<h>.finalmodel.2components.by.{status,trip,trip_season}.png` | Final PLS-DA sample plots |
| `05_sPLS_DA_<h>.error.rate.png` | sPLS-DA `keepX` tuning |
| `05_sPLS_DA_<h>.feature.stability.png` | Selection frequency across CV folds |
| `05_sPLS_DA_<h>.finalmodel.by.{status,trip,tripseason}.png` | Final sPLS-DA sample plots |
| `05_sPLS_DA_<h>.finalmodel.contribution.comp1.png` | Loadings, component 1 |
| `05_sPLS_DA_cim_<h>{bytrip,bytripseason,byfishstatus}.png` | Clustered Image Maps |
| `06_mint_sPLS_DA_cim_<h>error.rate.trip.png` | MINT component selection |
| `06_mint_sPLS_DA_<h>tuning.png` | MINT `keepX` tuning |
| `06_mint_sPLS_DA_<h>.finalmodel{.global,.bysector,}.png` | MINT sample plots |
| `06_sPLS_DA_cim_<h>{bytripseason,byfishstatus}.png` | MINT Clustered Image Maps |
| `06_sPLS_DA_<h>.training.auroc.png` | MINT training ROC |
| `06_sPLS_DA_<h>.LOGOCV.auroc.png` | MINT LOGOCV ROC |
| `06_sPLS_DA_<h>permutation.png` | Permutation null and observed accuracy |
| `06_sPLS_DA_permutation_test_results.rds` | Permutation results |
| `06_sPLS_DA_Supplementary_Table_All_Permutations.xlsx` | All permuted accuracies |
| `06_sPLS_DA_<h>_leave_one_reef_out_validation_splsda.png` | LORO summary figure |
| `06_sPLS_DA_loro_validation_results_splsda.{rds,xlsx}` | LORO results |

### Cached models (`OUT_dir`)

`perf_plsda_<h>_MAG.rds`, `plsda_final_<h>_MAG.rds`, `perf_plsda_final_<h>_MAG.rds` and the MINT equivalents. Delete them to force a refit (the filenames do not include `ncomp`).

---

## Notes on interpretation

- **Zeros and the pseudocount.** The smallest nonzero count is in the thousands, consistent with CoverM's covered-fraction filter, so a zero means "below detection breadth" rather than true absence. A pseudocount of 1 sets a large gap between zero and detected MAGs; check that key results hold when zeros are replaced with half the smallest nonzero count.
- **Sample-level cross-validation is optimistic.** In the `Mfold` error rates for PLS-DA and sPLS-DA, samples from the same reef can appear in both training and test folds. Report the reef- and sector-grouped estimates (LORO, MINT LOGOCV, permutation test) as the measures of performance.
- **Sequential PERMANOVA and dbRDA.** Protection status enters first and takes credit for variance it shares with trip and sector. Reef is nested within status and sector, so those terms are aliased once reef is included. Variance partitioning shows the unique and shared fractions.
- **Explained variance in PLS-based plots.** Axis percentages matter less than in PCA, because PLS-DA maximises covariance with Y rather than variance in X.
- **Direction of effects.** Use `zone_clr_diff` (positive means NTMR-enriched), not loading signs.

---

## References

- Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). mixOmics: an R package for 'omics feature selection and multiple data integration. *PLoS Computational Biology* 13(11): e1005752.
- Rohart F, Eslami A, Matigian N, Bougeard S, Lê Cao K-A (2017). MINT: a multivariate integrative method to identify reproducible molecular signatures across independent experiments and platforms. *BMC Bioinformatics* 18: 128.
- Oksanen J et al. vegan: Community Ecology Package. R package.
