# Tongue-Microbiome

Longitudinal SHIP cohort omics data analysis using machine learning and multivariate statistical methods in R.

## Overview

This repository contains code for analysing tongue microbiome data from the **Study of Health in Pomerania (SHIP)** cohort. The workflow links tongue microbiome count data with periodontal outcome data from **SHIP-1** and **SHIP-2** to examine microbial patterns associated with periodontal status and progression over time.

The analysis includes:

- data import and merging of SHIP clinical and microbiome datasets
- genus-level missingness assessment and filtering
- centered log-ratio (CLR) transformation for compositional microbiome data
- exploratory analysis using principal component analysis (PCA)
- supervised modelling using PLS-DA and sPLS-DA
- longitudinal genus-difference analysis between SHIP-1 and SHIP-2
- adjusted cross-sectional and longitudinal regression modelling

## Repository structure

- **SHIP 1 & 2.R**  
  Preprocessing pipeline for SHIP-1 and SHIP-2, including data merging, missingness checks, CLR transformation, PCA, and PLS-DA.

- **sPLS-DA.R**  
  Tuning and fitting of sparse PLS-DA models for SHIP-1 and SHIP-2 using CLR-transformed genus data.

- **Difference in Genus.R**  
  Comparison of SHIP-1 and SHIP-2 genus profiles, overlap of selected genera, and calculation of longitudinal CLR differences.

- **longitudinal modelling.R**  
  Genus-wise longitudinal modelling of change in CLR abundance between SHIP-1 and SHIP-2 with adjustment for baseline covariates.

- **Cross-sectional modelling.R**  
  Cross-sectional genus-wise modelling of CLR-transformed abundances within SHIP-1 and SHIP-2.

## Data requirements

This repository contains **analysis code only**. The underlying SHIP datasets are not included.

To run the scripts, you will need access to:

- SHIP-1 and SHIP-2 clinical datasets
- periodontal grouping / outcome dataset
- taxonomy-level tongue microbiome count data

Because file paths in the scripts are local, they should be updated before execution.

## Methods summary

### 1. Data preprocessing
- Load SHIP-1 and SHIP-2 clinical data.
- Merge periodontal outcome information with microbiome taxonomy tables.
- Focus analyses on genus-level abundance data.
- Treat zero counts as missing for missingness filtering where applicable.
- Retain genera passing the missingness threshold.

### 2. Compositional transformation
- Reshape long-format genus counts into subject-by-genus matrices.
- Add pseudocounts where needed.
- Apply centered log-ratio (CLR) transformation before downstream analyses.

### 3. Exploratory and supervised analyses
- **PCA** for unsupervised exploration of variance structure.
- **PLS-DA** for supervised discrimination of periodontal groups.
- **sPLS-DA** for sparse feature selection and identification of genera contributing most strongly to class separation.

### 4. Longitudinal analysis
- Align shared subjects across SHIP-1 and SHIP-2.
- Compute longitudinal change on the CLR scale:

```r
Delta CLR = SHIP-2 CLR - SHIP-1 CLR
```

- Model genus-specific change using adjusted linear regression.

### 5. Covariate adjustment
Longitudinal models include baseline covariates such as:

- age
- sex
- smoking
- diabetes
- BMI
- HbA1c
- dental visit reason
- previous periodontal treatment
- powered toothbrush use
- interdental cleaning aid use

## Main R packages

The analysis uses packages including:

- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `haven`
- `mixOmics`
- `compositions`
- `plotly`
- `patchwork`
- `broom`
- `stringr`
- `readr`
- `gridExtra`

Example installation:

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2", "haven",
  "mixOmics", "compositions", "plotly", "patchwork",
  "broom", "stringr", "readr", "gridExtra"
))
```

## Suggested run order

```r
source("SHIP 1 & 2.R")
source("sPLS-DA.R")
source("Difference in Genus.R")
source("longitudinal modelling.R")
source("Cross-sectional modelling.R")
```

## Outputs

Depending on the script, outputs may include:

- CLR-transformed genus matrices
- PCA plots
- PLS-DA and sPLS-DA score and loading plots
- selected genera from sPLS-DA
- longitudinal change summaries
- adjusted regression result tables with FDR correction
- forest plots and manuscript-ready summary tables

## Notes

- The scripts were developed for SHIP tongue microbiome analyses and may need adaptation for other datasets.
- Some scripts assume that objects created in earlier scripts are already present in the R environment.
- Local working directories and filenames should be updated before running the code.
- Interpretation should account for the compositional nature of microbiome count data and the observational design of the study.

## Author

**Stuti Sinha**
