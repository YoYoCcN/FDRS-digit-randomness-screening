# FDRS: Fabrication-risk Digit Randomness Screening

A statistical and machine-learning-assisted framework for screening non-random decimal digit-pattern irregularities in raw numerical research data.

**FDRS is an auxiliary screening and prioritization tool. It is not a standalone detector of data fabrication, falsification, or research misconduct.** Positive findings should be interpreted together with original instrument files, laboratory records, data-processing workflows, rounding rules, repeated experiments, and expert review.

## Overview

FDRS analyzes single-column numerical datasets and extracts decimal digit features to evaluate whether the fine-scale digit structure is consistent with expected or reference digit-randomness patterns.

The current scripts implement:

- single-decimal-digit multinomial testing;
- joint two-decimal-digit distribution analysis;
- chi-square statistics, P values, Cramer's V, Shannon entropy, normalized entropy, KL divergence, and standardized residuals;
- digit-preference indices;
- progressive subsampling stability analysis;
- semi-supervised machine-learning risk scoring using Random Forest, Elastic-net Logistic Regression, SVM radial, Isolation Forest, and ensemble risk integration;
- publication-style visualization using a blue/pink/red color palette.

## Repository structure

```text
FDRS_GitHub_repository/
├── R/
│   ├── 01_fdrs_full_statistics_progressive_visualization.R
│   ├── 02_fdrs_ml_realrawdata3_decimal_1_2.R
│   └── 03_fdrs_ml_realrawdata2_decimal_3_4.R
├── data/
│   └── README.md
├── examples/
│   └── example_single_column_input.txt
├── results/
│   └── .gitkeep
├── docs/
│   └── methodology_note.md
├── install_dependencies.R
├── DESCRIPTION
├── CITATION.cff
├── LICENSE
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
└── .gitignore
```

## Scripts

### `R/01_fdrs_full_statistics_progressive_visualization.R`

Full statistical screening and visualization for two target datasets.

Default input files in the script:

```r
INPUT_FILES <- c("RawData.txt", "ErrData.txt")
INPUT_NAMES <- c("RawData", "ErrData")
SINGLE_DIGIT_POS <- 3
JOINT_START_POS  <- 3
JOINT_LENGTH     <- 2
```

This script outputs full-sample digit statistics, progressive subsampling summaries, and visualizations.

### `R/02_fdrs_ml_realrawdata3_decimal_1_2.R`

Semi-supervised machine-learning FDRS pipeline for a target dataset analyzed using the first and second decimal positions.

Default target:

```r
TARGET_FILES <- c("RealRawData3.txt")
TARGET_NAMES <- c("RealRawData3")
PRIMARY_SINGLE_POS <- 1
JOINT_START_POS <- 1
JOINT_LENGTH <- 2
```

### `R/03_fdrs_ml_realrawdata2_decimal_3_4.R`

Semi-supervised machine-learning FDRS pipeline for a target dataset analyzed using the third and fourth decimal positions.

Default target:

```r
TARGET_FILES <- c("RealRawData2.txt")
TARGET_NAMES <- c("RealRawData2")
PRIMARY_SINGLE_POS <- 3
JOINT_START_POS <- 3
JOINT_LENGTH <- 2
```

## Input format

Input files should be plain-text, single-column numerical files, with one independent numerical value per row.

Example:

```text
0.3167
0.4281
0.5029
0.3915
```

The scripts automatically remove blank lines and non-numeric entries. Values in scientific notation are converted to fixed-decimal format before decimal digit extraction.

## Quick start

1. Install R packages:

```r
source("install_dependencies.R")
```

2. Put your input `.txt` files into the working directory used by the scripts.

3. Open the target script and edit the working directory:

```r
setwd("YOUR/LOCAL/PROJECT/PATH")
```

4. Edit the input file names and target decimal positions if needed.

5. Run the script in RStudio or from terminal:

```bash
Rscript R/01_fdrs_full_statistics_progressive_visualization.R
```

## Main outputs

Depending on the script, outputs may include:

- full-statistics summary CSV files;
- digit-frequency tables;
- joint digit-combination tables;
- progressive subsampling raw and summary tables;
- machine-learning performance metrics;
- model-specific and ensemble risk scores;
- risk grades;
- PCA, ROC, calibration, feature-importance, risk-score, and progressive subsampling visualizations.

## Interpretation boundary

FDRS detects numerical digit-structure irregularities. Such irregularities may arise from multiple sources, including:

- manual fabrication or selective modification;
- copy-and-paste reuse;
- rounding or formatting rules;
- instrument precision limits;
- bounded measurement ranges;
- data transformation or batch processing;
- repeated-measurement structures;
- benign recording habits.

Therefore, FDRS output should be used to prioritize further review, not to adjudicate misconduct.

## Recommended citation

Please cite this repository and the corresponding manuscript if you use FDRS. See `CITATION.cff`.

## License

This repository is distributed under the MIT License. See `LICENSE`.
