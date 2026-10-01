# Extended Data Figures Pipeline

This R Markdown document generates **extended data figures** for the variant classification manuscript. It produces Sankey diagrams, confusion matrices, heatmaps, and multi-ring donut plots to visualize classification performance across different variant categories and predictor combinations.

## Overview

The pipeline creates publication-quality figures that:
1. Compare ClinVar/ClinGen truth labels against predicted classifications
2. Visualize the flow of variants from original annotations to final classifications
3. Show reclassification rates for VUS
4. Display evidence source breakdowns (functional vs. predictive vs. both)
5. Highlight conflicting evidence patterns

## Requirements

### Software
- R 4.1.3
- RStudio (recommended for R Markdown)

### R Dependencies
```r
install.packages(c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "tidyverse",
  "readxl",
  "ggforce"
))

# ggsankey must be installed from GitHub
devtools::install_github("davidsjoberg/ggsankey")
```

## Input Data

### Required Files

| File | Location | Description |
|------|----------|-------------|
| `Supplementary_Data_5.xlsx` | `data/Supplemental_Data/` | Gene-specific calibration results (ExCALIBR + predictors) |
| `Supplementary_Data_6.xlsx` | `data/Supplemental_Data/` | OddsPath calibration results |

## Methods

### Figure Types

#### 1. Sankey Diagrams

Flow diagrams showing how variants move from their original classification (ClinVar, ClinGen, VUS, gnomAD) to their predicted classification.

**Functions:**
- `make_sankey_controls()` — ClinVar controls
- `make_sankey_clingen()` — ClinGen Evidence Repository variants
- `make_sankey_VUS()` — Variants of Uncertain Significance
- `make_sankey_gnomad()` — gnomAD population variants
- `make_sankey_preclass()` — Unobserved variants

#### 2. Confusion Matrices

Heatmap-style 2×2 matrices comparing true labels (ClinVar/ClinGen) against predicted classifications.

**Functions:**
- `make_confusion_matrix_controls_chek2()` — For ClinVar controls
- `make_confusion_matrix_clingen()` — For ClinGen variants

**Metrics calculated:**
- True Positives (TP), True Negatives (TN), False Positives (FP), False Negatives (FN)
- Sensitivity: TP / (TP + FN)
- Specificity: TN / (TN + FP)
- Concordance: (TP + TN) / Total
- Matthews Correlation Coefficient (MCC)
- Determinate call rate

**FP/FN breakdown statistics:**
- Functional evidence only
- Predictor evidence only
- Both evidence types
- Same vs. opposite direction evidence
- Direction of conflicting evidence

#### 3. VUS Reclassification Plots

Single-row confusion matrices showing how VUS are reclassified into determinate categories.

**Function:** `make_uncertain_confusion()`

#### 4. Per-Gene Analysis Plots

**Rank-based dot plot:**
Shows percentage of VUS remaining per gene, faceted by predictor (REVEL, AM, MP2).

**Heatmap:**
Shows distribution of total points across genes, faceted by predictor.

#### 5. Three-Ring Donut Plots

Nested donut charts with three concentric rings:

| Ring | Data | Description |
|------|------|-------------|
| Inner | Total Points bins | Classification outcome ranges |
| Middle | Point source | Functional only / Predictive only / Both |
| Outer | Conflict status | Conflicting / No conflict |

**Functions:**
- `build_three_rings()` — Computes ring geometry
- `make_three_ring_plot()` — Creates labeled plot
- `make_three_ring_plot_no_labels()` — Creates plot without labels

## Usage

### 1. Set Data Directories

Edit the paths at the top of the R Markdown file:

```r
DATA_DIR = "~/path/to/data"
OUT_DIR = "~/path/to/outputs"
```

### 2. Knit the Document

In RStudio:
```r
rmarkdown::render("Extended_data_figures.Rmd")
```

Or click "Knit" in the RStudio interface.

### 3. Output Files

Figures are saved to subdirectories within `OUT_DIR`:

## Extended Data Figure Mapping

| Figure | Description | Key Visualizations |
|--------|-------------|-------------------|
| **ED Fig 6alt** | Controls/Clingen performance (AM, MP2) | Sankey diagrams, confusion matrices |
| **ED Fig 8alt** | Controls/Clingen performance OddsPath | Sankey diagrams, confusion matrices |
| **ED Fig 9alt** | VUS reclassification (AM, MP2) | Sankey diagrams, reclassification bars |
| **ED Fig 12** | Per-gene VUS analysis | Rank plot, heatmap |
| **ED Fig 11** | Evidence source breakdown | Three-ring donut plots for gnomAD, VUS, unobserved |

## Reviewer-only gene-specific variant (ED Fig 8alt/8)

`Extended_data_figures_gene_specific.Rmd` rebuilds ED Fig 8alt and ED Fig 8
-- the only two figures here derived from `Supplementary_Data_6.xlsx` --
from `Supplementary_Data_6_gene_specific.xlsx` instead (the variant
`OddsPath_classifications.ipynb` produces with
`ONLY_ODDSPATH_CALIBRATED_DATASETS = True` and
`USE_GENE_SPECIFIC_PREDICTOR_CALIBRATIONS = True`; see
`notebooks/analysis/README_OddsPath_classifications.md`'s "Reviewer-only
variant" section). Not used in the paper. Render it the same way as the
main `.Rmd`:

```r
rmarkdown::render("Extended_data_figures_gene_specific.Rmd")
```

Writes to `extended_data_figure_8alt_gene_specific/` and
`extended_data_figure_8_gene_specific/` under `OUT_DIR` -- distinct
folders so nothing here can collide with or overwrite the standard ED
Fig 8alt/8 outputs. See [`docs/figures.md`](../../../docs/figures.md#extended-data-figure-8alt8-gene-specific-reviewer-only)
for the full write-up, including exact output filenames.

## Figure Dimensions

All figures are saved at 300 DPI with transparent backgrounds.

## Notes

- ClinVar significance is collapsed: P/LP → Pathogenic, B/LB → Benign
- The pipeline uses `ggsankey` for Sankey diagrams and `ggforce` for arc/donut plots
- Some plots require specific column types (text vs. numeric) — the script handles conversions

## References

- ggsankey package: https://github.com/davidsjoberg/ggsankey
- ggforce package: https://ggforce.data-imaginist.com/
