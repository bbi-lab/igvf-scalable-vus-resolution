# Biobank Carrier Odds-Ratio Analysis

Figure 2i, Figure 6b, and Extended Data Figure 2 each plot carrier odds
ratios computed against a biobank cohort -- functionally classified or
predictor-scored variant carriers vs. non-carriers, per gene/dataset. That
analysis can't run inside this repo: it depends on individual-level biobank
genotype/phenotype data under access restrictions this codebase doesn't
(and shouldn't) have a path to. The workflow instead spans three
environments:

1. **This repo** builds the input: the variant list + classification
   columns the biobank analysis needs.
2. **An external biobank environment** -- currently the [All of Us
   Research Program](https://www.researchallofus.org/)'s **Researcher
   Workbench**, a access-restricted data-exploration environment -- runs
   the actual carrier/non-carrier odds-ratio computation. Nothing in this
   repo can run or automate this step.
3. **This repo** (via `src/lib/figures/load_biobank_or_table.R`) consumes
   the result and feeds it to the three figure scripts above.

## 1. Building the input (this repo)

`src/build_variant_reclassification_dataset.py`, run with no `--dedup` flag
(README.md's main pipeline, step 6 -- see
[`docs/build_variant_reclassification_dataset.md`](build_variant_reclassification_dataset.md)),
produces
`data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz`
-- one row per surviving variant-effect measurement, with ACMG/AMP
points/classification columns. This is the file handed off to whoever runs
the biobank-side analysis.

## 2. The external analysis (All of Us Researcher Workbench)

There is currently no scripted upload step in this repo (contrast with
`src/submit_supplementary_data_1_to_igvf.py`, which does script a portal
upload elsewhere) -- the hand-off to the Researcher Workbench is manual:
someone with Researcher Workbench access uploads the file from step 1,
joins it against the All of Us cohort, and computes carrier log-odds-ratios
(with 95% CI) per Dataset/Predictor × Grouping/Calibration × Group/Range
bin, with carrier/control counts privacy-censored at ≤ 20. The result is
downloaded back out as a TSV (`AoU-OR-estimates_<date>.tsv.gz`) in the
Workbench's own column schema (`Dataset/Predictor`, `Grouping/Calibration`,
`Group/Range`, `Consequence`, `LogOR`, `LogOR_LI`, `LogOR_UI`, `Carrier
cases`, `Carrier controls`, `Total cases`, `Total controls`).

This repo cannot reproduce this step -- it's included here only so the full
chain is documented, not as something to automate or re-run casually.

## 3. Consuming the result (this repo)

The downloaded file is checked into `data/input/biobank/` -- currently
`AoU-OR-estimates_2026-10-07.tsv.gz`. `src/lib/figures/load_biobank_or_table.R`
translates that schema back into the column names/vocabulary
`Figure_2i.R`/`Figure_6b.R`/`Extended_Data_Figure_2.R` were originally
written against (`Dataset`, `Gene`, `Classifier`, `Classification`, `LogOR`,
`LogOR_LI`, `LogOR_UI`, `Cases with variants`), so those three scripts only
need a `load_biobank_or_table(path)` call in place of a plain `read_csv()`.
See that file's own comments for the exact vocabulary mapping and how
privacy-censored carrier counts (`"<= 20"`) are handled.

## Why this is a stand-in, not the long-term source

This same carrier odds-ratio analysis also has a stable, citable home:
**`IGVFFI3804AVJR.csv.gz`**, an IGVF data-portal export (intended download:
<https://data.igvf.org/tabular-files/IGVFFI3804AVJR/>) in a different
column layout, which the three figure scripts were originally written to
read directly. That export is currently unavailable pending a re-export (no
working download URL as of this writing); once one exists, it'll be
documented here, and the three figure scripts' `load_biobank_or_table()`
calls should be swapped back to a plain `read_csv()` against it (each
script's header comment says so). IGVF may publish its own re-export as a
separate, independently citable Zenodo record once it exists -- that would
be a different deposit from this project's own.

Until then, the All of Us file above is the only available copy of this
analysis and is checked into the repo rather than treated as something each
user must separately obtain. It's also the file this project distributes in
its own Zenodo data-prerequisites archive, since it's the one actually used
to produce the published figures.
