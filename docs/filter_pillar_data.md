# Filter Pillar Data

`src/filter_pillar_data.py` rebuilds the VUS / gnomAD / Unobserved category
datasets used by Extended Data Figure 10, replicating
`notebooks/analysis/Variant_Classification_analysis.ipynb`'s filtering,
conflict-tagging, category-selection, and REVEL-classification steps in the
same order. It is the first step of the Extended Data Figure 10 data-prep
workflow -- see `docs/gnomad_ancestry_resolution.md` for the full
walkthrough from obtaining gnomAD data through generating the figure.

## What it does

In order:

1. Computes functional points for every row (`excalibr_points`, with the
   Supp 4 ExCALIBR intervals filling in rows ExCALIBR left unscored, and
   OddsPath points for the F9/TP53 meta-analysis datasets); missing points
   default to 0.
2. Restricts F9/TP53 to the meta-analysis datasets (drops the individual
   assays that feed into that meta-analysis).
3. Tags (without yet removing) splice variants the assay doesn't measure,
   and start-lost variants in amino-acid-level assays.
4. Tags variants whose functional points conflict in sign between assays,
   separately within nucleotide-level and within amino-acid-level assays.
5. Spreads splice and REVEL-training flags to every nucleotide variant
   encoding the same amino-acid change.
6. Keeps only nucleotide- and amino-acid-level rows.
7. Tags sign conflicts between nucleotide- and amino-acid-level assays, and
   extends each tag to every row sharing that measurement
   (`mavedb_variant_urn`).
8. Removes SFPQ.
9. Applies the CHEK2 QC flag (redundant with upstream, included for parity
   with the notebook).
10. Removes tagged conflicts, unmeasured splice variants, and start-lost
    variants.
11. Removes rows with `Flag == '*'`.
12. Selects the requested category (VUS / gnomAD / Unobserved), keeps one
    row per variant (greatest absolute functional points, `Dataset` name as
    tie-break), and drops REVEL training variants into a separate output.
13. Classifies each variant: functional points + REVEL points
    (gene-specific threshold, genome-wide fallback, else 0) ->
    `Class_REVEL` (Pathogenic >=10, Likely Pathogenic 6-9, Uncertain 0-5,
    Likely Benign -6..-1, Benign <=-7); flags functional-vs-REVEL
    disagreement in `Conflicting_REVEL`.

## Inputs

- `--points-input` (default `data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz`)
  -- the ExCALIBR calibration run's output (see `docs/add_excalibr_columns.md`
  for its provenance).
- `--supp4-input` -- `Supplementary_Data_4.xlsx` (typically
  `data/output/supplementary_data/Supplementary_Data_4.xlsx`), used for the
  ExCALIBR interval backfill, F9/TP53 OddsPath points, and the final REVEL
  classification. Omit it to skip all three (statistics-only / ablation
  mode); the script logs when it does.
- `--supp3-input` -- `Supplementary_Data_3.xlsx` (typically
  `data/input/maves/Supplementary_Data_3.xlsx`), used for interval-based
  F9/TP53 functional classes.
- `--chek2-input` -- `CHEK2_Gebbia_2024.xlsx` (typically
  `data/input/maves/CHEK2_Gebbia_2024.xlsx`), optional; the CHEK2 QC flag it
  sets is already applied upstream, so this is redundant unless rerunning
  with different upstream data.

## Outputs

Written to `--outdir` (default `data/intermediate/figures/extended_data_figure_10/`,
a gitignored staging directory -- these feed `src/map_gnomad_ancestry.py`
next, not the manuscript directly):

- `replicated_<category>_category.csv` -- the deduplicated category
  selection.
- `replicated_<category>_category_REVEL.csv` -- the same, with REVEL
  training variants removed. This is the file to pass to
  `src/map_gnomad_ancestry.py --input`.

## Usage

```bash
poetry run python -m src.filter_pillar_data \
    --points-input data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz \
    --supp4-input data/output/supplementary_data/Supplementary_Data_4.xlsx \
    --supp3-input data/input/maves/Supplementary_Data_3.xlsx \
    --category VUS
```

Run once per category needed (`VUS`, `gnomAD`, or `Unobserved` -- Extended
Data Figure 10 needs `VUS` and `gnomAD`).

Ablation flags: `--no-excalibr-backfill` (leave rows ExCALIBR didn't score
at 0, the old behavior), `--no-urn-propagation` (cross-resolution conflict
flags only direct rows, not URN siblings), `--keep-bp4-very-strong` (apply a
`BP4_Very Strong` REVEL tier if Supp 4 still has one -- the current pipeline
drops it), `--validate-n` (rows sampled to check that Supp 4's intervals
reproduce `excalibr_points`; 0 skips the check -- a low agreement rate means
Supp 4 is stale).
