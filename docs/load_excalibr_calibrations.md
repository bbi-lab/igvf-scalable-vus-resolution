# Load exCALIBR Calibrations

`src/load_excalibr_calibrations.py` loads exCALIBR MAVE point-value
calibrations into Supplementary Data 4, so the workbook stays in sync with
the latest calibration run without hand-editing spreadsheet rows.

It reads every `*.json` calibration file in a directory (default:
`data/input/mave_calibration/excalibr/json/`)
and writes one row per file into the `ExCALIBR_calibrations` sheet of a
workbook (default: `data/output/supplementary_data/Supplementary_Data_4.xlsx`),
**overwriting that sheet's existing contents in place**. All other sheets,
the sheet's tab position, and its column widths are left untouched.

If `Supplementary_Data_4.xlsx` doesn't exist yet in a fresh checkout, start
from a copy of `data/input/mave_calibration/Supplementary_Data_4.seed.xlsx`
-- see [Where the gene-specific calibration sheets come from](#where-the-gene-specific-calibration-sheets-come-from)
below.

## Where the calibration JSON files come from

The per-dataset JSON files under `data/input/mave_calibration/excalibr/json/`
are produced by an external ExCALIBR calibration run that itself consumes
this repo's own `data/output/maves/integrated_variant_effect_dataset.tsv.gz`
as input -- so a fresh set is only obtainable after the variant-annotation
pipeline has already completed (see `docs/variant_annotation_pipeline.md`).
The same external run also produces
`data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz`
(see `docs/add_excalibr_columns.md`) and
`data/input/mave_calibration/excalibr/fig4_data_oob.tar.gz` (see
`docs/build_figure4_data.md`). A copy of the JSON files is checked into this
repo so the pipeline is reproducible without rerunning ExCALIBR; refresh
them by replacing this directory's contents the next time ExCALIBR
calibration is rerun.

## Where the gene-specific calibration sheets come from

`Supplementary_Data_4.xlsx` also carries `REVEL_gene_specific_calibration`,
`MP2_gene_specific_calibrations`, and `AM_gene_specific_calibrations` --
pre-populated, gene-specific predictor-score calibration thresholds that
**no script in this repo produces or overwrites**. Unlike
`ExCALIBR_calibrations`/`OddsPath_calibrations` above, these come from an
external source and must already be present in the workbook before this
script (or `load_oddspath_calibrations`) runs; this script never touches
them. `data/input/mave_calibration/Supplementary_Data_4.seed.xlsx` is a
checked-in snapshot of a workbook with these three sheets populated -- copy
it to `data/output/supplementary_data/Supplementary_Data_4.xlsx` as your
starting point in a fresh checkout (its own `ExCALIBR_calibrations`/
`OddsPath_calibrations` sheets are stale snapshots too, but this script and
`load_oddspath_calibrations` overwrite them in place on the next run
regardless).

## Column mapping

| Column | Source |
|---|---|
| `dataset` | The JSON filename, minus `.json` |
| `prior` | `prior` (float, copied through as-is) |
| `range_-8` .. `range_-1`, `range_1` .. `range_8` | `point_ranges` (see below) |
| `relax` | `relax` (0/1 -> Excel `FALSE`/`TRUE`) |
| `n_c` | `n_c` (string, copied through as-is) |
| `benign_method` | `benign_method` (string, copied through as-is) |
| `clinvar_2018` | `clinvar_2018` (0/1 -> Excel `FALSE`/`TRUE`) |
| `scoreset_flipped` | `scoreset_flipped` (0/1 -> Excel `FALSE`/`TRUE`) |

`point_ranges` is a JSON object keyed by point value (`"-8"` .. `"-1"`, `"1"`
.. `"8"`; there is no `"0"` key) whose value is an array of `[low, high]`
pairs -- usually at most one, but a point can have more than one disjoint
range (e.g. `DDX3X_Radford_2023`, whose +1 and -1 points each cover two
separate score intervals). A point with an empty array is written as a blank
cell; a point with one pair is written as `"<low> <high>"`, with `-Infinity`
and `Infinity` rendered as `-inf` and `inf`; a point with more than one pair
has each one formatted the same way and joined with `", "`, sorted by lower
endpoint (e.g. `"-0.56 -0.52, 0.46 inf"`) -- the sheet still has only one
column per point, so multiple ranges share a cell rather than a column.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.load_excalibr_calibrations \
  data/input/mave_calibration/excalibr/json \
  data/output/supplementary_data/Supplementary_Data_4.xlsx
```

Both positional arguments default to the paths above, so a bare invocation
works from the repo root.

Via Docker (same image as `flag_variants`/`mave_dataset_stats`, see
`compose.yaml`):

```bash
src/scripts/run_load_excalibr_calibrations.sh [json-dir] [workbook]
```

Like `run_mave_dataset_stats.sh`, this wrapper doesn't map paths against a
`/work` staging mount -- the script only reads/writes `data/input/`/
`data/output/` files, which are already available at `/usr/src/app` via the
whole-repo bind mount.
