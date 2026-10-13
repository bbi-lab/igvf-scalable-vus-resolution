# Add ExCALIBR Columns

`src/add_excalibr_columns.py` joins the ExCALIBR point/evidence columns
produced by the ExCALIBR calibration run onto the integrated variant effect
dataset.

This is a **manual step**, run after the main variant-annotation pipeline
(see `docs/variant_annotation_pipeline.md`) has produced
`data/output/maves/integrated_variant_effect_dataset.tsv.gz` -- the ExCALIBR
output this script reads isn't available until later in the process, so
this can't be folded into that pipeline run.

## Inputs

- `data/output/maves/integrated_variant_effect_dataset.tsv.gz` -- the
  pipeline's integrated dataset (95 columns).
- `data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz`
  -- the same dataset as produced by the ExCALIBR calibration run, with
  seven additional columns: `excalibr_clinvar_release`, `excalibr_points`,
  `excalibr_acmg_evidence_code`, `excalibr_lr_plus`, `excalibr_prior`,
  `excalibr_posterior`, `excalibr_filter_reason`.

  **Where this comes from:** an external ExCALIBR calibration run that
  itself consumes this repo's own
  `data/output/maves/integrated_variant_effect_dataset.tsv.gz` as input --
  so it's only obtainable after the variant-annotation pipeline above has
  already completed, same as `data/input/mave_calibration/excalibr/json/`
  (see `docs/load_excalibr_calibrations.md`) and
  `data/input/mave_calibration/excalibr/fig4_data_oob.tar.gz` (see
  `docs/build_figure4_data.md`), both produced by the same external run. A
  copy is checked into this repo at the path above (note the `.csv.gz` --
  distinct from this script's own `.tsv.gz` *output*, which shares the same
  stem) so the pipeline is reproducible without rerunning ExCALIBR; refresh
  it by replacing that file the next time ExCALIBR calibration is rerun on
  a newer integrated dataset.

The two files' other 95 columns carry the same values (modulo cosmetic
formatting differences -- e.g. `-1.0` vs `-1`, `False` vs `false`, trailing
zeros on floats -- introduced by the ExCALIBR run's own read/write path);
this script trusts the integrated dataset's own copies of those and only
pulls in the seven new `excalibr_*` columns.

## Join key

(`Dataset`, `mavedb_variant_urn`, `hgvs_c`) -- the only combination of
columns shared between the two files that is unique in both (neither
`Dataset` + `mavedb_variant_urn` alone, nor a genomic-coordinate-based key,
is unique: some datasets have several MaveDB score sets or several
transcript annotations sharing the same variant urn/coordinates). The
script raises `ValueError` (surfaced as a `ClickException`) if either input
has duplicate keys, or if the two inputs' key sets don't match exactly.

## Output

- `data/output/maves/integrated_variant_effect_dataset_with_excalibr.tsv.gz`
  -- the integrated dataset with the seven `excalibr_*` columns appended.
- `data/output/supplementary_data/Supplementary_Data_1.tsv.gz` -- a copy of
  the same file, written every run unless `--skip-supplementary-data-copy`
  is passed.

## Usage

```bash
poetry run python -m src.add_excalibr_columns
```

All four paths (integrated input, ExCALIBR input, output, supplementary
data copy) default to the locations above, so a bare invocation works from
the repo root once both inputs exist. Override any of them with
`--integrated-input`, `--excalibr-input`, `--output`, or
`--supplementary-data-copy`.
