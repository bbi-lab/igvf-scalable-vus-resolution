# Build Variant Reclassification Dataset

`src/build_variant_reclassification_dataset.py` turns a
`Variant_Classification_analysis.ipynb` checkpoint into either the
biobank analysis input file (one row per surviving measurement) or a
deduplicated reclassification dataset (one row per DNA variant), by
re-applying the notebook's own downstream exclusions and adding six
ACMG/AMP points columns (`ExCALIBR_points`, `OddsPath_points`,
`Functional_points`, `<predictor>_points`, `Conflict_<predictor>_*`,
`Combined_points` -- `REVEL_points`/`Conflict_REVEL_GeneSpecific` by
default, see "Choosing a predictor" below) -- see the script's own module
docstring for the exact column definitions and exclusion rules.

In README.md's main-pipeline "Running it" steps, this script runs twice
against the same `Variant_Classification_analysis.ipynb` checkpoint: once
per predictor with `--dedup` (step 5, one deduplicated reclassification
dataset per `REVEL`/`AM`/`MP2`), then once more with no arguments (step 6,
the default-predictor, non-deduplicated biobank analysis input file).

## Choosing a checkpoint

The positional `checkpoint_file` argument selects which notebook run to
build from. Three are currently in use, differing only in how F9/TP53 are
scored (every other gene is identical across all three):

| Checkpoint | Produced by | F9/TP53 scoring |
|---|---|---|
| `integrated_variant_effect_dataset_analysis.csv.gz` (default) | `Variant_Classification_analysis.ipynb` | Limited to their meta-analysis/model score sets (`F9_Popp_2025_model`, `TP53_Fayer_2021_meta`), scored via OddsPath |
| `integrated_variant_effect_dataset_analysis_excalibr_all_genes.csv.gz` | `Variant_Classification_analysis_ExCALIBR_all_genes.ipynb` | Every gene, F9/TP53 included, scored via ExCALIBR; restores F9/TP53's per-assay datasets (the only ones with an ExCALIBR calibration on file) |
| `integrated_variant_effect_dataset_analysis_excalibr_assays_oddspath_models.csv.gz` | `Variant_Classification_analysis_ExCALIBR_assays_OddsPath_models.ipynb` | Individual assay datasets scored via ExCALIBR; the model/meta score sets keep OddsPath |

All three live under `data/output/reclassification/` and are gitignored
build outputs, not checked in. See `docs/tp53_f9_dataset_inclusion.md` for
the full F9/TP53 dataset-inclusion matrix across them.

A checkpoint must be generated (or refreshed, if its notebook changed)
before this script can read it:

```bash
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=550 \
  --output executed_<notebook-name>.ipynb \
  notebooks/analysis/<notebook-name>.ipynb
```

## Output modes

- **Biobank analysis input file** (default, no `--dedup`): one row per
  surviving dataset/measurement -- a variant scored by more than one
  dataset keeps one row per dataset. This is the form collaborators use
  for downstream biobank analysis. 866,359 rows on the official checkpoint.
- **Deduplicated reclassification dataset** (`--dedup`): collapses to one
  row per DNA variant via `src.lib.dedup.dedup_by_max_abs_points` (greatest
  `abs(Combined_points)` wins, ties broken by `Dataset` name), and
  additionally excludes bare `conflicting_fxn_data` rows outright --
  picking a "winner" among datasets that genuinely disagree on a variant's
  effect would silently misrepresent that disagreement as a clear result
  (see the module docstring). 635,028 rows on the official checkpoint.

Both modes work with any of the three checkpoints above. `--output`
defaults to `integrated_variant_effect_biobank_input_data.tsv.gz`
regardless of `--dedup` -- **always pass `--output` explicitly when using
`--dedup`**, or it overwrites the per-measurement file at that default
path with the deduplicated one.

## Choosing a predictor

`Combined_points` (`Functional_points` plus a computational predictor's
points) defaults to REVEL, gene-specific with genome-wide fallback --
`--predictor AM`/`--predictor MP2` swap in AlphaMissense or MutPred2
instead, and `--predictor-scope genome-wide` drops the gene-specific
contribution entirely. The output's predictor-named columns follow suit
(`<predictor>_points`, `Conflict_<predictor>_GeneSpecific` or
`Conflict_<predictor>_GenomeWide`), reducing to today's `REVEL_points`/
`Conflict_REVEL_GeneSpecific` for the default predictor/scope.

Two things change with the predictor, not just the column values:

- The training-circularity exclusion (dropping variants used to train the
  predictor, to avoid the file scoring a variant partly by how well it was
  predicted during that predictor's own training) switches from
  `revel_train_amino` to `mp2_train_amino` for `--predictor MP2`. `AM` has
  no such column in the checkpoint, so no row is excluded on this basis
  when `--predictor AM` is selected.
- `--dedup`'s winning measurement for a variant scored by more than one
  dataset can change, since it's chosen by greatest `abs(Combined_points)`
  and `Combined_points` is predictor-dependent. On the official checkpoint,
  switching `REVEL` -> `AM` changes the winning dataset for 1,634 of
  635,028 variants; `REVEL` -> `MP2` changes it for 10,117. Two real
  examples, both candidates sharing `Functional_points`:

  - **BRCA1 `chr17:43049179 A>G` (`p.Met1783Thr`)**, two datasets:
    `BRCA1_Findlay_2018` (`Functional_points=-5`) and
    `BRCA1_Adamovich_2022_Cisplatin_Resistance` (`Functional_points=0`).
    Both share `REVEL_points=-1` and `AM_points=3`. Under `REVEL`,
    `Combined_points` is -6 and -1 -- `Findlay_2018` wins (`abs`=6). Under
    `AM`, it's -2 and 3 -- `Adamovich_2022_Cisplatin_Resistance` wins
    instead (`abs`=3 > 2): the predictor swap alone flips which dataset's
    row survives `--dedup`, with no change to either dataset's own
    functional evidence.
  - **PAX6 `chr11:31801588 CCC>TAT` (`p.Gly124Ile`)**, two datasets:
    `PAX6_McDonnell_2024_LE9_geneticin` (`Functional_points=0`) and
    `PAX6_McDonnell_2024_LE9_no_geneticin` (`Functional_points=-8`). Both
    share `REVEL_points=0` and `MP2_points=4`. Under `REVEL`,
    `Combined_points` is 0 and -8 -- `no_geneticin` wins outright
    (`abs`=8). Under `MP2`, it's 4 and -4 -- an exact tie in `abs`, broken
    by `Dataset` name ascending (`"...geneticin"` sorts before
    `"...no_geneticin"`), so `geneticin` wins instead. This one shows the
    tie-break rule actually deciding a real case, not just the predictor
    switch.

## Usage

```bash
# Official deduplicated reclassification dataset, one per predictor (README step 5)
for predictor in REVEL AM MP2; do
  poetry run python -m src.build_variant_reclassification_dataset --dedup --predictor "$predictor" \
    --output "data/output/reclassification/integrated_variant_effect_reclassification_${predictor}.tsv.gz"
done

# Official biobank analysis input file (one row per measurement; README step 6)
poetry run python -m src.build_variant_reclassification_dataset

# ExCALIBR-for-all-genes biobank analysis input file
poetry run python -m src.build_variant_reclassification_dataset \
  data/output/reclassification/integrated_variant_effect_dataset_analysis_excalibr_all_genes.csv.gz \
  --output data/output/reclassification/integrated_variant_effect_biobank_input_data_excalibr_all_genes.tsv.gz

# ExCALIBR-for-all-genes deduplicated reclassification dataset
poetry run python -m src.build_variant_reclassification_dataset \
  data/output/reclassification/integrated_variant_effect_dataset_analysis_excalibr_all_genes.csv.gz \
  --dedup \
  --output data/output/reclassification/integrated_variant_effect_reclassification_excalibr_all_genes.tsv.gz

# F9/TP53-assays-ExCALIBR/models-OddsPath biobank analysis input file
poetry run python -m src.build_variant_reclassification_dataset \
  data/output/reclassification/integrated_variant_effect_dataset_analysis_excalibr_assays_oddspath_models.csv.gz \
  --output data/output/reclassification/integrated_variant_effect_biobank_input_data_excalibr_assays_oddspath_models.tsv.gz

# ...and its deduplicated reclassification dataset
poetry run python -m src.build_variant_reclassification_dataset \
  data/output/reclassification/integrated_variant_effect_dataset_analysis_excalibr_assays_oddspath_models.csv.gz \
  --dedup \
  --output data/output/reclassification/integrated_variant_effect_reclassification_excalibr_assays_oddspath_models.tsv.gz
```

For the official checkpoint, `src/scripts/run_build_variant_reclassification_dataset.sh`
is the Docker-wrapped equivalent (same CLI, e.g. `... --dedup --output ...`)
-- see README.md's pipeline-stage-2 section. It only bind-mounts the
default checkpoint/CHEK2 paths, so pass an alternate checkpoint
positionally the same way as above if using it for one of the ExCALIBR
variants.

## `EXCALIBR_VINTAGE_OVERRIDE_GENES`

`BRCA1`/`PTEN`/`MSH2`/`TP53` prefer `ExC_points_2018` over `ExC_points_2025`
for the standalone `ExCALIBR_points` column, since their on-file ExCALIBR
calibrations are all 2018-vintage. This has no effect on the official
checkpoint -- its only TP53 dataset, `TP53_Fayer_2021_meta`, has no
ExCALIBR calibration in either vintage -- but matters for the two ExCALIBR
checkpoints above, several of whose restored TP53 datasets do have a
2018-vintage calibration.
