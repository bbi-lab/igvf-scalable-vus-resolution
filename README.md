# A scalable approach to resolving variants of uncertain significance

This repository contains the data-processing pipeline and figure-generation
code behind the paper *"A scalable approach to resolving variants of
uncertain significance."* It harmonizes variant-effect measurements from 90
multiplexed assays of variant effect (MAVE) datasets, maps them to genomic
coordinates, annotates them with ClinVar/gnomAD/predictor evidence, runs the
OddsPath/ACMG-AMP classification analysis described in the paper, and
generates the manuscript's main and extended data figures.

> **Preprint:** *A scalable approach to resolving variants of uncertain
> significance.* bioRxiv (2026).
> https://www.biorxiv.org/content/10.64898/2026.02.14.705848

Before any of that, see [Preparation](#preparation) for the tools and data
files every stage depends on.

The pipeline runs in three stages, each covered in its own section below --
with two external interim hand-offs in between, where data leaves this
repository for a process this codebase can't run itself, and the result
comes back before the next stage continues. A third, manual interim step
follows Stage 3, needed only when assembling a paper submission:

1. [Data preparation and variant annotation](#1-data-preparation-and-variant-annotation)
   - *Interim:* [ExCALIBR calibration](#interim-excalibr-calibration-external)
2. [Data analysis / variant classification and table preparation](#2-data-analysis--variant-classification-and-table-preparation)
   - *Interim:* [Biobank carrier odds-ratio analysis](#interim-biobank-carrier-odds-ratio-analysis-external)
3. [Figure preparation](#3-figure-preparation)
   - *Interim:* [Figure assembly](#interim-figure-assembly-manual-adobe-illustrator)

All three numbered stages run in Docker -- see
[Software Requirements](#software-requirements) below. The first two interim
steps are external and manual; see each section for why. The third (figure
assembly) is manual too, but requires a local Adobe Illustrator install
rather than Docker.

**Verifying your setup:** `scripts/smoke_test_variant_annotation.sh` and
`scripts/smoke_test_analysis.sh` run Stage 1 and Stage 2 end-to-end against
tiny real-data fixtures, so you can check your environment setup without a
multi-hour, multi-GB full run. See
[`docs/smoke_tests.md`](docs/smoke_tests.md).

**Unit tests:** `poetry run pytest` runs the unit test suite under `tests/`
(one `test_<module>.py` per `src/` script, exercising its core logic
directly rather than through the pipeline) -- faster and narrower in scope
than the smoke tests above, and intended for use during development rather
than as an environment check.

---

## Preparation

Do this once, before running any of the three numbered stages below.

### Environmental prerequisites

- **Docker and Docker Compose** — both `variant-annotation`'s steps and this
  project's own annotation steps run as Docker containers, not directly on
  the host.
- **`tabix` and `bgzip`** (both from htslib) — needed on the host itself,
  not in Docker, to prepare the AlphaMissense/REVEL predictor-score files
  below (`tabix` to index both; `bgzip` as part of REVEL's reformatting
  step). Install via your package manager (e.g. `apt install tabix`,
  `brew install htslib`) or alongside `samtools`.
- **A `variant-annotation` checkout** (Stage 1) — either the vendored
  submodule (`git submodule update --init vendor/variant-annotation`) or
  your own existing checkout that already has the large reference caches
  downloaded (see below), pointed to via
  `export VARIANT_ANNOTATION_DIR=/path/to/checkout`. See
  [`docs/variant_annotation_pipeline.md`](docs/variant_annotation_pipeline.md)
  for the full architecture and the submodule/override convention.
- **Poetry** (see [Software requirements](#software-requirements) below) —
  needed to build this project's own Docker image (`Dockerfile`/`compose.yaml`)
  for the steps it contributes to the pipeline.
- **A gnomAD Hail table cache** (Stage 1), built once per
  `variant-annotation` checkout via
  `scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache` (can
  take 6-7 hours) — not needed if `VARIANT_ANNOTATION_DIR` already points
  at a checkout where this has been done.

### Data-file prerequisites

Committed to this repo (no action needed):

- `data/input/maves/score_sets.tsv`, `Supplementary_Data_3.xlsx`,
  `Supplementary_Data_1_schema.xlsx` — this project's own dataset-name
  mapping, curation metadata, and Supplementary Data 1 column
  documentation (copy the latter to
  `data/output/supplementary_data/Supplementary_Data_1_schema.xlsx` before
  running `submit_supplementary_data_1_to_igvf.py`).
- `data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz` — RefSeq/Ensembl
  transcript mapping.
- `data/input/predictors/mp2_annotations.csv.gz`,
  `data/input/predictors/mp2_gene_symbol_map.tsv` — MutPred2 scores; our
  team ran MutPred2 to obtain annotations for this dataset's coding
  variants (see [`docs/data.md`](docs/data.md#where-mp2_annotationscsvgz-comes-from)).
- `data/input/predictors/supplementary_data_3_gene_entrez_ids.tsv` —
  upstream source `build_training_variant_files` consumes to build the
  MutPred2 side of the predictor training-overlap flags (see
  [`docs/build_training_variant_files.md`](docs/build_training_variant_files.md)).

(The exCALIBR calibration files and the biobank carrier odds-ratio table
are also already committed, but aren't needed until after Stage 1 and
Stage 2 respectively — see their own sections,
[Interim: ExCALIBR calibration](#interim-excalibr-calibration-external) and
[Interim: Biobank carrier odds-ratio analysis](#interim-biobank-carrier-odds-ratio-analysis-external),
for what each file is.)

Not committed — you must populate, download, or build these yourself, all
for Stage 1:

1. **`data/input/maves/cvfg_variants.0.tsv`** — this project's own variant
   list, Step 1's input. Without access to a local MaveDB Postgres mirror,
   use the public-API script, which writes directly to
   `data/input/maves/cvfg_variants.0.tsv` by default:

   ```bash
   poetry run python -m src.fetch_mavedb_scores
   ```

   With direct access to the MaveDB database (internal to our team — most
   end users won't have this, but it's documented here for completeness),
   `src/mavedb_scores.sql` is the source of truth instead (run with
   `psql -f src/mavedb_scores.sql your_mavedb_db`, then export the final
   `SELECT` to `data/input/maves/cvfg_variants.0.tsv`).

   Either way, once the file is at `data/input/maves/cvfg_variants.0.tsv`,
   `scripts/run_variant_annotation_pipeline.sh` stages it into
   `data/intermediate/variant_annotation/data/` automatically — no manual
   copy step needed. See
   [Where `cvfg_variants.0.tsv` comes from](docs/data.md#where-cvfg_variants0tsv-comes-from)
   in `docs/data.md` for full details, known differences between the two
   routes, and why the API version isn't yet validated against the SQL
   version's output. (A Zenodo-hosted copy may be made available later for
   convenience; none exists yet.)

2. **AlphaMissense predictor scores** (`data/input/predictors/AlphaMissense_hg38.tsv.gz`
   + `.tbi`) — download and tabix-index as-is, no reformatting needed:

   ```bash
   curl -o data/input/predictors/AlphaMissense_hg38.tsv.gz \
     https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
   tabix -s 1 -b 2 -e 2 data/input/predictors/AlphaMissense_hg38.tsv.gz
   ```

3. **REVEL predictor scores** (`data/input/predictors/revel_hg38.tsv.gz` +
   `.tbi`) — unlike AlphaMissense, REVEL's own download needs reformatting
   first: this pipeline's `--revel-mode aa` needs `aaref`/`aaalt`/
   `ensembl_transcriptid` columns laid out in a specific order, which the
   raw REVEL distribution doesn't do directly:

   ```bash
   # Download revel_with_transcript_ids.csv.zip from
   # https://sites.google.com/site/revelgenomics/downloads, then:
   unzip revel_with_transcript_ids.csv.zip
   tail -n +2 revel_with_transcript_ids.csv \
     | awk -F',' 'NF>=9 && $3!="" && $3!="." \
                  {print $1"\t"$3"\t"$4"\t"$5"\t"$8"\t"$6"\t"$7"\t"$9}' \
     | (printf '#chr\tpos\tref\talt\trevel_score\taaref\taaalt\tensembl_transcriptid\n'; sort -k1,1V -k2,2n) \
     | bgzip > data/input/predictors/revel_hg38.tsv.gz
   tabix -s 1 -b 2 -e 2 -S 1 data/input/predictors/revel_hg38.tsv.gz
   ```

   See `vendor/variant-annotation/docs/annotate_predictors.md`'s "Data file
   preparation" and "REVEL disambiguation modes" sections for why, and for
   the other `--revel-mode` options.

4. **MutPred2 training-variant source**
   (`data/input/predictors/mp2_actual_training_data.txt.gz`) — feeds
   `build_training_variant_files`, which Step 16 runs automatically with no
   path override, so the file must land at exactly this path. The source
   itself isn't gzip-compressed, so pipe it through `gzip` on the way in:

   ```bash
   curl -s https://mutpred.mutdb.org/wo_exclusive_hgmd_mp2_training_data.txt \
     | gzip > data/input/predictors/mp2_actual_training_data.txt.gz
   ```

   This is a publicly-shareable subset of MutPred2's own training set (it
   omits some variants MutPred2 can't redistribute due to a third-party
   license). If you separately have access to the full training dataset,
   use that instead, at the same destination path — see
   [`docs/build_training_variant_files.md`](docs/build_training_variant_files.md#getting-the-mutpred2-training-variant-source)
   for details and for what the substitute above is missing.

5. **A REVEL training-overlap source** (`build_training_variant_files.py
   --revel-source`) — *optional.* This project's own copy has unconfirmed
   provenance (see [`docs/data.md`](docs/data.md#where-the-predictor-training-variant-sources-come-from)),
   so it isn't distributed. If you don't provide one, the pipeline skips
   this one circularity check (zero variants get flagged as REVEL-training
   overlaps) rather than failing — see
   [`docs/build_training_variant_files.md`](docs/build_training_variant_files.md#running-without-a-revel-training-overlap-source)
   for what to provide if you have a suitable source and want the check
   re-enabled.

6. **A `variant-annotation` checkout with the large, generic reference
   caches already downloaded**: SpliceAI VCFs, dbNSFP, `clinvar_cache/`,
   and the gnomAD Hail table cache above. These are tens of GB and are
   deliberately not duplicated into this repo — follow `variant-annotation`'s
   own setup to populate them in whichever checkout `VARIANT_ANNOTATION_DIR`
   points at.

---

## 1. Data preparation and variant annotation

The first stage turns curated MAVE inputs into the integrated variant-effect
dataset (Supplementary Data 1) that everything downstream reads from. It used
to be a single notebook (`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`,
now archived); that notebook has been replaced by a Dockerized pipeline built
on the sibling **[variant-annotation](https://github.com/bbi-lab/variant-annotation)**
project.

### Relation to the variant-annotation pipeline

`variant-annotation` is a general-purpose, Dockerized variant-mapping and
annotation pipeline (reverse translation, ClinVar/gnomAD/SpliceAI/ClinGen
evidence repository/VEP/MaveDB/predictor annotation) that is developed and
versioned independently of this project. This repo depends on it rather than
reimplementing it:

- `vendor/variant-annotation` is a git submodule pinned to a known-good
  commit of `variant-annotation`.
- `scripts/variant_annotation_pipeline.sh` drives that pipeline's steps
  against this project's MAVE data, interleaving a handful of this project's
  own steps (`postprocess_mavedb_functional_classifications`,
  `add_mavedb_active_calibration_columns`, `derive_score_set_urn`,
  `annotate_simplified_consequence`, `recalculate_clingen_classification`,
  `flag_variants`) between `variant-annotation`'s numbered steps.
- `scripts/run_variant_annotation_pipeline.sh` is the orchestrator you
  actually run: it stages this project's inputs, points the
  `variant-annotation` checkout at them via `VARIANT_DATA_DIR`, runs the
  pipeline, and collects the final gzipped output.

Full architecture, the submodule/override convention, and a non-obvious
`VARIANT_DATA_DIR` path-mapping subtlety are documented in
[`docs/variant_annotation_pipeline.md`](docs/variant_annotation_pipeline.md)
— read that before touching `scripts/variant_annotation_pipeline.sh` or its
wrapper scripts.

Environmental and data-file prerequisites for this stage are covered in
[Preparation](#preparation) above.

### Running it

```bash
git submodule update --init vendor/variant-annotation   # if using the vendored submodule
poetry install --all-extras

# One-time per variant-annotation checkout (skip if already done):
scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache

# Full run:
scripts/run_variant_annotation_pipeline.sh
```

This produces `data/output/maves/integrated_variant_effect_dataset.tsv.gz`
(expanded, one row per DNA-level variant) and its `.condensed.tsv.gz`
sibling (one row per measurement). Pass `--step N` to (re)run a single
numbered step for debugging — see
[`docs/variant_annotation_pipeline.md`](docs/variant_annotation_pipeline.md)
for the full step list and every per-step doc it links to.

---

## Interim: ExCALIBR calibration (external)

Between Stage 1 and Stage 2, `data/output/maves/integrated_variant_effect_dataset.tsv.gz`
is shared outside this repository with the ExCALIBR calibration process --
a separate analysis pipeline, maintained elsewhere, that fits per-dataset
functional-evidence calibrations. Nothing in this repo can run that step
itself.

ExCALIBR calibration returns four artifacts, all checked into
`data/input/mave_calibration/`:

- Per-dataset calibration summaries (`excalibr/json/`), loaded into
  Supplementary Data 4's `ExCALIBR_calibrations` sheet by Stage 2's
  `src/load_excalibr_calibrations.py` below.
- Its own copy of the integrated dataset, with `excalibr_*` evidence
  columns added (`integrated_variant_effect_dataset_with_excalibr.csv.gz`)
  -- used just below to finish Supplementary Data 1, and by
  `src/filter_pillar_data.py` (Extended Data Figure 10 prep; see
  [`docs/gnomad_ancestry_resolution.md`](docs/gnomad_ancestry_resolution.md)).
- Gene-specific predictor calibration thresholds, pre-populated into
  `Supplementary_Data_4.seed.xlsx`'s three externally-sourced sheets --
  copy this to `data/output/supplementary_data/Supplementary_Data_4.xlsx`
  as your starting point in a fresh checkout.
- Figure 4's out-of-bag export and seed cache (`excalibr/fig4_data_oob.tar.gz`,
  `excalibr/figure4_data.seed.json.gz`) -- see
  [`docs/build_figure4_data.md`](docs/build_figure4_data.md).

See [`docs/load_excalibr_calibrations.md`](docs/load_excalibr_calibrations.md)
and [`docs/data.md`](docs/data.md) for exactly where each file lives and how
it's produced.

### Adding ExCALIBR evidence columns (Supplementary Data 1)

A separate, manual step folds in per-variant ExCALIBR likelihood-ratio
evidence — `excalibr_prior`, `excalibr_posterior`, `excalibr_lr_plus`,
`excalibr_points`, `excalibr_acmg_evidence_code`, `excalibr_filter_reason`,
`excalibr_clinvar_release` — onto the integrated dataset above. Publishing
the raw prior/posterior/likelihood-ratio values (not just the points ExCALIBR
resolved them to) lets anyone recompute ACMG evidence strength under a prior
probability of their own choosing, rather than being limited to the one this
pipeline's own classification used.

This step isn't part of `scripts/run_variant_annotation_pipeline.sh` because
its input — the ExCALIBR calibration run's own copy of the integrated
dataset,
`data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz`
above (checked in; see [`docs/add_excalibr_columns.md`](docs/add_excalibr_columns.md)
for where it comes from) — only exists once the interim hand-off above has
returned. Once both that file and
`data/output/maves/integrated_variant_effect_dataset.tsv.gz` exist, run:

```bash
poetry run python -m src.add_excalibr_columns
```

This writes
`data/output/maves/integrated_variant_effect_dataset_with_excalibr.tsv.gz`
and copies it to `data/output/supplementary_data/Supplementary_Data_1.tsv.gz`
— the version of Supplementary Data 1 actually distributed with the paper.
See [`docs/add_excalibr_columns.md`](docs/add_excalibr_columns.md).

---

## 2. Data analysis / variant classification and table preparation

The second stage takes the integrated variant-effect dataset from Stage 1
and produces the manuscript's classification tables: OddsPath calibrations,
ACMG/AMP evidence-point assignment, and the final per-gene, per-predictor
classification files (Supplementary Data 5).

This stage runs as Jupyter notebooks under `notebooks/analysis/`, executed
inside Docker via the `analysis-notebooks` Compose service (same approach as
Stage 1 and the R figures) — no local Poetry/Jupyter setup needed to
reproduce results; see [Software Requirements](#software-requirements)
below only if you want to develop/debug a notebook interactively. Each notebook has a
companion `README_*.md` in the same directory documenting its inputs,
methods, and outputs in detail — read those before running or modifying
one.

Run in this order:

1. **`OddsPath_calculations.ipynb`** — OddsPath likelihood-ratio
   calculations per dataset, following Brnich et al. (2019).
2. **`src/load_oddspath_calibrations.py`** — loads that notebook's output
   (`data/output/mave_calibration/OddsPath_calibrations.csv.gz`) into the
   `OddsPath_calibrations` sheet of Supplementary Data 4, analogous to
   `src/load_excalibr_calibrations.py` below. See
   [`docs/load_oddspath_calibrations.md`](docs/load_oddspath_calibrations.md).
3. **`OddsPath_classifications.ipynb`** — turns those OddsPath values into
   ACMG/AMP evidence points.
4. **`Variant_Classification_analysis.ipynb`** — combines functional
   evidence (via exCALIBR calibrations or OddsPath) with computational
   predictor evidence (REVEL, AlphaMissense, MutPred2; genome-wide and
   gene-specific), assigns final classifications, deduplicates variants
   seen in multiple assays, and exports the per-category classification
   files that become Supplementary Data 5. See
   [`notebooks/analysis/README_Variant_Classification_analysis.md`](notebooks/analysis/README_Variant_Classification_analysis.md).
5. **`src/build_variant_reclassification_dataset.py --dedup --predictor
   <REVEL|AM|MP2>`** — reads that notebook's own cell-69 checkpoint,
   re-applies its downstream exclusions, adds six ACMG/AMP points columns,
   and collapses to one row per DNA variant (greatest `abs(Combined_points)`
   wins, ties broken by `Dataset` name). Run once per computational
   predictor to produce a deduplicated reclassification dataset scored by
   each in turn (`Combined_points` is predictor-dependent, so the
   surviving measurement for a given variant can differ by predictor —
   see [`docs/build_variant_reclassification_dataset.md`](docs/build_variant_reclassification_dataset.md)).
   Produces `data/output/reclassification/integrated_variant_effect_reclassification_<REVEL|AM|MP2>.tsv.gz`.

(The same script's default invocation -- no `--dedup` -- builds the biobank
analysis input file; see
[Interim: Biobank carrier odds-ratio analysis](#interim-biobank-carrier-odds-ratio-analysis-external)
below, after this stage.)

Two more `src/` scripts support this stage (beyond `load_oddspath_calibrations.py`
above) and are already converted:

- **`src/load_excalibr_calibrations.py`** — loads exCALIBR JSON calibrations
  (`data/input/mave_calibration/excalibr/json/`) into Supplementary Data 4
  so the workbook stays in sync with the latest calibration run. See
  [`docs/load_excalibr_calibrations.md`](docs/load_excalibr_calibrations.md).
- **`src/mave_dataset_stats.py`** — reports summary statistics (dataset,
  measurement, and variant counts; predictor score coverage; clinical
  attribute breakdowns; ExCALIBR calibration coverage; reclassification
  agreement) over the integrated dataset. See
  [`docs/mave_dataset_stats.md`](docs/mave_dataset_stats.md).

### Running it

```bash
# One-time: build the tools image (adds the "notebooks" Poetry extra --
# ipykernel/jupyterlab -- and Arial for matplotlib; see Dockerfile)
docker compose build mave-dataset-stats

# Refresh Supplementary Data 4 from the latest exCALIBR calibrations
src/scripts/run_load_excalibr_calibrations.sh

# 1. OddsPath likelihood-ratio calculations
src/scripts/run_notebook.sh --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=600 \
  --output executed_OddsPath_calculations.ipynb \
  notebooks/analysis/OddsPath_calculations.ipynb

# 2. Refresh Supplementary Data 4's OddsPath_calibrations sheet from that run
src/scripts/run_load_oddspath_calibrations.sh

# 3-4. Remaining notebooks above, in order
for nb in Variant_Classification_analysis OddsPath_classifications; do
  src/scripts/run_notebook.sh --to notebook --execute \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    notebooks/analysis/${nb}.ipynb
done

# 5. Build a deduplicated reclassification dataset per computational predictor
for predictor in REVEL AM MP2; do
  src/scripts/run_build_variant_reclassification_dataset.sh --dedup --predictor "$predictor" \
    --output "data/output/reclassification/integrated_variant_effect_reclassification_${predictor}.tsv.gz"
done
```

Each `nbconvert --execute` leaves a side-effect `executed_<name>.ipynb` next
to the original (nbconvert's copy with outputs attached) — gitignored, so no
need to delete them. `kernel_name=python3` is `ipykernel`'s default kernel
inside the container, not a registered project-specific one — nothing in
the repo depends on a particular kernel name.

Outputs land under `data/output/supplementary_data/` (`Supplementary_Data_4.xlsx`,
`Supplementary_Data_5.xlsx`), `data/output/predictor_calibration/` (the
per-gene control files also used by Extended Data Figure 4), and
`data/output/reclassification/` (the per-predictor deduplicated
reclassification datasets from step 5).

---

## Interim: Biobank carrier odds-ratio analysis (external)

Between Stage 2 and Stage 3, a subset of figures (Figure 2i, Figure 6b,
Extended Data Figure 2) need carrier odds ratios computed against a
biobank cohort -- something this repo can't do itself, since it requires
individual-level biobank genotype/phenotype data under access
restrictions this codebase doesn't (and shouldn't) have a path to.

1. Build the biobank analysis input file from Stage 2's
   `Variant_Classification_analysis.ipynb` checkpoint:

   ```bash
   src/scripts/run_build_variant_reclassification_dataset.sh
   ```

   Produces `data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz`
   (one row per surviving variant-effect measurement, not deduplicated to
   one row per DNA variant -- see
   [`docs/build_variant_reclassification_dataset.md`](docs/build_variant_reclassification_dataset.md)).
2. That file is handed off manually to the biobank-side analysis --
   currently the All of Us Research Program's Researcher Workbench, an
   access-restricted environment -- which computes carrier log-odds-ratios
   per gene/dataset against its cohort. Nothing in this repo can run or
   automate this step.
3. The result is checked back into `data/input/biobank/` (currently
   `AoU-OR-estimates_2026-10-07.tsv.gz`, a stand-in pending a long-term,
   citable re-export of the same analysis) and consumed by the three
   figures above via `src/lib/figures/load_biobank_or_table.R`.

See [`docs/biobank_analysis.md`](docs/biobank_analysis.md) for the full
workflow.

---

## 3. Figure preparation

The third stage generates the manuscript's main and extended data figures
from Stage 2's outputs. Each figure directory under `notebooks/figures/`
contains its own plotting notebook(s)/script(s) and any small supporting
inputs.

Two toolchains are in play, depending on the figure, both Dockerized:

- **Python notebooks/scripts** (Altair/matplotlib), run via the
  `analysis-notebooks` Compose service (`jupyter nbconvert --execute`) or a
  dedicated per-script service (e.g. `build-figure3-data`) -- same tools
  image as Stage 2, no local Poetry/Jupyter setup needed.
- **R scripts/`.Rmd` files** (tidyverse, ggplot2 extensions like `ggsankey`,
  `patchwork`, `ggh4x`), run via the `r-figures` Docker Compose service
  (`Dockerfile.r`) rather than a local R install, since they depend on
  system fonts (Arial) and `cairo_pdf` output that are easiest to reproduce
  in a container.

A few figures (Figure 2i, Figure 6b, Extended Data Figure 2) additionally
need a biobank carrier odds-ratio table -- see
[Interim: Biobank carrier odds-ratio analysis](#interim-biobank-carrier-odds-ratio-analysis-external)
above.

### Running it

One-time setup:

```bash
docker compose build mave-dataset-stats r-figures
```

#### Figure 2

```bash
# 1. Prep notebook (run first; everything else in this figure depends on it)
src/scripts/run_notebook.sh --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=600 \
  --output executed_PP_ProcessBigDataFrame.ipynb \
  notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb

# 2. Panel notebooks (independent of each other; all read step 1's output)
for nb in PP_ClinVarPrecisionRecall PP_Fig2_Heatmaps PP_ResolutionOverview PP_StackedHistograms; do
  src/scripts/run_notebook.sh --to notebook --execute \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    notebooks/figures/figure_2/${nb}.ipynb
done

# 3. Figure_2i.R -- reads the biobank carrier odds-ratio table (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_2 \
  r-figures Figure_2i.R
```

#### Figure 3

```bash
# Rebuild Figure3a/c/d.csv.gz (gitignored intermediates, not committed) --
# see docs/build_figure3_data.md.
src/scripts/run_build_figure3_data.sh

# Writes executed_curation_summary_figure3.html next to the .Rmd (gitignored)
# and PNGs/an SVG to data/output/figures/assets/figure_3/.
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_3 \
  r-figures -e 'rmarkdown::render("curation_summary_figure3.Rmd", output_file = "executed_curation_summary_figure3.html")'
```

See [`docs/figures.md`](docs/figures.md#figure-3) for what each panel shows
and [`docs/build_figure3_data.md`](docs/build_figure3_data.md) for the data
pipeline behind it.

#### Figure 4

```bash
# 1. Rebuild figure4_data.json.gz, carrying forward fields that can't be
#    regenerated (see docs/build_figure4_data.md)
src/scripts/run_build_figure4_data.sh \
  --cached-json data/input/mave_calibration/excalibr/figure4_data.seed.json.gz

# 2. Execute the notebook to produce data/output/figures/assets/figure_4/
#    (figure_4.png, figure_4.pdf, and the 6.5in-wide figure_4_calibrated.pdf
#    placed into figures_2-6.ai)
src/scripts/run_notebook.sh --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=600 \
  --output executed_figure4.ipynb \
  notebooks/figures/figure_4/figure4.ipynb
```

#### Figure 5/6

```bash
# Figure_6b.R -- reads the biobank carrier odds-ratio table (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures Figure_6b.R

# Figure5_6.Rmd: writes executed_Figure5_6.html next to the .Rmd (gitignored)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures -e 'rmarkdown::render("Figure5_6.Rmd", output_file = "executed_Figure5_6.html")'

# Optional: same Figure5_6.Rmd, restricted to missense_variant rows only.
# Writes to figure_5_missense/ and figure_6_missense/ instead of
# figure_5/figure_6, so it doesn't overwrite the default run above.
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures -e 'rmarkdown::render("Figure5_6.Rmd", output_file = "executed_Figure5_6_missense.html", params = list(consequence_filter = "missense"))'
```

#### Extended Data Figure 2

```bash
# Extended_Data_Figure_2.R -- reads the biobank carrier odds-ratio table (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_2 \
  r-figures Extended_Data_Figure_2.R
```

#### Extended Data Figure 4

```bash
src/scripts/run_notebook.sh --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=600 \
  --output executed_extended_data_figure_4.ipynb \
  notebooks/figures/extended_data_figure_4/Extended_Data_Figure_4.ipynb
```

#### Extended Data Figures 3, 3alt, 3alt missense, 5, 5alt, 6alt, 8, 9 (`Extended_data_figures.Rmd`)

```bash
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_3_5_6_8_9 \
  r-figures -e 'rmarkdown::render("Extended_data_figures.Rmd")'
```

#### Extended Data Figure 6

```bash
src/scripts/run_make_extended_data_figure_6.sh

# Optional: same chart, restricted to missense_variant rows only. Writes
# new_classification_heatmap_missense.pdf into the same
# extended_data_figure_6/ dir, so it doesn't overwrite the all-consequences
# run above.
src/scripts/run_make_extended_data_figure_6.sh --consequence-filter missense
```

OddsPath-calibrated variant, visualizing `Supplementary_Data_6.xlsx`'s
VUS/gnomAD/Unobserved sheets instead (see
[`docs/figures.md`](docs/figures.md#extended-data-figure-6-oddspath-variant-srcmake_extended_data_figure_6_oppy)):

```bash
src/scripts/run_make_extended_data_figure_6_op.sh
```

#### Extended Data Figure 7 (functional-vs-predictor evidence ablation)

```bash
# ablation_missense.pdf is the version actually used in the manuscript figure.
src/scripts/run_ablation_variant_reclassification.sh \
  --consequence missense_only \
  --calibrated-figure data/output/figures/assets/extended_data_figure_7/ablation_missense.pdf

# All-variants version -- generated but not used in the manuscript. Writes to
# ablation.pdf instead, so it doesn't overwrite the missense-only run above.
src/scripts/run_ablation_variant_reclassification.sh \
  --calibrated-figure data/output/figures/assets/extended_data_figure_7/ablation.pdf
```

#### Extended Data Figure 10 (classification resolution rate by gnomAD ancestry group)

Reads `data/input/gnomAD/vus_with_ancestry_af.csv` and
`gnomad_with_ancestry_af.csv` directly -- if these aren't already present (a
fresh checkout won't have them; they aren't produced by any other step
above), prepare them first per
[`docs/gnomad_ancestry_resolution.md`](docs/gnomad_ancestry_resolution.md).
That requires a GRCh38 reference FASTA and gnomAD v4.1 sites VCFs, which
aren't included in this repo and must be downloaded separately (see that doc
for sources).

```bash
src/scripts/run_make_extended_data_figure_10.sh
```

Every `nbconvert --execute` call above leaves a side-effect `executed_<name>.ipynb`
copy next to the source, and `Figure5_6.Rmd`'s `rmarkdown::render()` call
does the same with `executed_Figure5_6.html` — both gitignored, so no need
to delete them. `Extended_data_figures.Rmd` above is the one exception: it
still renders the plain, non-gitignored `Extended_data_figures.html` next to
itself — delete that manually if you don't want it in the working tree. See
[`docs/figures.md`](docs/figures.md) for output paths for every
individual panel and several figure-specific gotchas (date-stamped
filenames, figures that only reconstruct from a cached intermediate, output
directories that aren't created automatically).

---

## Interim: Figure assembly (manual, Adobe Illustrator)

Only needed by the project team when preparing a paper submission -- not
part of a normal pipeline run, and the only step in this repo that isn't
Dockerized or Poetry-managed. Two checked-in Illustrator documents,
`data/output/figures/illustrator/figures_2-6.ai` and
`extended_data_figures.ai`, assemble every figure Stage 3 generates assets
for, by placing those assets onto per-figure artboards;
`export_artboards.jsx` in the same directory then exports each artboard to
its own PNG or PDF.

Figures 1 and 7, and Extended Data Figures 1, 11, and 12, aren't generated
by this repo at all -- they're checked into
`data/input/external_figure_assets/` as provided externally. Figures 1 and 7
are each their own standalone Illustrator document, exported the same way as
the two above; Extended Data Figures 1, 11, and 12 have no `.ai` of their
own here -- their PDFs are placed directly into `extended_data_figures.ai`
and come out as part of its own export.

See [`docs/figure_assembly.md`](docs/figure_assembly.md) for the full
process: resolving the `.ai` files' absolute-path links (Illustrator has no
relative-link mechanism), the PNG/PDF export options and which PDF standard
to use, and the two artboards (with and without a caption) each Extended
Data Figure exports.

---

## Data Availability

The integrated variant effect dataset (Supplementary Data 1) is too large
for GitHub and is hosted externally on Zenodo. Every other Supplementary
Data file is generated by running this pipeline (see
[Section 2](#2-data-analysis--variant-classification-and-table-preparation))
rather than distributed pre-built.

**Zenodo record:** https://doi.org/10.5281/zenodo.23308563

Files included in this GitHub repository are:

- Analysis scripts and notebooks
- Figure-generation scripts and notebooks
- Small supporting inputs

## License

This project's own code is released under the [BSD 3-Clause
License](LICENSE). See below for the licenses of third-party data files
bundled in the repository.

## Third-Party Data & Licenses

`data/input/genes/` bundles three external gene-level reference files used to
build Figure 3's bubble plot (see [`docs/data.md`](docs/data.md) for exact
download dates/URLs and [`docs/build_figure3_data.md`](docs/build_figure3_data.md)
for how they're used and how to refresh them):

- **GenCC** gene-disease validity submissions
  (`data/input/genes/gencc-submissions.csv.gz`) — [CC0 1.0 Universal Public
  Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/);
  GenCC requests attribution to GenCC and its contributing sources.
- **UniProtKB/Swiss-Prot** reviewed human proteome
  (`data/input/genes/uniprotkb_9606_reviewed.tsv.gz`) — [CC BY
  4.0](https://creativecommons.org/licenses/by/4.0/). Attribution: this work
  uses data from UniProt, which is distributed under a Creative Commons
  Attribution 4.0 International License.
- **NCBI Genetic Testing Registry (GTR)** test/condition/gene export
  (`data/input/genes/test_condition_gene.txt.gz`) — public domain (US
  government work); NLM requests attribution to GTR as a data source.

This pipeline also reads predictor scores and reference annotations from
several other third-party sources, each under its own license:

- **AlphaMissense** (`data/input/predictors/AlphaMissense_hg38.tsv.gz`) —
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) as of the
  database's March 2024 relicensing (an earlier, non-commercial-only release
  has since been superseded). Attribution: Cheng et al., "Accurate
  proteome-wide missense variant effect prediction with AlphaMissense,"
  *Science* (2023). ⚠️ The file's own in-content header comment still says
  `Licensed under CC BY-NC-SA 4.0 license` — stale boilerplate DeepMind
  never updated in the file itself when they relicensed it; confirmed
  byte-identical to the file currently served from their CC-BY-4.0-licensed
  [Zenodo record](https://zenodo.org/records/10813168) (same size, matching
  MD5 on a byte-range check), so the file's own comment is wrong, not the
  license. Trust the Zenodo record's license, not the file's text.
- **REVEL** (`data/input/predictors/revel_hg38.tsv.gz`,
  `clustering_variants_revel_training_overlap.csv.gz`) — free for
  non-commercial use only; commercial use requires contacting the developer
  (Weiva Sieh) directly. Attribution: Ioannidis et al., "REVEL: An Ensemble
  Method for Predicting the Pathogenicity of Rare Missense Variants," *AJHG*
  (2016).
- **MutPred2** (`data/input/predictors/mp2_annotations.csv.gz`) — no
  separate data license is published by the MutPred2 team beyond the MIT
  license on their source code; cite the source publication when using
  these scores.
  Attribution: Pejaver et al., "Inferring the molecular and phenotypic impact
  of amino acid variants with MutPred2," *Nature Communications* (2020).
- **MANE** (`data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz`) — a
  joint NCBI/EMBL-EBI product with no separately published license; both
  contributing agencies make their own content freely available (NCBI
  content is public domain as a US government work, EMBL-EBI/Ensembl content
  is CC BY 4.0). Attribution: Morales et al., "A joint NCBI and EMBL-EBI
  transcript set for clinical genomics and research," *Nature* (2022).
- **Ensembl** (`data/input/reference/extended_ensembl_consequence.csv.gz`,
  derived from Ensembl VEP consequence terms) — [CC BY
  4.0](https://creativecommons.org/licenses/by/4.0/).
- **gnomAD** (ancestry-group allele frequencies used for Extended Data
  Figure 10 and as a control population throughout) — [CC0 1.0 Universal
  Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/);
  gnomAD requests attribution and a link back to
  [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/) wherever
  possible.
- **ClinVar** (used throughout for variant pathogenicity annotation and as a
  classification truth set) — public domain (NCBI/US government work); no
  restriction on use or redistribution.

This section covers only the sources listed above. Other third-party data
this pipeline reads (MaveDB score sets, the IGVF data portal) should be used
and cited per each source's own terms; consult the upstream provider for
current license terms before redistributing this repository's derived
outputs.

## Software Requirements

- **Docker and Docker Compose** — the only requirement to reproduce
  results. All three stages run in Docker: the variant-annotation pipeline
  and this project's own steps (Stage 1), the analysis notebooks and
  per-script figure builders (Stages 2-3, `python:3.12-slim`-based tools
  image), and the R figure scripts (Stage 3, `r-figures` service, R 4.4.2
  via `rocker/r-ver`). Neither a local Python nor R install is required.
- **Python 3.12** (optional) — only needed to develop/debug a notebook or
  script interactively outside Docker. Dependencies are Poetry-managed
  (`pyproject.toml`):
  ```bash
  poetry env use python3.12   # one-time, only if `poetry env use` picks the wrong interpreter
  poetry install --all-extras
  ```
  This creates an in-project `.venv/` with the pipeline dependencies plus
  the `dev` (Ruff, nbstripout, pre-commit), `tests` (pytest), and
  `notebooks` (JupyterLab/ipykernel) extras.
- Additional per-script package dependencies are declared in
  `pyproject.toml` (Python) or the top of each `.R`/`.Rmd` file (R).
