# gnomAD Ancestry Resolution Workflow

The full, reproducible path from raw gnomAD v4.1 data to Extended Data
Figure 10 (classification resolution rate by gnomAD ancestry group). Four
steps:

1. **Obtain the gnomAD reference data** -- a GRCh38 FASTA and gnomAD v4.1
   sites VCFs, placed under `data/input/gnomAD/`.
2. **`src/filter_pillar_data.py`** ([docs](filter_pillar_data.md)) -- selects
   the VUS and gnomAD category datasets from the ExCALIBR calibration
   output and classifies each variant (`Class_REVEL`).
3. **`src/map_gnomad_ancestry.py`** ([docs](map_gnomad_ancestry.md)) --
   annotates each category's variants with gnomAD v4.1 exome
   ancestry-group allele frequencies.
4. **`src/make_extended_data_figure_10.py`** ([docs](figures.md#extended-data-figure-10-srcmake_extended_data_figure_10py))
   -- computes the resolution-rate statistics and draws the figure.

Steps 2 and 3 are manual (not Dockerized, like `src/add_excalibr_columns.py`
-- see `docs/add_excalibr_columns.md`): step 3 needs the large external
gnomAD data from step 1, which doesn't fit this project's `docker compose`
bind-mount convention (the committed repo tree plus an optional
`VARIANT_DATA_DIR` staging mount), and both are infrequent, one-off
data-refresh steps rather than something run on every pipeline invocation.

## Where files live

- **`data/input/gnomAD/`** -- gnomAD-related data that's either
  user-downloaded or a checked-in snapshot, never a script's working
  output:
  - `reference/`, `vcf/` -- the GRCh38 FASTA and gnomAD sites VCFs you
    download in step 1 (gitignored -- hundreds of GB, never committed; see
    `.gitignore`'s gnomAD section).
  - `vus_with_ancestry_af.csv`, `gnomad_with_ancestry_af.csv` -- the
    checked-in snapshot `src/make_extended_data_figure_10.py` reads by
    default. These are reproducible via steps 2-3 below, but deliberately
    stay a committed snapshot rather than a rebuilt-every-run pipeline
    output (same reasoning as `data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz`
    -- see `docs/add_excalibr_columns.md`); refresh them by copying over the
    result of step 3 when you want to update the figure's input data.
- **`data/intermediate/figures/extended_data_figure_10/`** -- gitignored
  working directory for everything steps 2 and 3 themselves produce
  (`src/filter_pillar_data.py`'s category CSVs, `src/map_gnomad_ancestry.py`'s
  ancestry-annotated CSVs). Nothing here is committed or read directly by
  `src/make_extended_data_figure_10.py`.

## 1. Obtain the gnomAD reference data

`src/map_gnomad_ancestry.py` needs two things not checked into this repo,
both conventionally placed under `data/input/gnomAD/` (its CLI defaults
point there, so once these are in place you don't need to pass
`--ref-fasta`/`--exomes-dir` explicitly):

- **An indexed GRCh38 reference FASTA**, at
  `data/input/gnomAD/reference/GRCh38.fa` (`.fai` alongside it). Any
  standard GRCh38 primary-assembly FASTA works, as long as its contig names
  match gnomAD's (`chr1`, `chr2`, ... `chrX`, `chrY`, `chrM`) -- e.g. the
  Broad's `Homo_sapiens_assembly38.fasta` reference bundle, or UCSC's
  `hg38.fa`. Index it once with:

  ```bash
  samtools faidx data/input/gnomAD/reference/GRCh38.fa
  ```

- **gnomAD v4.1 sites VCFs**, bgzipped and tabix-indexed, one per
  chromosome, at `data/input/gnomAD/vcf/exomes/` (required) and
  `data/input/gnomAD/vcf/genomes/` (optional -- recorded for reference, not
  used by the ancestry-group analysis), from gnomAD's public downloads
  (<https://gnomad.broadinstitute.org/downloads> -> v4.1 -> "Sites" ->
  GRCh38/hg38). Either download through the browser/`gsutil`, or pull
  directly from the public bucket, e.g.:

  ```bash
  gsutil -m cp \
    'gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr*.vcf.bgz*' \
    data/input/gnomAD/vcf/exomes/
  ```

  (repeat for `vcf/genomes/` if genome frequencies are wanted). Each
  `.vcf.bgz` needs its `.tbi` alongside it -- the public bucket already
  provides these; if building your own subset, index with
  `tabix -p vcf <file>.vcf.bgz`.

These are large (the full exomes release is hundreds of GB); a partial
download covering only the chromosomes your variant set touches is fine --
`src/map_gnomad_ancestry.py` looks up one chromosome's VCF at a time and
logs a warning (not an error) for any chromosome file it can't open.

## 2. Select and classify the category datasets

Run once per category needed (`VUS` and `gnomAD` for Extended Data Figure
10):

```bash
poetry run python -m src.filter_pillar_data \
    --points-input data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz \
    --supp4-input data/output/supplementary_data/Supplementary_Data_4.xlsx \
    --supp3-input data/input/maves/Supplementary_Data_3.xlsx \
    --category VUS

poetry run python -m src.filter_pillar_data \
    --points-input data/input/mave_calibration/integrated_variant_effect_dataset_with_excalibr.csv.gz \
    --supp4-input data/output/supplementary_data/Supplementary_Data_4.xlsx \
    --supp3-input data/input/maves/Supplementary_Data_3.xlsx \
    --category gnomAD
```

Writes `data/intermediate/figures/extended_data_figure_10/replicated_VUS_category_REVEL.csv`
and `replicated_gnomAD_category_REVEL.csv` (see
[`docs/filter_pillar_data.md`](filter_pillar_data.md) for the full column/flag
reference).

## 3. Annotate with gnomAD ancestry allele frequencies

```bash
poetry run python -m src.map_gnomad_ancestry \
    --input data/intermediate/figures/extended_data_figure_10/replicated_VUS_category_REVEL.csv \
    --output data/intermediate/figures/extended_data_figure_10/vus_with_ancestry_af.csv

poetry run python -m src.map_gnomad_ancestry \
    --input data/intermediate/figures/extended_data_figure_10/replicated_gnomAD_category_REVEL.csv \
    --output data/intermediate/figures/extended_data_figure_10/gnomad_with_ancestry_af.csv
```

(`--ref-fasta`/`--exomes-dir`/`--genomes-dir` default to the
`data/input/gnomAD/` paths from step 1; pass them explicitly if your copies
live elsewhere.) See [`docs/map_gnomad_ancestry.md`](map_gnomad_ancestry.md)
for the full column reference and matching algorithm.

To update the figure's checked-in input data with these results:

```bash
cp data/intermediate/figures/extended_data_figure_10/vus_with_ancestry_af.csv data/input/gnomAD/
cp data/intermediate/figures/extended_data_figure_10/gnomad_with_ancestry_af.csv data/input/gnomAD/
```

## 4. Generate the figure

```bash
src/scripts/run_make_extended_data_figure_10.sh
```

See [`docs/figures.md`](figures.md#extended-data-figure-10-srcmake_extended_data_figure_10py)
for outputs and CLI flags.
