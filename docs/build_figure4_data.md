# build_figure4_data

Rebuilds `notebooks/figures/figure_4/figure4_data.json.gz`, the cache
`notebooks/figures/figure_4/figure4.ipynb` loads before calling `plot_figure4()`.
Nothing in the repo previously generated that cache — it was committed
pre-built. This script reconstructs the parts a current pipeline run can
support, and carries the rest forward unchanged from a prior cache.

## What gets rebuilt from scratch

| Fields | Source |
|---|---|
| `scoreset_2018`, `scoreset` | `Scoreset` objects for `MSH2_Jia_2021` (2018 vs. 2025 ClinVar release), built from `data/output/maves/integrated_variant_effect_dataset.tsv.gz` |
| `n_c`, `scoreset_flipped`, `indv_summary["prior"]`, `indv_summary["point_ranges"]` | the per-dataset exCALIBR calibration JSON at `data/input/mave_calibration/excalibr/json/MSH2_Jia_2021_clinvar_2018.json` |
| `gene_4e`, `dist_4e`, `labdat_4e`, `snvdf_4e`, `sorted_thresholds_4e`, `oldsorted_thresholds_4e` | MSH2 REVEL gene-specific-vs-genome-wide calibration, from the integrated dataset plus the `REVEL_gene_specific_calibration` sheet of `Supplementary_Data_4.xlsx` |
| `dist_4f`, `finalout_4f` | the same REVEL comparison across a fixed 6-gene set (`FIGURE_4F_GENES`) |

## What can't be rebuilt by this script (must be carried forward)

Pass `--cached-json` pointing at a prior `figure4_data.json.gz` (e.g.
`notebooks/figures/figure_4/old_figure4_data.json.gz`) to carry these fields
forward unchanged -- **this script** has no other source for them:

- **Panel a's density-band fit** — `fits`, and `indv_summary`'s
  `priors`/`log_lr_plus`/`score_range`/`C`, plus the top-level `score_range`.
  This needs exCALIBR's raw per-bootstrap fit output (1000 bootstrap replicate
  fits), which this script doesn't read from anywhere.
- **Panel c's out-of-bag confusion matrices** — `danzs_oob`, `auths_oob`,
  `datasets`. Needs exCALIBR's out-of-bag evidence-direction predictions per
  MAVE dataset, which this script also doesn't read from anywhere.
- **Panel d's cartoon** — `prior`, `Post_p`, `Post_b`, `p_data`, `b_data`.
  This panel is illustrative: the posterior curves it actually plots
  (`1 - x**2` and `exp(3x)/exp(3)`) are hardcoded in
  `notebooks/figures/figure_4/plot_utils.py`'s `plot_panel_d`, not derived from
  these inputs. There's no recorded recipe for the original synthetic
  `p_data`/`b_data`/`Post_p`/`Post_b` values -- this one is genuinely
  unrecoverable, not just unsupported by this script.

Running without `--cached-json` raises immediately and lists these fields.

**Note:** `figure4.ipynb` itself no longer actually needs this script's cache
for panels a, b, c, e, or f -- it reads
`data/input/mave_calibration/excalibr/fig4_data_oob.tar.gz` (a full
out-of-bag export from the same external ExCALIBR calibration run as the
JSON calibration files above; see
`docs/load_excalibr_calibrations.md#where-the-calibration-json-files-come-from`)
directly, which carries the raw bootstrap fits (`bootstrap_fits.json.gz`)
and the out-of-bag confusion matrices (`panel_c.json`'s `danz_agg`/`auth_agg`)
that this script can't reconstruct, plus the MSH2 scoreset and REVEL data
this script *can*. Only panel d's cartoon and part of panel e's legacy
genome-wide REVEL thresholds (`genome_wide_thresholds_4e`) still come from
`--cached-json`'s blob -- see `figure4.ipynb`'s first cell. This script's
own output (`figure4_data.json.gz`) is therefore mostly superseded for this
notebook already; it's kept for the fields `fig4_data_oob.tar.gz` doesn't
carry, and as a `--cached-json` source for anyone who wants to reconstruct
`danzs_oob`/`auths_oob` into *this script's* schema rather than reading the
tarball directly.

## Known discrepancies vs. the original figure

- **Schema drift**: the integrated dataset no longer has a plain `ID` column
  (`Scoreset`/`Variant` group rows by it); it now ships
  `mavedb_variant_urn` instead. This script renames the column before
  constructing `Scoreset` objects. Verified against the previous cache: panel
  b's MSH2 ClinVar BLB count goes from 230 (published figure) to 229 today —
  the residual difference is ClinVar's ongoing updates since the figure was
  made, not a reconstruction bug.
- **`FIGURE_4F_GENES`** (`BRCA1`, `BRCA2`, `F9`, `MSH2`, `TP53`, `TSC2`) is the
  fixed gene list the *original* figure compared in panel f. It isn't
  recoverable from any single data-driven rule — e.g. "genes with a
  fully-populated gene-specific REVEL calibration" matches 30 genes, not these
  6, and doesn't even contain all 6 of these. Treat it as an editorial choice;
  confirm with the paper authors before changing it.
- **`finalout_4f` row counts don't match the original 1:1** (e.g. today's
  MSH2 subset is 6,184 rows vs. ~10,535 in the original). Panel e's
  single-gene REVEL data (`snvdf_4e`, deduplicated by genomic position) does
  match the original almost exactly, so the REVEL thresholds and methodology
  here are verified correct — but the original `finalout_4f` was evidently
  built at a different row granularity (e.g. including both nucleotide- and
  amino-acid-level variant representations) that could not be reverse
  engineered from the current, already-deduplicated integrated dataset. This
  script deliberately does not try to inflate the row count to match; the
  category proportions and threshold assignments are correct.

## Usage

```bash
src/scripts/run_build_figure4_data.sh \
    --cached-json notebooks/figures/figure_4/old_figure4_data.json.gz
```

Optional flags: `--integrated-dataset`, `--excalibr-json-dir`,
`--supplementary-data-4`, `--output` (defaults to
`notebooks/figures/figure_4/figure4_data.json.gz`).
