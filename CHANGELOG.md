# Changelog

## v2.0.1 — 2026-10-16

Small fixes and figure polish since v2.0.0. 7 commits.

- **Figure 2g**: fixed a VAMP-seq codon substitution count inflated by 2
  (29,725 instead of 29,723) — `len(set(...))` collapsed all NaN
  `amino_acid_change` values (from in-frame deletions with no single codon
  substitution) into one spurious "unique" entry per TSC2 domain subset;
  switched to `.nunique()`, which excludes NaN by default and matches the
  SGE side's existing counting method.
- **Figure 4 panel d**: dropped a spurious "±8" threshold line/label — the
  panel drew 5 evidence-strength threshold lines per side but only ±1
  through ±4 are meaningful; the 5th jumped straight to a "+8"/"-8" label
  with no such tier.
- **Biobank OR plots (Extended Data Figure 2, Figure 6b)**: significance
  indicators now use Holm-adjusted p-values; added a pathogenic/benign
  divider; the VEP plot now uses missense groups; Figure 6b reverted to
  fill-mapped (white/grey/black) significance points instead of varying
  point shape, per feedback that the color-coded circles read better; and
  its D3 "Cardiovascular" facet header is now hyphenated ("Cardio-
  vascular") to match D2 so it doesn't overflow when only KCNH2 remains
  under that category after filtering.
- **docs(figure-assembly)**: noted that Figure 2's panels are embedded
  directly in `figures_2-6.ai` rather than linked, since Illustrator can't
  place or link SVG files; the assembled Illustrator documents were
  re-exported to pick up all of the above.

## v2.0.0 — 2026-10-13

Everything since `v2.0.0-beta.1`. Promotes the beta to a final release:
no breaking changes, but closes out several reproducibility gaps the beta
left open (externally-produced inputs that weren't checked in or
documented, a predictor/training-data provenance audit, and the manual
figure-assembly step), plus a handful of MaveDB-fetching correctness
fixes. 24 commits.

### Reproducibility: externally-produced inputs now documented and checked in

- **ExCALIBR calibration artifacts** (per-dataset JSON calibrations, the
  Figure 4 out-of-bag export, the gene-specific predictor calibration seed
  sheet, and ExCALIBR's own copy of the integrated dataset) are now checked
  in and documented, correcting docs that previously claimed some of this
  data had no source anywhere in the repo.
- **Extended Data Figure 10** gained a full, reproducible pipeline
  (`src/map_gnomad_ancestry.py`, `src/filter_pillar_data.py`) from raw
  gnomAD v4.1 data through to the figure, replacing a checked-in snapshot
  with no documented way to regenerate it.
- **The biobank carrier odds-ratio workflow** is now documented end-to-end
  (`docs/biobank_analysis.md`): this repo builds the input, a manual
  All of Us Researcher Workbench hand-off computes the odds ratios, and the
  result is checked back in. README's three-stage overview now shows both
  this and the ExCALIBR hand-off as explicit "Interim" steps instead of
  oddly-placed sub-items of Stages 1 and 2.
- **Predictor and predictor-training-data provenance** was audited:
  documented where AlphaMissense, REVEL, and the MutPred2 training-variant
  source actually come from (including the REVEL reformatting step the
  pipeline needs that wasn't documented anywhere), consolidated into a new
  README "Preparation" section with working commands for every external
  file the pipeline needs before a run. `build_training_variant_files.py`
  now auto-detects and accepts MutPred2's publicly-downloadable
  training-variant format directly. Added `.gitignore` entries for the
  pinned AlphaMissense/REVEL files so they can't be committed by accident,
  and removed stale documentation referencing raw MutPred2/DDX3X source
  files that had since been deleted.

### MaveDB data-fetching fixes

- Synced the score-set URN list with `score_sets.tsv` (a stale TSC2_IGVF
  URN was causing 404s; a dropped LDLR dataset's URN and a duplicate
  BARD1_IGVF entry were cleaned up).
- `fetch_mavedb_scores.py` now tolerates a 404 from MaveDB's
  `/mapped-variants` endpoint instead of crashing, and its default output
  path matches what the pipeline actually reads with no manual rename.
- Documented the full validation of the API-based fetch route against the
  SQL/DB version's output, and why each remaining difference doesn't
  affect pipeline results.

### Figure assembly (manual, Adobe Illustrator)

- Documented and checked in the manual figure-assembly step for the first
  time (`docs/figure_assembly.md`): the assembled Illustrator documents,
  externally-sourced figure assets, and two reusable ExtendScript utilities
  (`export_artboards.jsx`, `relink_extended_data_figures.jsx`).
- Figure 4 gained a calibrated PDF export and was folded into the main
  Illustrator assembly alongside Figures 2, 3, 5, and 6.
- Extended Data Figure 2's panel a now facets by assay instead of gene
  (so genes with multiple assays get one panel each), added
  CHEK2/KCNE1/LDLR/PTEN datasets, and switched to `cairo_pdf` to fix
  silently-substituted tick-label glyphs.

### Other

- Added `CITATION.cff` and a script to submit updated Supplementary Data 1
  to the IGVF data portal.
- Corrected README's MAVE dataset count (91 → 90), Zenodo citation, and
  third-party data license/attribution coverage.

## v2.0.0-beta.1 — 2026-10-12

Everything since the `v1.0` tag (the state behind the original bioRxiv
preprint, submitted 2026-02-14). This release accompanies the paper
resubmission and represents a near-total rebuild of the pipeline: the
single-notebook workflow behind v1 has been replaced with a three-stage,
Dockerized pipeline, the classification analysis has grown a calibration and
ablation framework that didn't exist in v1, and most manuscript figures have
been rebuilt to read directly from that pipeline's output instead of static
caches. 252 commits.

### Breaking changes

- **Extended Data Figure numbering has changed** — figures were renumbered
  more than once as panels were added, dropped, or merged (most recently to
  close the gap left by dropped panels). Any reference to an Extended Data
  Figure number from the preprint should be re-checked against the current
  manuscript draft.
- **Pipeline steps were renumbered** as new steps were inserted (e.g. the old
  Step 13 is now Step 16). Scripts/docs referencing step numbers from v1 are
  stale.
- **The archived `Integrated_variant_effect_dataset_pipeline.ipynb` notebook
  has been removed.** Stage 1 (variant annotation) now runs exclusively
  through the Dockerized `scripts/variant_annotation_pipeline.sh`, built on
  the sibling `variant-annotation` project (see below).
- Tracked data paths were reconciled from mixed `Data/`/`data/` casing onto
  lowercase `data/` throughout, fixing a macOS/Linux case-sensitivity
  mismatch that broke on Linux CI and Docker.
- Figure/output directory layout changed: figure outputs now live under
  `notebooks/figures/<figure>/assets/` (previously scattered), and several
  figure scripts/notebooks/directories were renamed to match the new
  numbering.

### Pipeline & infrastructure

- **All three stages now run in Docker.** Stage 1 (variant annotation) was
  Dockerized first; Stage 2 (analysis notebooks) and Stage 3 (Python figures)
  followed, plus a dedicated `r-figures` Docker service for the R-based
  figure scripts. Poetry remains for interactive notebook/script development.
- Stage 1 is now built on the sibling **variant-annotation** project,
  vendored as a git submodule (`vendor/variant-annotation`) and driven by
  `scripts/run_variant_annotation_pipeline.sh`, which interleaves this
  project's own steps between `variant-annotation`'s numbered ones. Bumped
  repeatedly this cycle for upstream fixes (ClinVar download retries,
  identity-allele trust, ClinGen/identity-HGVS passthrough, identity-delins
  ref/alt, intron-spanning deletion ref handling).
- Adopted Poetry + Ruff tooling (PEP 621 `pyproject.toml`, pinned Python
  3.11+/3.12 dev env), a `pre-commit` hook that strips notebook outputs via
  `nbstripout`, and a `click`-based CLI convention for standalone scripts.
- Added smoke tests (`scripts/smoke_test_variant_annotation.sh`,
  `scripts/smoke_test_analysis.sh`) that run Stage 1 and Stage 2 end-to-end
  against small real-data fixtures, so environment setup can be verified
  without a multi-hour full run.
- Added a **BSD 3-Clause license**, crediting Jeremy Stone, Malvika Tejura,
  Ross Stewart, Yuriy Sverchkov, Ivan Woo, and Haneen Abderrazzaq.
- README rewritten as a complete, non-interactive-commands pipeline overview;
  `docs/` grew substantially (variant-annotation pipeline architecture, smoke
  tests, figure regeneration, data provenance, and several methodology notes
  — splice-aware filtering motivation, dedup strategy investigations, OddsPath
  calibration comparisons, ClinGen/AlphaMissense discordance).

### New pipeline steps and tooling

- `postprocess_mavedb_functional_classifications`,
  `add_mavedb_active_calibration_columns`, `recalculate_clingen_classification`,
  `flag_variants`, `build_training_variant_files`, and `derive_score_set_urn`
  joined the Dockerized pipeline between `variant-annotation`'s numbered
  steps.
- `annotate_simplified_consequence` (condensed consequence, splice-variant,
  splice-amino-acid columns) migrated in-project.
- `build_variant_reclassification_dataset.py`: builds a deduplicated,
  DNA-level variant export for reclassification analysis, with ExCALIBR- and
  OddsPath-scoring variants and a `--predictor` option (F9/TP53 use cases
  documented in `docs/`).
- `add_excalibr_columns`: merges ExCALIBR point/evidence columns into the
  integrated dataset.
- Vectorized the variant-classification row-wise scoring step and split out
  a slow Excel-write cell for performance.

### Classification & calibration analysis

- Added ExCALIBR and OddsPath calibration loaders (`load_excalibr_calibrations`,
  `load_oddspath_calibrations`), supporting multiple ExCALIBR score ranges
  per evidence-strength point.
- Supplementary Data 6 now defaults to OddsPath-calibrated datasets.
- Closed out numerous correctness gaps in the OddsPath/gene-specific controls
  dedup logic (genomic-coordinate casting, remaining controls gaps,
  CHEK2/LDLR/TSC2 fixes below).
- Added a VUS/gnomAD/Unobserved dedup-strategy parameter and adopted the
  decided default after investigation (see `docs/` for the sort-order
  comparisons behind that decision).
- Added a **functional-vs-predictor evidence ablation analysis**
  (Extended Data Figure 10): class-upgrade accounting, ClinVar discordance
  computation for AlphaMissense/MutPred2, PASS-only exome filtering, and a
  finalized calibrated by-scope figure.
- `mave_dataset_stats` grew from a single report into the project's main
  stats surface: filter-funnel, ClinGen, and control-concordance sections;
  chi-squared tests for gnomAD/VUS/Unobserved PLP and BLB rates;
  missense-only and ClinVar-2025 variants of each; per-predictor ablation
  synergy tables (split into PLP/BLB); VUS-missense and VAMP-seq coverage;
  registered-test and ablation-sufficiency stats; a paired Markdown report
  output with manuscript-cited figures highlighted.

### Data updates

- Retained splice-aware SGE assay variants that were previously caught by a
  blanket splice-variant exclusion (PR #4); documented the motivation in
  `docs/`.
- LDLR: excluded LA module 1 (aa 25–65) and all six LA modules from
  calibration; prioritized the +VLDL uptake assay for modules 2 and 6.
- TSC2_IGVF updated to MaveDB v3; dropped `LDLR_Tabet_2025_presence_VLDL`.
- CHEK2 QC workbook fixes (transcript-prefix stripping for flag matching).
- REVEL scores are now matched by amino-acid change rather than genomic
  coordinate, and REVEL-training variants are excluded/exposed as points in
  the reclassification export.
- Refreshed the gnomAD/All of Us (AoU) biobank export; the AoU-only input is
  now the stand-in for the (gitignored, separately downloaded) biobank data.
- Final pipeline output is now gzipped (`data/output/maves/`).

### Figures

- **Figure 2i**: added "Functionally Abnormal (Truncating)" and
  missense-specific rows (the latter behind a flag), dropped TSC2.
- **Figure 3**: now fully reproducible from checked-in/rebuildable inputs;
  added a calibrated bar chart, Euler diagram (deduped consistently by
  genomic variant), configurable variant-counting scope, AA-onto-NT
  projection by default, and bubble-plot legibility fixes (label size,
  callout crowding, mixed-provenance gene distinction).
- **Figure 4**: cache is now rebuilt from the live pipeline instead of a
  static file; added REVEL/GeneSpecific panels sourced from the
  `fig4_data_oob` export, fixed histogram bin alignment, gene-title
  bold-italic, and panel-c statistics.
- **Figures 5/6**: added calibrated confusion-matrix, Sankey, and
  real-print-size donut charts; restricted the Fig 5c/5f metrics bar charts
  to sensitivity/specificity; added a missense-only render mode and optional
  consequence filtering.
- **Extended Data Figures**: added/finalized calibrated donut, Sankey,
  confusion-matrix, rank-plot, and heatmap variants across Figures 4, 6, 7,
  8, and 9; added a missense-only Extended Data Figure 4 and an "alt" variant
  of Figures 6 and 7; added Extended Data Figure 5 (ClinVar discordance for
  AlphaMissense/MutPred2); finalized Extended Data Figure 10 as the
  calibrated ablation-by-scope figure; added a new, Dockerized gnomAD
  ancestry resolution-rate figure (Extended Data Figure 16, with a calibrated
  6.5in print version), with PASS-only exome filtering and uppercase
  ancestry labels by default.
- Numerous figure-correctness fixes: Arial font enforcement, gene-symbol
  facet-strip italicization, legend/axis/panel sizing, raw
  `Total_Points_REVEL` plotted instead of a zero-split copy, stray
  `Rplots.pdf`/untracked HTML output cleanup, and missing output-directory
  creation before writes.
