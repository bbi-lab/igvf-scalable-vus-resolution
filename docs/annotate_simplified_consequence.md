# Annotate Simplified Consequence

`src/annotate_simplified_consequence.py` adds a `simplified_consequence`
column to a CVFG variants TSV, collapsing each DNA candidate's VEP
consequence term down to a single, coarser "SO summary term" (e.g.
`missense_variant` rather than `missense_variant^splice_region_variant`).

## Source logic

Ported from the `get_simplified_consequence` cell in
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`, which
built a `VEP output term -> SO summary term` mapping from a curated
Ensembl-consequence table (`extended_ensembl_consequence.csv.gz`, extracted
from Ensembl's [predicted data
page](https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html)
on 2025-03-18) -- picking, for any VEP term listed more than once, the row
with the highest `Importance` -- and applied it to a single flat
`consequence` column (comma-separated VEP terms, one row per variant),
taking the first comma-separated term with a mapping-table entry.

The variant-annotation pipeline's `annotate_vep.py` (step 10) instead
produces a **pipe-delimited** `vep.most_severe_mutational_consequence`
column, aligned one segment per DNA candidate (the same convention
`flag_variants.py`'s `Flag` column and `recalculate_clingen_classification.py`
use), with each segment already reduced to VEP's own single most-severe term
for that candidate. This script maps each segment independently through the
same mapping table and re-joins the results with `|`, rather than picking
among multiple comma-separated terms itself the way the notebook did.

A term with no mapping-table entry (e.g. `no_change`, a
`variant-annotation`-specific value not in Ensembl's own vocabulary) maps to
an empty string, same as the notebook's `NaN` case for an unmapped term.

## `condensed_consequence` and `splice_var_amino` are per-assayed-variant, not per-candidate

Besides `simplified_consequence`, this script adds three more columns.
`condensed_consequence` collapses `simplified_consequence`'s pipe-delimited
per-DNA-candidate values to one value: that value if every candidate agrees,
else `conflict_consequence` (`compute_condensed_consequence`). `splice_variant`
is a pipe-delimited, per-candidate Yes/No flag (from `simplified_consequence`
and the SpliceAI delta-score columns); `splice_var_amino` collapses it the
same way `condensed_consequence` collapses `simplified_consequence` -- blank
if no candidate has a value, that value if all non-blank candidates agree,
else `Yes` (`compute_splice_var_amino`).

`simplified_consequence` and `splice_variant` stay pipe-delimited, one
segment per DNA candidate, all the way through the rest of the pipeline.
`condensed_consequence` and `splice_var_amino` don't -- they're already
single values *per assayed variant* by the time this script runs. That
distinction matters once step 21 (`flatten_dna_variants.py`, see
`docs/variant_annotation_pipeline.md`'s "Step 21's explicit
`--dna-variant-columns`" section) explodes a multi-candidate row into one row
per DNA candidate: `simplified_consequence`/`splice_variant` are on that
step's explicit re-split list and get sliced down to each candidate's own
segment, but `condensed_consequence`/`splice_var_amino` aren't (they were
never pipe-delimited to begin with) and are simply copied unchanged onto
every exploded row of the group.

The practical effect shows up most for variants assayed at amino-acid
resolution (`nucleotide_or_aa == "aa"`): one measured functional score can map
to several different DNA-level candidates (different nucleotide changes
producing the same amino-acid substitution), and those candidates don't
always share a consequence. For example, TSC2's `urn:mavedb:00001267-0-3#56`
maps to one missense SNV (`NM_000548.5:c.4849G>C`) and three
splice-site-affecting delins candidates (`c.4849_4851delinsCCA`/`CCG`/`CCT`)
at the same codon. After the explode, each of those four candidates' own
`simplified_consequence` is a clean single value (`missense_variant` or
`splice_site_variant`), but all four rows' `condensed_consequence` still
reads `conflict_consequence`, inherited from the whole group.

**This is intentional, not a staleness bug.** `condensed_consequence` answers
an assayed-variant-group-level question -- "did every DNA-level
interpretation of this measurement agree?" -- not a per-candidate one, and
downstream code relies on exactly that: `src/ablation_variant_reclassification.py`'s
`is_missense_only`/`is_missense_or_start_loss` filters use strict
`condensed_consequence == "missense_variant"` (rather than
`simplified_consequence`) to build an unambiguous missense-only population
for predictor-performance comparisons (REVEL/AlphaMissense/MutPred2),
deliberately excluding any candidate whose assayed-variant group doesn't
unanimously agree -- see `docs/ablation_variant_reclassification.md`.
`splice_var_amino` is the mirror-image "did *any* candidate in this group
affect splicing?" flag, propagated across the same exploded rows the same
way.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.annotate_simplified_consequence \
  data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
```

Via Docker (same image as `flag_variants`/`recalculate_clingen_classification`,
see `compose.yaml`):

```bash
src/scripts/run_annotate_simplified_consequence.sh \
  data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
```

Like `run_flag_variants.sh`, this wrapper maps its input/output paths against
the `/work` staging mount (`${VARIANT_DATA_DIR:-./data}`), since it reads
from the same staged pipeline data those files live in rather than a
committed repo file. `--consequence-map-file` is remapped the same way
`run_flag_variants.sh` remaps `--filtering-dir` -- its default resolves
against this project's own bind mount (`/usr/src/app`) instead, since
`extended_ensembl_consequence.csv.gz` is a committed repo file, not staged
pipeline data.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--consequence-col` | `vep.most_severe_mutational_consequence` | Column with pipe-delimited VEP consequence term(s) per DNA candidate |
| `--consequence-map-file` | `data/input/reference/extended_ensembl_consequence.csv.gz` | CSV mapping `VEP output term` to `SO summary term`, with an `Importance` column for resolving duplicates |

Raises `click.ClickException` if `--consequence-col` isn't present in the
input file.

## Stdout report

After simplifying, the script prints per-candidate `simplified_consequence`
value counts (a candidate with an empty segment -- no VEP consequence, or an
unmapped term -- is reported as `(no consequence)`).
