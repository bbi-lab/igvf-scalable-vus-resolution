#!/usr/bin/env python3
"""Build the biobank-input export: one row per surviving variant-effect
measurement from `Variant_Classification_analysis.ipynb`'s exclusion rules,
not deduplicated to one row per DNA variant.

Reads the notebook's own intermediate checkpoint
(`data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`,
written at its cell 69, *before* category split), which already has the
LDLR LA-module-1 exclusion and the F9/TP53 restricted-dataset filter baked
in. Three further exclusions the notebook applies later, downstream of that
checkpoint, are re-applied here:

- `SFPQ` is dropped entirely (insufficient ClinVar controls).
- The CHEK2 QC flag (`CHEK2_Gebbia_2024.xlsx`, `Filter_CI == 1`, merged on
  (`hgvs_p` with its transcript prefix stripped, `auth_reported_score`) /
  (`hgvs_pro`, `score`)) is applied and then filtered out, along with every
  other `Flag == '*'` row. The prefix-stripping is a deliberate deviation
  from the notebook's own merge, which joins unstripped `hgvs_p` (e.g.
  `"NP_009125.1:p.Ser2Ala"`) directly against `hgvs_pro` (e.g.
  `"p.Ser2Ala"`, no prefix) -- confirmed against real data to never
  actually match, an apparently unintentional no-op in production that's
  out of scope to fix there.
- Rows tagged `splice_variant_not_measured` or `start_lost_variant_not_measured`
  in `VariantNotes`, or `splice_var_amino == 'Yes'`, are dropped -- unless
  `splice_measure == 'Yes'` (the dataset is curated as able to detect
  splicing effects), in which case the row is kept regardless of its splice
  flags. Rows tagged bare `conflicting_fxn_data` are dropped only when
  `dedup=True` (see below) -- with no deduplication, each row is one
  dataset's measurement of a variant rather than a single per-variant
  record, so more than one dataset disagreeing on a variant's effect is
  expected and kept, not treated as disqualifying. Deduplicating rows that
  are known to conflict would instead silently pick one dataset's value as
  the winner, so `dedup=True` excludes them outright rather than doing that.
- Rows with `<predictor>_train_amino == 'Yes'` are dropped -- variants used
  to train the chosen predictor, matching the notebook's own
  predictor-specific category sheets (e.g. `VUS_REVEL`, `Unobserved_REVEL`,
  `gnomAD_REVEL`, `controls_REVEL_GeneSpecific` for REVEL), which exclude
  them to avoid circularity with this file's predictor-based
  `Combined_points` (see `--predictor` below). `AM` has no such column in
  the checkpoint -- AlphaMissense isn't trained on this dataset's ClinVar
  controls the way REVEL/MutPred2 are -- so no row is excluded on this
  basis when `--predictor AM` is selected.

Six points columns are added:

- `ExCALIBR_points`: the literal ExCALIBR score-interval calibration value
  (`ExC_points_2018`, falling back to `ExC_points_2025`, for `BRCA1`/
  `PTEN`/`MSH2`/`TP53`; `ExC_points_2025` for every other gene) --
  independent of whether ExCALIBR is what's actually used for a given
  gene's evidence. TP53's on-file ExCALIBR calibrations are all
  2018-vintage, same reason as BRCA1/PTEN/MSH2 -- it has no effect on the
  official checkpoint today (the only TP53 dataset there with any ExCALIBR
  calibration, `TP53_Funk_2025`, was moved into the "not the meta analysis"
  exclusion list below alongside its true F9/TP53 assay peers), but matters
  for checkpoints that restore those per-assay TP53 datasets (see below).
- `OddsPath_points`: `OP_points` verbatim.
- `Functional_points`: `Fxn_points` verbatim -- the pipeline's own choice of
  `ExCALIBR_points` or `OddsPath_points` per gene (`F9`/`TP53` use
  `OddsPath_points`; every other gene uses `ExCALIBR_points`).
- `<predictor>_points` (named `REVEL_points` by default -- see `--predictor`
  below): `Points_<predictor>_GeneSpecific_GenomeWide` verbatim
  (gene-specific points, falling back to genome-wide; `--predictor-scope
  genome-wide` uses `Points_<predictor>_GenomeWide` instead, with no
  gene-specific contribution at all).
- `Conflict_<predictor>_GeneSpecific` or `Conflict_<predictor>_GenomeWide`
  (named `Conflict_REVEL_GeneSpecific` by default), matching
  `--predictor-scope`: `Conflicting_<predictor>_GeneSpecific` or
  `Conflicting_<predictor>_GenomeWide` verbatim -- the notebook's
  `split_zero`-derived column: `"Conflicting evidence"` when
  `Functional_points` and `<predictor>_points` disagree in sign, `"No
  evidence"` when both are missing, otherwise `Combined_points`'s value
  restated.
- `Combined_points`: `Functional_points + <predictor>_points` -- matches
  the notebook's own `Total_Points_GeneSpecific_REVEL` by default.

`--predictor` (`REVEL` by default; `AM`/`MP2` also available) selects which
computational predictor's points feed `<predictor>_points`/
`Conflict_<predictor>_*`/`Combined_points` -- the checkpoint carries the
same shape of columns (`Points_*_GeneSpecific_GenomeWide`,
`Points_*_GenomeWide`, `Conflicting_*_GeneSpecific`,
`Conflicting_*_GenomeWide`) for all three. Changing it can change
`--dedup`'s deduplication winner for a variant scored by more than one
dataset, since `dedup_by_max_abs_points` picks the candidate with the
greatest `abs(Combined_points)`, and `Combined_points` depends on the
chosen predictor: on the official checkpoint, switching `REVEL` -> `AM`
changes the winning dataset for 1,634 of 635,028 variants; `REVEL` -> `MP2`
changes it for 10,117.

By default (`dedup=False`), no deduplication is applied: a variant scored by
more than one dataset/assay keeps one row per measurement (this is the
difference from the pipeline's other per-category exports, which collapse to
one row per DNA variant). Passing `dedup=True` (`--dedup` on the CLI)
collapses to one row per DNA variant via
`src.lib.dedup.dedup_by_max_abs_points` (keyed on `Gene`/`Chrom`/
`hg38_start`/`ref_allele`/`alt_allele`: the candidate with the greatest
`abs(Combined_points)` wins, ties broken by `Dataset` name), and additionally
excludes bare `conflicting_fxn_data` rows -- see above.

Output columns are `integrated_variant_effect_dataset.tsv`'s full schema,
in its column order, with the six new points columns appended at the end.

To score F9/TP53 off ExCALIBR instead of OddsPath (so `Functional_points`
is ExCALIBR-based for every gene), pass a *checkpoint_file* built by
`notebooks/analysis/Variant_Classification_analysis_ExCALIBR_all_genes.ipynb`
(a copy of the main notebook that stops overriding F9/TP53's `Fxn_points`
with `OP_points`, and restores the per-assay F9/TP53 datasets the main
notebook drops -- `TP53_Fayer_2021_meta` and `F9_Popp_2025_model`, the
datasets the main notebook keeps, have no ExCALIBR calibration on file at
all).
"""

from pathlib import Path

import click
import pandas as pd

from src.lib.dedup import GENOMIC_KEY_COLS, dedup_by_max_abs_points

DEFAULT_CHECKPOINT_FILE = Path("data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz")
DEFAULT_CHEK2_FILE = Path("data/input/maves/CHEK2_Gebbia_2024.xlsx")
DEFAULT_OUTPUT_FILE = Path("data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz")

# ExC_points vintage override for ExCALIBR_points: these genes' on-file
# ExCALIBR calibrations are all 2018-vintage, so ExCALIBR_points prefers
# ExC_points_2018 over ExC_points_2025 for them (every other gene uses
# ExC_points_2025 only). TP53 is included even though its Functional_points
# come entirely from OddsPath on the official checkpoint -- ExCALIBR_points
# is documented above as the literal calibration value regardless of what's
# actually used, and harmless to include here regardless: the only TP53
# dataset on the official checkpoint with any ExCALIBR calibration,
# TP53_Funk_2025, is excluded from it entirely (see
# Variant_Classification_analysis.ipynb's "not the meta analysis" cell).
EXCALIBR_VINTAGE_OVERRIDE_GENES = frozenset({"BRCA1", "PTEN", "MSH2", "TP53"})

DISALLOWED_VARIANT_NOTES = frozenset({
    "splice_variant_not_measured",
    "splice_variant_not_measured;conflicting_fxn_data",
    "start_lost_variant_not_measured",
})

PREDICTOR_CHOICES = ("REVEL", "AM", "MP2")
PREDICTOR_SCOPE_CHOICES = ("gene-specific", "genome-wide")

# Training-circularity exclusion column per predictor (see
# apply_notebook_exclusions) -- None where the checkpoint carries no such
# column (AM has no training-circularity concept here).
PREDICTOR_TRAIN_EXCLUSION_COLUMN = {
    "REVEL": "revel_train_amino",
    "AM": None,
    "MP2": "mp2_train_amino",
}


def points_output_columns(predictor: str, predictor_scope: str) -> list[str]:
    """The six points columns' names for `predictor`/`predictor_scope` --
    always `ExCALIBR_points`, `OddsPath_points`, `Functional_points`, and
    `Combined_points`, plus `<predictor>_points` and
    `Conflict_<predictor>_GeneSpecific`/`Conflict_<predictor>_GenomeWide`
    (matching `predictor_scope`). Reduces to the historical
    `REVEL_points`/`Conflict_REVEL_GeneSpecific` names for the default
    `predictor="REVEL"`, `predictor_scope="gene-specific"`.
    """
    conflict_scope_label = "GeneSpecific" if predictor_scope == "gene-specific" else "GenomeWide"
    return [
        "ExCALIBR_points",
        "OddsPath_points",
        "Functional_points",
        f"{predictor}_points",
        f"Conflict_{predictor}_{conflict_scope_label}",
        "Combined_points",
    ]


BASE_OUTPUT_COLUMNS = [
    "Dataset", "Gene", "HGNC_id", "mavedb_variant_urn", "Chrom", "Strand", "hg38_start", "hg38_end",
    "ref_allele", "alt_allele", "auth_transcript_id", "transcript_pos", "transcript_ref", "transcript_alt",
    "aa_pos", "aa_ref", "aa_alt", "hgvs_c", "hgvs_p", "consequence", "most_severe_mutational_consequence",
    "auth_reported_score", "rna_score", "auth_reported_func_class", "auth_reported_func_class_category",
    "splice_measure", "gnomad_MAF",
    "clinvar_sig_2026", "clinvar_star_2026", "clinvar_date_last_reviewed_2026",
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025",
    "clinvar_sig_2018", "clinvar_star_2018", "clinvar_date_last_reviewed_2018",
    "nucleotide_or_aa", "Ensembl Transcript ID", "RefSeq Transcript ID",
    "Interval 1 Name", "Interval 1 Range", "Interval 1 Class",
    "Interval 2 Name", "Interval 2 Range", "Interval 2 Class",
    "Interval 3 Name", "Interval 3 Range", "Interval 3 Class",
    "Interval 4 Name", "Interval 4 Range", "Interval 4 Class",
    "Interval 5 Name", "Interval 5 Range", "Interval 5 Class",
    "Interval 6 Name", "Interval 6 Range", "Interval 6 Class",
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL",
    "spliceAI_DP_AG", "spliceAI_DP_AL", "spliceAI_DP_DG", "spliceAI_DP_DL",
    "ClinVar Variation Id_ClinGen_repo", "Allele Registry Id_ClinGen_repo", "Disease_ClinGen_repo",
    "Mondo Id_ClinGen_repo", "Mode of Inheritance_ClinGen_repo", "Assertion_ClinGen_repo",
    "Applied Evidence Codes (Met)_ClinGen_repo", "Applied Evidence Codes (Not Met)_ClinGen_repo",
    "Summary of interpretation_ClinGen_repo", "PubMed Articles_ClinGen_repo", "Expert Panel_ClinGen_repo",
    "Guideline_ClinGen_repo", "Approval Date_ClinGen_repo", "Published Date_ClinGen_repo",
    "Retracted_ClinGen_repo", "Evidence Repo Link_ClinGen_repo", "Uuid_ClinGen_repo",
    "Updated_Classification_ClinGen_repo", "Updated_Evidence Codes_ClinGen_repo",
    "REVEL", "REVEL_train", "AM_score", "AM_class", "MutPred2", "MP2_train",
    "simplified_consequence", "condensed_consequence", "splice_variant", "splice_var_amino", "Flag",
]

OUTPUT_COLUMNS = BASE_OUTPUT_COLUMNS + points_output_columns("REVEL", "gene-specific")


def apply_notebook_exclusions(
    df: pd.DataFrame, chek2_file: Path, dedup: bool = False, predictor: str = "REVEL"
) -> pd.DataFrame:
    """Re-apply the checkpoint-downstream exclusions from
    `Variant_Classification_analysis.ipynb` cells 71-77 and 95: `SFPQ`, the
    CHEK2 QC flag, unmeasured-splice/start-lost `VariantNotes` tags (bare
    `conflicting_fxn_data` is additionally excluded when `dedup=True` -- see
    module docstring), any other `Flag == '*'` row, and `predictor`-training
    variants (via `PREDICTOR_TRAIN_EXCLUSION_COLUMN`; no-op for `AM`).

    The CHEK2 merge key is normalized (`hgvs_p`'s transcript prefix, e.g.
    `"NP_009125.1:"`, stripped before matching against `hgvs_pro`) rather
    than reproducing the notebook's own merge verbatim: the notebook joins
    unstripped `hgvs_p` directly against `hgvs_pro`, which never actually
    matches (`hgvs_pro` carries no transcript prefix) -- a pre-existing,
    apparently unintentional no-op in production, flagged separately rather
    than fixed there.
    """
    df = df[df["Gene"] != "SFPQ"].copy()

    chek2 = pd.read_excel(chek2_file, header=0)
    df["_hgvs_p_no_transcript"] = df["hgvs_p"].str.replace(r"^[^:]+:", "", regex=True)
    df = df.merge(
        chek2[["hgvs_pro", "score", "Filter_CI"]],
        left_on=["_hgvs_p_no_transcript", "auth_reported_score"],
        right_on=["hgvs_pro", "score"],
        how="left",
    )
    df["Flag"] = df["Flag"].where(df["Filter_CI"] != 1, "*")
    df = df.drop(columns=["_hgvs_p_no_transcript", "hgvs_pro", "score", "Filter_CI"])

    disallowed_variant_notes = DISALLOWED_VARIANT_NOTES | ({"conflicting_fxn_data"} if dedup else set())
    df = df[
        ~df["VariantNotes"].isin(disallowed_variant_notes)
        & ((df["splice_var_amino"] != "Yes") | (df["splice_measure"] == "Yes"))
    ]
    df = df[df["Flag"] != "*"]
    train_exclusion_col = PREDICTOR_TRAIN_EXCLUSION_COLUMN[predictor]
    if train_exclusion_col is not None:
        df = df[df[train_exclusion_col] != "Yes"]
    return df


def add_points_columns(df: pd.DataFrame, predictor: str = "REVEL", predictor_scope: str = "gene-specific") -> pd.DataFrame:
    """Add `ExCALIBR_points`, `OddsPath_points`, `Functional_points`,
    `<predictor>_points`, `Conflict_<predictor>_GeneSpecific`/
    `Conflict_<predictor>_GenomeWide`, and `Combined_points` (see
    `points_output_columns` for the exact names) -- see module docstring
    for the exact definitions.
    """
    df = df.copy()

    df["ExCALIBR_points"] = df["ExC_points_2025"]
    vintage_mask = df["Gene"].isin(EXCALIBR_VINTAGE_OVERRIDE_GENES)
    excalibr_2018 = pd.to_numeric(df.loc[vintage_mask, "ExC_points_2018"], errors="coerce")
    excalibr_2025 = pd.to_numeric(df.loc[vintage_mask, "ExC_points_2025"], errors="coerce")
    df.loc[vintage_mask, "ExCALIBR_points"] = excalibr_2018.fillna(excalibr_2025)

    df["OddsPath_points"] = df["OP_points"]
    df["Functional_points"] = df["Fxn_points"]

    if predictor_scope == "gene-specific":
        points_source = f"Points_{predictor}_GeneSpecific_GenomeWide"
        conflict_source = f"Conflicting_{predictor}_GeneSpecific"
    else:
        points_source = f"Points_{predictor}_GenomeWide"
        conflict_source = f"Conflicting_{predictor}_GenomeWide"
    points_col, conflict_col = points_output_columns(predictor, predictor_scope)[3:5]

    df[points_col] = df[points_source]
    df[conflict_col] = df[conflict_source]
    df["Combined_points"] = (
        pd.to_numeric(df["Functional_points"], errors="coerce").fillna(0)
        + pd.to_numeric(df[points_col], errors="coerce").fillna(0)
    )
    return df


def build_reclassification_dataset(
    checkpoint_file: Path,
    chek2_file: Path,
    dedup: bool = False,
    predictor: str = "REVEL",
    predictor_scope: str = "gene-specific",
) -> pd.DataFrame:
    df = pd.read_csv(checkpoint_file)
    df = apply_notebook_exclusions(df, chek2_file, dedup=dedup, predictor=predictor)
    df = add_points_columns(df, predictor=predictor, predictor_scope=predictor_scope)
    if dedup:
        df = dedup_by_max_abs_points(df, points_col="Combined_points", genomic_key_cols=GENOMIC_KEY_COLS)
    return df[BASE_OUTPUT_COLUMNS + points_output_columns(predictor, predictor_scope)]


@click.command(help=__doc__)
@click.argument(
    "checkpoint_file",
    required=False,
    default=DEFAULT_CHECKPOINT_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--chek2-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CHEK2_FILE,
    help=f"Path to the CHEK2 QC workbook (default {DEFAULT_CHEK2_FILE}).",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_FILE,
    help=f"Output path (default {DEFAULT_OUTPUT_FILE}), written as gzip-compressed TSV.",
)
@click.option(
    "--dedup",
    is_flag=True,
    default=False,
    help=(
        "Collapse to one row per DNA variant (greatest abs(Combined_points) wins) and additionally "
        "exclude bare conflicting_fxn_data rows. Default is one row per surviving measurement, for "
        "biobank-input use."
    ),
)
@click.option(
    "--predictor",
    type=click.Choice(PREDICTOR_CHOICES),
    default="REVEL",
    show_default=True,
    help="Computational predictor to combine with Functional_points into Combined_points.",
)
@click.option(
    "--predictor-scope",
    type=click.Choice(PREDICTOR_SCOPE_CHOICES),
    default="gene-specific",
    show_default=True,
    help=(
        "gene-specific falls back to genome-wide when no gene-specific calibration exists "
        "(Points_<predictor>_GeneSpecific_GenomeWide); genome-wide ignores gene-specific "
        "calibration entirely (Points_<predictor>_GenomeWide)."
    ),
)
def main(
    checkpoint_file: Path,
    chek2_file: Path,
    output: Path,
    dedup: bool,
    predictor: str,
    predictor_scope: str,
) -> None:
    result = build_reclassification_dataset(
        checkpoint_file, chek2_file, dedup=dedup, predictor=predictor, predictor_scope=predictor_scope
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False, compression="gzip")
    click.echo(f"Wrote {len(result)} rows to {output}")


if __name__ == "__main__":
    main()
