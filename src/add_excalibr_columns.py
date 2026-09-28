#!/usr/bin/env python3
"""Add ExCALIBR point/evidence columns to the integrated variant effect dataset.

Manual pipeline step, run *after* the main variant-annotation pipeline
(`flag_variants`/`recalculate_clingen_classification` included) has produced
`data/output/maves/integrated_variant_effect_dataset.tsv.gz` -- the ExCALIBR
calibration output this script reads
(`data/input/maves/dataframe_with_points.csv.gz`) isn't available until
later, so this step can't be folded into that pipeline run.

`dataframe_with_points.csv.gz` carries every column of the integrated
dataset plus seven new `excalibr_*` columns (`excalibr_clinvar_release`,
`excalibr_points`, `excalibr_acmg_evidence_code`, `excalibr_lr_plus`,
`excalibr_prior`, `excalibr_posterior`, `excalibr_filter_reason`). This
script joins just those seven columns onto the integrated dataset by
(`Dataset`, `mavedb_variant_urn`, `hgvs_c`) -- the only combination of
shared columns that is unique in both files -- and writes the result to
`data/output/maves/integrated_variant_effect_dataset_with_excalibr.tsv.gz`.
A copy is then written to
`data/output/supplementary_data/Supplementary_Data_1.tsv.gz` (pass
`--skip-supplementary-data-copy` to suppress this).

`hgvs_c` is matched with one special case: variant-annotation's
map_variants.py fix for reference-identical alleles (see
`_reformat_confirmed_identity_hgvs` there) rewrites e.g.
`NM_022445.4:c.612C=` as the explicit `NM_022445.4:c.612delCinsC` --
ExCALIBR's own pipeline still emits the original single-anchor `=` form for
the same alleles. Both forms are normalized to the `=` form for key-matching
purposes only (see `_canonicalize_identity_hgvs_c`); the integrated
dataset's own `hgvs_c` value is kept, unmodified, in the output.
"""

import re
import shutil
from pathlib import Path

import click
import pandas as pd

DEFAULT_INTEGRATED_INPUT = Path("data/output/maves/integrated_variant_effect_dataset.tsv.gz")
DEFAULT_EXCALIBR_INPUT = Path("data/input/maves/dataframe_with_points.csv.gz")
DEFAULT_OUTPUT = Path("data/output/maves/integrated_variant_effect_dataset_with_excalibr.tsv.gz")
DEFAULT_SUPPLEMENTARY_DATA_COPY = Path("data/output/supplementary_data/Supplementary_Data_1.tsv.gz")

KEY_COLUMNS = ["Dataset", "mavedb_variant_urn", "hgvs_c"]
EXCALIBR_COLUMNS = [
    "excalibr_clinvar_release",
    "excalibr_points",
    "excalibr_acmg_evidence_code",
    "excalibr_lr_plus",
    "excalibr_prior",
    "excalibr_posterior",
    "excalibr_filter_reason",
]

# Matches the explicit del<bases>ins<bases> identity-allele hgvs_c that
# variant-annotation's map_variants.py emits for reference-identical alleles
# (ref == alt == bases -- see _reformat_confirmed_identity_hgvs there), e.g.
# "NM_022445.4:c.612delCinsC" or, for a multi-base identity allele,
# "NM_022445.4:c.593_595delCCTinsCCT".
_EXPLICIT_IDENTITY_HGVS_C_RE = re.compile(
    r"^(?P<prefix>.*:c\.)(?P<start>\d+)(?:_\d+)?del(?P<bases>[ACGTN]+)ins(?P<bases2>[ACGTN]+)$"
)


def _canonicalize_identity_hgvs_c(hgvs_c):
    """Canonicalize an explicit identity-form `hgvs_c` to the single-anchor `=` form.

    ExCALIBR's own pipeline hasn't picked up variant-annotation's identity-
    allele notation change (see module docstring), so its `hgvs_c` for these
    alleles is still e.g. `NM_022445.4:c.612C=` rather than
    `NM_022445.4:c.612delCinsC`. This makes the two forms compare equal for
    key-matching -- it does not touch any value that ends up in the output.

    Returns *hgvs_c* unchanged if it isn't the explicit identity form (e.g.
    it's already the `=` form, or it's a real, non-identity variant where
    the deleted/inserted bases differ).
    """
    m = _EXPLICIT_IDENTITY_HGVS_C_RE.match(hgvs_c)
    if not m or m.group("bases") != m.group("bases2"):
        return hgvs_c
    return f"{m.group('prefix')}{m.group('start')}{m.group('bases')}="


def add_excalibr_columns(integrated_df, excalibr_df):
    """Return `integrated_df` with `EXCALIBR_COLUMNS` joined in from `excalibr_df`.

    Joins on `KEY_COLUMNS`, with `hgvs_c` normalized via
    `_canonicalize_identity_hgvs_c` for matching purposes only -- the output
    keeps `integrated_df`'s own `hgvs_c` values. Raises ValueError if
    `KEY_COLUMNS` is missing from either input, if `EXCALIBR_COLUMNS` is
    missing from `excalibr_df`, if `KEY_COLUMNS` (with `hgvs_c` normalized)
    isn't unique in either input, or if the two inputs' normalized
    `KEY_COLUMNS` values don't match exactly (every row on each side must
    have a corresponding row on the other).
    """
    missing_integrated = [c for c in KEY_COLUMNS if c not in integrated_df.columns]
    missing_excalibr = [c for c in KEY_COLUMNS + EXCALIBR_COLUMNS if c not in excalibr_df.columns]
    if missing_integrated:
        raise ValueError(f"integrated dataset is missing required column(s): {', '.join(missing_integrated)}")
    if missing_excalibr:
        raise ValueError(f"ExCALIBR dataset is missing required column(s): {', '.join(missing_excalibr)}")

    join_key_columns = ["Dataset", "mavedb_variant_urn", "_hgvs_c_key"]
    integrated_join_df = integrated_df[["Dataset", "mavedb_variant_urn"]].assign(
        _hgvs_c_key=integrated_df["hgvs_c"].map(_canonicalize_identity_hgvs_c)
    )
    excalibr_join_df = excalibr_df[["Dataset", "mavedb_variant_urn"]].assign(
        _hgvs_c_key=excalibr_df["hgvs_c"].map(_canonicalize_identity_hgvs_c)
    )

    if integrated_join_df.duplicated().any():
        raise ValueError(f"integrated dataset has duplicate rows for key columns {KEY_COLUMNS}")
    if excalibr_join_df.duplicated().any():
        raise ValueError(f"ExCALIBR dataset has duplicate rows for key columns {KEY_COLUMNS}")

    integrated_keys = set(integrated_join_df.itertuples(index=False, name=None))
    excalibr_keys = set(excalibr_join_df.itertuples(index=False, name=None))
    if integrated_keys != excalibr_keys:
        only_integrated = len(integrated_keys - excalibr_keys)
        only_excalibr = len(excalibr_keys - integrated_keys)
        raise ValueError(
            f"key columns {KEY_COLUMNS} don't match between inputs: {only_integrated} key(s) only in the "
            f"integrated dataset, {only_excalibr} key(s) only in the ExCALIBR dataset"
        )

    excalibr_to_merge = excalibr_join_df.assign(
        **{col: excalibr_df[col] for col in EXCALIBR_COLUMNS}
    )
    return (
        integrated_df.assign(_hgvs_c_key=integrated_join_df["_hgvs_c_key"])
        .merge(excalibr_to_merge, on=join_key_columns, how="left")
        .drop(columns=["_hgvs_c_key"])
    )


@click.command(help=__doc__)
@click.option(
    "--integrated-input",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_INTEGRATED_INPUT,
    show_default=True,
    help="Integrated variant effect dataset TSV to add ExCALIBR columns to.",
)
@click.option(
    "--excalibr-input",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_EXCALIBR_INPUT,
    show_default=True,
    help="ExCALIBR CSV to pull the excalibr_* columns from.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT,
    show_default=True,
    help="Where to write the merged TSV.",
)
@click.option(
    "--supplementary-data-copy",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_SUPPLEMENTARY_DATA_COPY,
    show_default=True,
    help="Where to also copy the merged TSV, as Supplementary Data 1.",
)
@click.option(
    "--skip-supplementary-data-copy",
    is_flag=True,
    help="Don't write --supplementary-data-copy.",
)
def main(integrated_input, excalibr_input, output, supplementary_data_copy, skip_supplementary_data_copy):
    integrated_df = pd.read_csv(
        integrated_input, sep="\t", dtype=str, keep_default_na=False, engine="c", low_memory=False
    )
    excalibr_df = pd.read_csv(excalibr_input, dtype=str, keep_default_na=False, engine="c", low_memory=False)

    try:
        merged = add_excalibr_columns(integrated_df, excalibr_df)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {len(merged)} row(s) with {len(EXCALIBR_COLUMNS)} added ExCALIBR column(s) to {output}.")

    if not skip_supplementary_data_copy:
        supplementary_data_copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(output, supplementary_data_copy)
        click.echo(f"Copied to {supplementary_data_copy}.")


if __name__ == "__main__":
    main()
