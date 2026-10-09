#!/usr/bin/env python3
"""Total and per-gene count of registered clinical genetic tests (NCBI GTR)
for this study's 41 genes.

Reads NCBI's Genetic Testing Registry bulk export
(`data/input/genes/test_condition_gene.txt.gz`, the same file
`src/build_figure3_data.py` uses for Figure 3a's `gene_test_count`) and counts
rows with `test_type == "Clinical"` and `object == "gene"`, restricted to this
study's gene symbols. The study's gene list comes from
`Supplementary_Data_3.xlsx`'s `Curation` sheet `Gene` column (41 unique
entries; `CALM1, CALM2, CALM3` is one combined entry there, since it's a
single gene target assayed across three paralogs -- its count is the sum of
the three symbols' individual GTR counts, reported under the combined label).

Unlike the sibling `igvf-coding-variants-portal` project's
`populate-derived-data.ts`, which scrapes NCBI's live GTR search page per gene
(and silently logs+continues on any HTTP/parse failure), this reads the
already-checked-in bulk export and raises rather than swallowing errors: a
missing input file, an unreadable/malformed bulk export, or a study gene
symbol absent from it are all treated as failures, not zero counts -- see
`load_gene_test_counts` and `count_registered_genetic_tests`.

Usage:
    python -m src.count_registered_genetic_tests \\
        [--curation-sheet PATH] [--testing-registry PATH]
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_CURATION_SHEET = Path("data/input/maves/Supplementary_Data_3.xlsx")
DEFAULT_TESTING_REGISTRY_PATH = Path("data/input/genes/test_condition_gene.txt.gz")

CURATION_SHEET_NAME = "Curation"
CURATION_GENE_COL = "Gene"

GTR_TEST_TYPE_COL = "test_type"
GTR_OBJECT_COL = "object"
GTR_GENE_SYMBOL_COL = "gene_symbol"
GTR_CLINICAL_TEST_TYPE = "Clinical"
GTR_GENE_OBJECT = "gene"


def load_study_genes(curation_sheet_path):
    """Return the study's 41 `Gene` entries (as they appear in the Curation
    sheet, e.g. "CALM1, CALM2, CALM3" as one combined entry) mapped to their
    individual gene symbols.

    Returns {curation entry: [symbol, ...]}. Raises ValueError if the
    Curation sheet or its `Gene` column is missing.
    """
    try:
        curation = pd.read_excel(curation_sheet_path, sheet_name=CURATION_SHEET_NAME)
    except (FileNotFoundError, ValueError) as error:
        raise ValueError(f"couldn't read {CURATION_SHEET_NAME!r} sheet from {curation_sheet_path}: {error}") from error

    if CURATION_GENE_COL not in curation.columns:
        raise ValueError(
            f"{curation_sheet_path} {CURATION_SHEET_NAME!r} sheet is missing a {CURATION_GENE_COL!r} column"
        )

    entries = sorted(curation[CURATION_GENE_COL].dropna().unique())
    return {entry: [symbol.strip() for symbol in entry.split(",") if symbol.strip()] for entry in entries}


def load_gene_test_counts(testing_registry_path):
    """Return a Series of clinical gene-level GTR test counts indexed by gene symbol.

    Raises ValueError if the bulk export is missing its expected columns.
    """
    testing_registry = pd.read_csv(testing_registry_path, sep="\t")

    required_cols = {GTR_TEST_TYPE_COL, GTR_OBJECT_COL, GTR_GENE_SYMBOL_COL}
    missing_cols = required_cols - set(testing_registry.columns)
    if missing_cols:
        raise ValueError(f"{testing_registry_path} is missing expected column(s): {sorted(missing_cols)}")

    testing_registry_gene = testing_registry[
        testing_registry[GTR_TEST_TYPE_COL].eq(GTR_CLINICAL_TEST_TYPE)
        & testing_registry[GTR_OBJECT_COL].eq(GTR_GENE_OBJECT)
    ]
    return testing_registry_gene.groupby(GTR_GENE_SYMBOL_COL).size()


def count_registered_genetic_tests(
    curation_sheet_path=DEFAULT_CURATION_SHEET, testing_registry_path=DEFAULT_TESTING_REGISTRY_PATH
):
    """Return (total_count, per_gene_counts) for this study's genes.

    `per_gene_counts` is a {Curation `Gene` entry: test count} dict, in
    Curation-sheet order, with combined entries (e.g. "CALM1, CALM2, CALM3")
    summed across their individual symbols. Raises ValueError if any study
    gene symbol has no matching rows in the testing registry bulk export --
    every one of this study's 41 genes is expected to have at least one
    registered clinical test, so a zero count signals a symbol mismatch
    (renamed/aliased gene) rather than a genuine absence of testing.
    """
    study_genes = load_study_genes(curation_sheet_path)
    test_counts = load_gene_test_counts(testing_registry_path)

    all_symbols = sorted({symbol for symbols in study_genes.values() for symbol in symbols})
    missing_symbols = [symbol for symbol in all_symbols if symbol not in test_counts.index]
    if missing_symbols:
        raise ValueError(
            f"no clinical gene-level test records found in {testing_registry_path} for gene symbol(s): "
            f"{missing_symbols} -- check for a renamed/aliased gene symbol"
        )

    per_gene_counts = {entry: int(test_counts.loc[symbols].sum()) for entry, symbols in study_genes.items()}
    total_count = sum(per_gene_counts.values())
    return total_count, per_gene_counts


def format_report(total_count, per_gene_counts):
    lines = [f"Total registered clinical genetic tests across {len(per_gene_counts)} genes: {total_count}", ""]
    for entry, count in per_gene_counts.items():
        lines.append(f"{entry}: {count}")
    return "\n".join(lines)


@click.command(help=__doc__)
@click.option(
    "--curation-sheet",
    "curation_sheet_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CURATION_SHEET,
    show_default=True,
    help="Supplementary_Data_3.xlsx (or equivalent), read for its Curation sheet's Gene column.",
)
@click.option(
    "--testing-registry",
    "testing_registry_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_TESTING_REGISTRY_PATH,
    show_default=True,
    help="NCBI GTR bulk export (test_condition_gene.txt.gz or equivalent).",
)
def main(curation_sheet_path, testing_registry_path):
    try:
        total_count, per_gene_counts = count_registered_genetic_tests(curation_sheet_path, testing_registry_path)
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    click.echo(format_report(total_count, per_gene_counts))


if __name__ == "__main__":
    main()
