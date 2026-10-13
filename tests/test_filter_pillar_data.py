import pandas as pd
import pytest
from click.testing import CliRunner

from src.filter_pillar_data import (
    add_revel_classification,
    find_matching_class,
    group_any,
    main,
    opposite_signs,
    parse_bound,
    points_from_intervals,
)
from src.filter_pillar_data import filter_pillar_data as filter_pillar_data_fn

DEFAULTS = {
    "excalibr_filter_reason": "",
    "aa_pos": "",
    "aa_ref": "",
    "aa_alt": "",
    "RefSeq Transcript ID": "NM_000001.1",
    "spliceAI_DS_AG": 0.0,
    "spliceAI_DS_AL": 0.0,
    "spliceAI_DS_DG": 0.0,
    "spliceAI_DS_DL": 0.0,
    "simplified_consequence": "missense_variant",
    "splice_measure": "No",
    "REVEL_train": "False",
    "mavedb_variant_urn": "urn1",
    "clinvar_sig_2025": "",
    "gnomad_MAF": "",
    "Flag": "",
    "REVEL": "0.5",
    "auth_reported_score": "1",
    "auth_reported_func_class_category": "",
}


def _points_df(rows):
    """rows: list of dicts overriding DEFAULTS; required keys must be given explicitly."""
    out = []
    for row in rows:
        merged = dict(DEFAULTS)
        merged.update(row)
        out.append(merged)
    return pd.DataFrame(out)


# ---------------- pure helper functions ----------------
def test_group_any_respects_required_column():
    df = pd.DataFrame({"k": ["a", "a", "b"], "req": [1.0, None, 1.0]})
    values = pd.Series([True, False, False])
    out = group_any(df, ["k"], values, required=["req"])
    # index 1 has NaN in the required column -> excluded from grouping, stays False
    assert list(out) == [True, False, False]


def test_opposite_signs_detects_sign_conflict_within_group():
    df = pd.DataFrame({"k": ["a", "a", "b"]})
    pts = pd.Series([1.0, -1.0, 1.0])
    out = opposite_signs(df, ["k"], pts)
    assert list(out) == [True, True, False]


def test_points_from_intervals_parses_matching_range():
    row = pd.Series({"_score": 0.5, "range_2": "0 1, -inf 0"})
    assert points_from_intervals(row, ["range_2"]) == 2


def test_points_from_intervals_no_match_returns_nan():
    row = pd.Series({"_score": 5.0, "range_2": "0 1"})
    assert pd.isna(points_from_intervals(row, ["range_2"]))


def test_parse_bound():
    assert parse_bound("-inf") == float("-inf")
    assert parse_bound("Infinity") == float("inf")
    assert parse_bound("1.5") == 1.5


def test_find_matching_class_matches_bracket_interval():
    row = pd.Series({"Interval 1 Range": "[0,1]", "Interval 1 Class": "Normal"})
    assert find_matching_class(row, 0.5) == "Normal"
    assert find_matching_class(row, 1.5) == "Not specified"


# ---------------- REVEL classification ----------------
def test_add_revel_classification_genome_wide_fallback(tmp_path):
    supp4 = tmp_path / "supp4.xlsx"
    pd.DataFrame({"Gene": []}).to_excel(supp4, sheet_name="REVEL_gene_specific_calibration", index=False)
    d = pd.DataFrame({"Gene": ["GENE1", "GENE1"], "REVEL": [0.01, 0.95], "points": [0, 0]})

    out = add_revel_classification(d, supp4)

    # REVEL 0.01 -> genome-wide BP4_Strong (-4 pts); no functional points -> total -4 -> Likely Benign
    assert out.loc[0, "Class_REVEL"] == "Likely Benign"
    # REVEL 0.95 -> genome-wide PP3_Strong (+4 pts); total 4 -> still Uncertain (needs >=10 for Pathogenic)
    assert out.loc[1, "Class_REVEL"] == "Uncertain"


def test_add_revel_classification_gene_specific_threshold_overrides_genome_wide(tmp_path):
    supp4 = tmp_path / "supp4.xlsx"
    pd.DataFrame({"Gene": ["GENE1"], "PP3_Strong": [0.5]}).to_excel(
        supp4, sheet_name="REVEL_gene_specific_calibration", index=False
    )
    d = pd.DataFrame({"Gene": ["GENE1"], "REVEL": [0.6], "points": [6]})

    out = add_revel_classification(d, supp4)

    # gene-specific PP3_Strong threshold (REVEL >= 0.5) adds 4 points; 6 (functional) + 4 = 10 -> Pathogenic
    assert out.loc[0, "Class_REVEL"] == "Pathogenic"
    assert out.loc[0, "Total_Points_REVEL"] == 10


def test_add_revel_classification_flags_functional_vs_revel_conflict(tmp_path):
    supp4 = tmp_path / "supp4.xlsx"
    pd.DataFrame({"Gene": []}).to_excel(supp4, sheet_name="REVEL_gene_specific_calibration", index=False)
    d = pd.DataFrame({"Gene": ["GENE1"], "REVEL": [0.01], "points": [5]})  # functional +, REVEL-derived -

    out = add_revel_classification(d, supp4)

    assert out.loc[0, "Conflicting_REVEL"] == "Conflicting evidence"


# ---------------- full pipeline ----------------
def test_filter_pillar_data_dedups_and_selects_category(tmp_path):
    rows = _points_df(
        [
            {
                "Dataset": "dsA",
                "Gene": "GENE1",
                "excalibr_points": "2",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "100",
                "ref_allele": "A",
                "alt_allele": "T",
                "clinvar_sig_2025": "Uncertain significance",
            },
            {  # same variant, smaller |points| -- should lose the dedup tie-break
                "Dataset": "dsB",
                "Gene": "GENE1",
                "excalibr_points": "1",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "100",
                "ref_allele": "A",
                "alt_allele": "T",
                "clinvar_sig_2025": "Uncertain significance",
            },
            {  # different variant, not VUS -- excluded from the VUS category
                "Dataset": "dsC",
                "Gene": "GENE1",
                "excalibr_points": "3",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "200",
                "ref_allele": "C",
                "alt_allele": "G",
                "clinvar_sig_2025": "Benign",
            },
        ]
    )
    points_input = tmp_path / "points.csv"
    rows.to_csv(points_input, index=False)

    dedup, revel = filter_pillar_data_fn(points_input, "VUS")

    assert len(dedup) == 1
    assert dedup.iloc[0]["Dataset"] == "dsA"
    assert "Class_REVEL" not in dedup.columns  # no --supp4-input -> REVEL classification skipped
    assert len(revel) == 1  # REVEL_train is "False" for both rows


def test_filter_pillar_data_removes_sfpq_and_flagged_rows(tmp_path):
    rows = _points_df(
        [
            {
                "Dataset": "dsA",
                "Gene": "SFPQ",
                "excalibr_points": "2",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "100",
                "ref_allele": "A",
                "alt_allele": "T",
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "dsB",
                "Gene": "GENE1",
                "excalibr_points": "2",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "200",
                "ref_allele": "C",
                "alt_allele": "G",
                "clinvar_sig_2025": "Uncertain significance",
                "Flag": "*",
            },
        ]
    )
    points_input = tmp_path / "points.csv"
    rows.to_csv(points_input, index=False)

    dedup, _ = filter_pillar_data_fn(points_input, "VUS")

    assert len(dedup) == 0


def test_filter_pillar_data_rejects_missing_columns(tmp_path):
    points_input = tmp_path / "points.csv"
    pd.DataFrame({"Dataset": ["ds1"]}).to_csv(points_input, index=False)

    with pytest.raises(ValueError, match="missing required column"):
        filter_pillar_data_fn(points_input, "VUS")


def test_cli_writes_category_and_revel_csvs(tmp_path):
    rows = _points_df(
        [
            {
                "Dataset": "dsA",
                "Gene": "GENE1",
                "excalibr_points": "2",
                "nucleotide_or_aa": "nt",
                "Chrom": "1",
                "hg38_start": "100",
                "ref_allele": "A",
                "alt_allele": "T",
                "clinvar_sig_2025": "Uncertain significance",
            }
        ]
    )
    points_input = tmp_path / "points.csv"
    rows.to_csv(points_input, index=False)
    outdir = tmp_path / "out"

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["--points-input", str(points_input), "--category", "VUS", "--outdir", str(outdir), "--validate-n", "0"],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "replicated_VUS_category.csv").exists()
    assert (outdir / "replicated_VUS_category_REVEL.csv").exists()
