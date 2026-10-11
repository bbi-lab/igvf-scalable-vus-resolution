import math

import pandas as pd
import pytest

from src.make_extended_data_figure_16 import GROUPS, analyze, bootstrap_ci, load, make_figures, ztest


def _write_csv(path, rows):
    """rows: list of {"Class_REVEL": ..., "AF_<group>": ...} dicts. Any
    AF_<group> not given for a row defaults to 0."""
    df = pd.DataFrame(rows)
    for g in GROUPS:
        col = f"AF_{g}"
        if col not in df.columns:
            df[col] = 0
    df.to_csv(path, index=False)
    return df


def test_load_rejects_missing_columns(tmp_path):
    path = tmp_path / "data.csv"
    pd.DataFrame({"Class_REVEL": ["Benign"]}).to_csv(path, index=False)

    with pytest.raises(ValueError, match="missing columns"):
        load(path)


def test_load_coerces_af_columns_to_numeric(tmp_path):
    path = tmp_path / "data.csv"
    rows = [{"Class_REVEL": "Benign", "AF_nfe": "not_a_number", "AF_afr": "0.01"}]
    _write_csv(path, rows)

    d = load(path)

    assert d["AF_nfe"].iloc[0] == 0
    assert d["AF_afr"].iloc[0] == pytest.approx(0.01)


def test_analyze_drops_variants_with_no_ancestry_observation(tmp_path):
    path = tmp_path / "data.csv"
    rows = [
        {"Class_REVEL": "Benign", "AF_nfe": 0.01},  # observed in nfe
        {"Class_REVEL": "Uncertain", "AF_nfe": 0, "AF_afr": 0},  # unobserved everywhere -- dropped
    ]
    _write_csv(path, rows)
    d = load(path)

    results = []
    analyze(d, "gnomAD", n_bootstrap=100, results=results)

    all_present = [r for r in results if r["framing"] == "present" and r["group"] == "ALL"]
    assert len(all_present) == 1
    assert all_present[0]["n"] == 1


def test_analyze_present_vs_unique_framing(tmp_path):
    path = tmp_path / "data.csv"
    rows = [
        # shared by afr and nfe -- counts in "present" for both, "unique" for neither
        {"Class_REVEL": "Benign", "AF_afr": 0.01, "AF_nfe": 0.01},
        # unique to nfe
        {"Class_REVEL": "Pathogenic", "AF_nfe": 0.01},
        {"Class_REVEL": "Uncertain", "AF_nfe": 0.01},
    ]
    _write_csv(path, rows)
    d = load(path)

    results = []
    analyze(d, "gnomAD", n_bootstrap=100, results=results)
    results_df = pd.DataFrame(results)

    present_nfe = results_df[(results_df.framing == "present") & (results_df.group == "nfe")].iloc[0]
    unique_nfe = results_df[(results_df.framing == "unique") & (results_df.group == "nfe")].iloc[0]
    # afr has zero unique-only variants -- analyze() omits n=0 group rows entirely.
    unique_afr = results_df[(results_df.framing == "unique") & (results_df.group == "afr")]

    assert present_nfe["n"] == 3
    assert unique_nfe["n"] == 2
    assert unique_afr.empty


def test_ztest_matches_known_two_proportion_result():
    # Equal rates -> p should be close to 1 (not exactly 1 due to pooling).
    p_equal = ztest(50, 100, 50, 100)
    assert p_equal == pytest.approx(1.0)

    # Clearly different rates -> small p-value.
    p_diff = ztest(90, 100, 10, 100)
    assert p_diff < 0.001


def test_ztest_handles_zero_denominator():
    assert math.isnan(ztest(0, 0, 5, 10))


def test_bootstrap_ci_brackets_the_point_estimate():
    values = [1, 1, 1, 0, 0]  # rate = 0.6
    lo, hi = bootstrap_ci(values, n_boot=500)
    assert lo <= 0.6 <= hi


def test_bootstrap_ci_empty_input_returns_nan():
    lo, hi = bootstrap_ci([], n_boot=100)
    assert math.isnan(lo)
    assert math.isnan(hi)


def test_make_figures_writes_a_png_per_variant_set(tmp_path):
    vus_path = tmp_path / "vus.csv"
    gnomad_path = tmp_path / "gnomad.csv"
    _write_csv(
        vus_path,
        [
            {"Class_REVEL": "Benign", "AF_nfe": 0.01},
            {"Class_REVEL": "Uncertain", "AF_afr": 0.01},
        ],
    )
    _write_csv(
        gnomad_path,
        [
            {"Class_REVEL": "Pathogenic", "AF_nfe": 0.01},
            {"Class_REVEL": "Uncertain", "AF_sas": 0.01},
        ],
    )

    results = []
    analyze(load(vus_path), "VUS", n_bootstrap=50, results=results)
    analyze(load(gnomad_path), "gnomAD", n_bootstrap=50, results=results)
    results_df = pd.DataFrame(results)

    outdir = tmp_path / "out"
    make_figures(results_df, outdir, width=10, height=6, no_stars=False, ci_label="")

    assert (outdir / "ancestry_resolution_VUS.png").exists()
    assert (outdir / "ancestry_resolution_gnomAD.png").exists()
