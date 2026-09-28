import pandas as pd
import pytest
from click.testing import CliRunner

from src.add_excalibr_columns import (
    EXCALIBR_COLUMNS,
    _canonicalize_identity_hgvs_c,
    add_excalibr_columns,
    main,
)


def _integrated_df():
    return pd.DataFrame(
        {
            "Dataset": ["ds1", "ds1", "ds2"],
            "mavedb_variant_urn": ["urn1", "urn1", "urn2"],
            "hgvs_c": ["c.1A>T", "c.2A>T", "c.3A>T"],
            "Gene": ["GENE_A", "GENE_A", "GENE_B"],
        }
    )


def _excalibr_df(**overrides):
    data = {
        "Dataset": ["ds1", "ds1", "ds2"],
        "mavedb_variant_urn": ["urn1", "urn1", "urn2"],
        "hgvs_c": ["c.1A>T", "c.2A>T", "c.3A>T"],
        "Gene": ["GENE_A", "GENE_A", "GENE_B"],
        "excalibr_clinvar_release": ["2026", "2026", "2026"],
        "excalibr_points": ["1", "-1", "0"],
        "excalibr_acmg_evidence_code": ["PS3", "BS3", ""],
        "excalibr_lr_plus": ["10.0", "0.1", "1.0"],
        "excalibr_prior": ["0.1", "0.1", "0.1"],
        "excalibr_posterior": ["0.5", "0.01", "0.1"],
        "excalibr_filter_reason": ["", "", "low_coverage"],
    }
    data.update(overrides)
    return pd.DataFrame(data)


def test_add_excalibr_columns_joins_on_key():
    merged = add_excalibr_columns(_integrated_df(), _excalibr_df())

    assert list(merged["excalibr_points"]) == ["1", "-1", "0"]
    assert list(merged.columns) == ["Dataset", "mavedb_variant_urn", "hgvs_c", "Gene", *EXCALIBR_COLUMNS]


def test_add_excalibr_columns_missing_key_column_in_integrated_raises():
    integrated = _integrated_df().drop(columns=["hgvs_c"])
    with pytest.raises(ValueError, match="integrated dataset is missing required column"):
        add_excalibr_columns(integrated, _excalibr_df())


def test_add_excalibr_columns_missing_excalibr_column_raises():
    excalibr = _excalibr_df().drop(columns=["excalibr_points"])
    with pytest.raises(ValueError, match="ExCALIBR dataset is missing required column"):
        add_excalibr_columns(_integrated_df(), excalibr)


def test_add_excalibr_columns_duplicate_key_in_integrated_raises():
    integrated = pd.concat([_integrated_df(), _integrated_df().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="integrated dataset has duplicate rows"):
        add_excalibr_columns(integrated, _excalibr_df())


def test_add_excalibr_columns_duplicate_key_in_excalibr_raises():
    excalibr = pd.concat([_excalibr_df(), _excalibr_df().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="ExCALIBR dataset has duplicate rows"):
        add_excalibr_columns(_integrated_df(), excalibr)


def test_add_excalibr_columns_key_mismatch_raises():
    excalibr = _excalibr_df()
    excalibr.loc[2, "hgvs_c"] = "c.4A>T"
    with pytest.raises(ValueError, match="don't match between inputs"):
        add_excalibr_columns(_integrated_df(), excalibr)


# ---------------------------------------------------------------------------
# _canonicalize_identity_hgvs_c / identity-allele notation special case
# ---------------------------------------------------------------------------


def test_canonicalize_identity_hgvs_c_single_base():
    assert _canonicalize_identity_hgvs_c("NM_022445.4:c.612delCinsC") == "NM_022445.4:c.612C="


def test_canonicalize_identity_hgvs_c_multi_base():
    assert (
        _canonicalize_identity_hgvs_c("NM_022445.4:c.593_595delCCTinsCCT")
        == "NM_022445.4:c.593CCT="
    )


def test_canonicalize_identity_hgvs_c_already_equals_form_unchanged():
    assert _canonicalize_identity_hgvs_c("NM_022445.4:c.612C=") == "NM_022445.4:c.612C="


def test_canonicalize_identity_hgvs_c_real_variant_unchanged():
    # ref != alt -- a genuine delins, not a reference-identical allele.
    assert (
        _canonicalize_identity_hgvs_c("NM_022445.4:c.510_511delinsGG")
        == "NM_022445.4:c.510_511delinsGG"
    )
    assert (
        _canonicalize_identity_hgvs_c("NM_022445.4:c.510_511delCAinsGG")
        == "NM_022445.4:c.510_511delCAinsGG"
    )


def test_canonicalize_identity_hgvs_c_ordinary_substitution_unchanged():
    assert _canonicalize_identity_hgvs_c("NM_022445.4:c.10A>G") == "NM_022445.4:c.10A>G"


def test_add_excalibr_columns_matches_integrated_identity_allele_notation():
    """The integrated dataset's explicit del<bases>ins<bases> identity-allele
    hgvs_c (variant-annotation's map_variants.py fix) must still match
    ExCALIBR's un-updated single-anchor "=" form for the same allele -- and
    the integrated dataset's own hgvs_c value must survive in the output,
    unmodified.
    """
    integrated = _integrated_df()
    integrated.loc[2, "hgvs_c"] = "NM_000000.1:c.612delCinsC"
    excalibr = _excalibr_df()
    excalibr.loc[2, "hgvs_c"] = "NM_000000.1:c.612C="

    merged = add_excalibr_columns(integrated, excalibr)

    assert list(merged["hgvs_c"]) == ["c.1A>T", "c.2A>T", "NM_000000.1:c.612delCinsC"]
    assert list(merged["excalibr_points"]) == ["1", "-1", "0"]
    assert list(merged.columns) == ["Dataset", "mavedb_variant_urn", "hgvs_c", "Gene", *EXCALIBR_COLUMNS]


def test_add_excalibr_columns_real_variant_mismatch_still_raises():
    """A genuine hgvs_c difference (not the identity-allele special case)
    must still be caught as a real key mismatch."""
    integrated = _integrated_df()
    integrated.loc[2, "hgvs_c"] = "c.510_511delinsGG"
    excalibr = _excalibr_df()
    excalibr.loc[2, "hgvs_c"] = "c.510_511delinsCC"

    with pytest.raises(ValueError, match="don't match between inputs"):
        add_excalibr_columns(integrated, excalibr)


def test_main_cli_writes_output_and_supplementary_data_copy(tmp_path):
    integrated_path = tmp_path / "integrated.tsv.gz"
    excalibr_path = tmp_path / "excalibr.csv.gz"
    output_path = tmp_path / "output.tsv.gz"
    supplementary_data_path = tmp_path / "Supplementary_Data_1.tsv.gz"

    _integrated_df().to_csv(integrated_path, sep="\t", index=False)
    _excalibr_df().to_csv(excalibr_path, index=False)

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--integrated-input",
            str(integrated_path),
            "--excalibr-input",
            str(excalibr_path),
            "--output",
            str(output_path),
            "--supplementary-data-copy",
            str(supplementary_data_path),
        ],
    )

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["excalibr_points"]) == ["1", "-1", "0"]
    copy_out = pd.read_csv(supplementary_data_path, sep="\t", dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(out, copy_out)


def test_main_cli_skip_supplementary_data_copy(tmp_path):
    integrated_path = tmp_path / "integrated.tsv.gz"
    excalibr_path = tmp_path / "excalibr.csv.gz"
    output_path = tmp_path / "output.tsv.gz"
    supplementary_data_path = tmp_path / "Supplementary_Data_1.tsv.gz"

    _integrated_df().to_csv(integrated_path, sep="\t", index=False)
    _excalibr_df().to_csv(excalibr_path, index=False)

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--integrated-input",
            str(integrated_path),
            "--excalibr-input",
            str(excalibr_path),
            "--output",
            str(output_path),
            "--supplementary-data-copy",
            str(supplementary_data_path),
            "--skip-supplementary-data-copy",
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.exists()
    assert not supplementary_data_path.exists()


def test_main_cli_raises_click_exception_on_key_mismatch(tmp_path):
    integrated_path = tmp_path / "integrated.tsv.gz"
    excalibr_path = tmp_path / "excalibr.csv.gz"
    output_path = tmp_path / "output.tsv.gz"

    _integrated_df().to_csv(integrated_path, sep="\t", index=False)
    excalibr = _excalibr_df()
    excalibr.loc[2, "hgvs_c"] = "c.4A>T"
    excalibr.to_csv(excalibr_path, index=False)

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--integrated-input",
            str(integrated_path),
            "--excalibr-input",
            str(excalibr_path),
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code != 0
    assert "don't match between inputs" in result.output
