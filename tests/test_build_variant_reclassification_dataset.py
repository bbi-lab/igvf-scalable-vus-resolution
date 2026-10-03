from pathlib import Path

import pandas as pd
import pytest

from src.build_variant_reclassification_dataset import (
    BASE_OUTPUT_COLUMNS,
    OUTPUT_COLUMNS,
    add_points_columns,
    apply_notebook_exclusions,
    build_reclassification_dataset,
    points_output_columns,
)

CHECKPOINT_COLS = BASE_OUTPUT_COLUMNS + [
    "VariantNotes", "ExC_points_2025", "ExC_points_2018", "OP_points", "Fxn_points",
    "revel_train_amino", "mp2_train_amino",
    "Points_REVEL_GeneSpecific_GenomeWide", "Points_REVEL_GenomeWide",
    "Points_AM_GeneSpecific_GenomeWide", "Points_AM_GenomeWide",
    "Points_MP2_GeneSpecific_GenomeWide", "Points_MP2_GenomeWide",
    "Conflicting_REVEL_GeneSpecific", "Conflicting_REVEL_GenomeWide",
    "Conflicting_AM_GeneSpecific", "Conflicting_AM_GenomeWide",
    "Conflicting_MP2_GeneSpecific", "Conflicting_MP2_GenomeWide",
]


def _checkpoint_row(**overrides):
    row = {col: None for col in CHECKPOINT_COLS}
    row.update({
        "mavedb_variant_urn": "urn:mavedb:1",
        "Dataset": "G1_Dataset",
        "Gene": "G1",
        "Chrom": "1",
        "hg38_start": 1000,
        "ref_allele": "A",
        "alt_allele": "T",
        "hgvs_p": "p.NoMatch1Xxx",
        "auth_reported_score": 99.0,
        "Flag": None,
        "VariantNotes": None,
        "splice_var_amino": "No",
        "ExC_points_2025": 5,
        "ExC_points_2018": None,
        "OP_points": None,
        "Fxn_points": 5,
        "Points_REVEL_GeneSpecific_GenomeWide": 2,
        "Conflicting_REVEL_GeneSpecific": 7,
        "revel_train_amino": "No",
        "mp2_train_amino": "No",
    })
    row.update(overrides)
    return row


def _checkpoint_frame(rows):
    return pd.DataFrame([_checkpoint_row(**r) for r in rows])


@pytest.fixture
def chek2_file(tmp_path):
    path = tmp_path / "CHEK2_Gebbia_2024.xlsx"
    pd.DataFrame(
        {
            "hgvs_pro": ["p.Val1Ala", "p.Val2Ala"],
            "score": [1.0, 1.0],
            "Filter_CI": [1, 0],
        }
    ).to_excel(path, index=False)
    return path


# --- apply_notebook_exclusions -----------------------------------------------------------


def test_sfpq_dropped(chek2_file):
    df = _checkpoint_frame([{"Gene": "SFPQ"}, {"Gene": "G1"}])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["Gene"]) == ["G1"]


def test_chek2_flagged_row_dropped(chek2_file):
    """hgvs_p (transcript-prefixed, as in the real checkpoint)/
    auth_reported_score matching a Filter_CI==1 CHEK2 row (unprefixed
    hgvs_pro) gets Flag='*' and is then removed by the Flag!='*' filter."""
    df = _checkpoint_frame([
        {"Gene": "CHEK2", "hgvs_p": "NP_009125.1:p.Val1Ala", "auth_reported_score": 1.0},
        {"Gene": "CHEK2", "hgvs_p": "NP_009125.1:p.Val2Ala", "auth_reported_score": 1.0},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert len(out) == 1
    assert out["hgvs_p"].iloc[0] == "NP_009125.1:p.Val2Ala"


def test_conflicting_fxn_data_kept_by_default(chek2_file):
    """Bare `conflicting_fxn_data` is not excluded by default: each row of
    this export is one dataset's measurement, so more than one dataset
    disagreeing on a variant's effect is expected and kept, not treated as
    disqualifying."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "VariantNotes": "conflicting_fxn_data"},
        {"mavedb_variant_urn": "urn:mavedb:2", "VariantNotes": "max_fxn_pts"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["VariantNotes"]) == ["conflicting_fxn_data", "max_fxn_pts"]


def test_conflicting_fxn_data_dropped_when_dedup_enabled(chek2_file):
    """With `dedup=True`, deduplication would otherwise silently pick one
    dataset's value as the winner for a variant with genuinely conflicting
    measurements, so bare `conflicting_fxn_data` rows are excluded instead."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "VariantNotes": "conflicting_fxn_data"},
        {"mavedb_variant_urn": "urn:mavedb:2", "VariantNotes": "max_fxn_pts"},
    ])
    out = apply_notebook_exclusions(df, chek2_file, dedup=True)
    assert list(out["VariantNotes"]) == ["max_fxn_pts"]


def test_splice_variant_not_measured_dropped(chek2_file):
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "VariantNotes": "splice_variant_not_measured"},
        {"mavedb_variant_urn": "urn:mavedb:2", "VariantNotes": "max_fxn_pts"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["VariantNotes"]) == ["max_fxn_pts"]


def test_splice_variant_not_measured_conflicting_fxn_data_dropped(chek2_file):
    """The compound tag is still dropped even though bare `conflicting_fxn_data`
    is not -- it also carries `splice_variant_not_measured`, an assay
    limitation rather than a cross-dataset conflict."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "VariantNotes": "splice_variant_not_measured;conflicting_fxn_data"},
        {"mavedb_variant_urn": "urn:mavedb:2", "VariantNotes": "max_fxn_pts"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["VariantNotes"]) == ["max_fxn_pts"]


def test_splice_var_amino_dropped(chek2_file):
    df = _checkpoint_frame([
        {"splice_var_amino": "Yes"},
        {"splice_var_amino": "No"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["splice_var_amino"]) == ["No"]


def test_splice_var_amino_kept_when_splice_measure_yes(chek2_file):
    """A dataset curated as able to detect splicing effects
    (`splice_measure == 'Yes'`) keeps its splice-flagged rows instead of
    having them dropped unconditionally."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "splice_var_amino": "Yes", "splice_measure": "Yes"},
        {"mavedb_variant_urn": "urn:mavedb:2", "splice_var_amino": "Yes", "splice_measure": "No"},
        {"mavedb_variant_urn": "urn:mavedb:3", "splice_var_amino": "Yes", "splice_measure": None},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["mavedb_variant_urn"]) == ["urn:mavedb:1"]


def test_pre_existing_flag_still_removed(chek2_file):
    df = _checkpoint_frame([
        {"Flag": "*"},
        {"Flag": None},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert len(out) == 1
    assert pd.isna(out["Flag"].iloc[0])


def test_revel_train_amino_dropped(chek2_file):
    """Variants used to train REVEL are excluded by default (predictor=
    REVEL), matching the notebook's REVEL-specific category sheets
    (VUS_REVEL, Unobserved_REVEL, etc.), which drop them to avoid
    circularity with this file's REVEL-based Combined_points."""
    df = _checkpoint_frame([
        {"revel_train_amino": "Yes"},
        {"revel_train_amino": "No"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["revel_train_amino"]) == ["No"]


def test_mp2_train_amino_dropped_when_predictor_is_mp2(chek2_file):
    """With predictor=MP2, the circularity concern shifts to MP2's own
    training flag -- REVEL's is no longer relevant, since Combined_points
    is no longer REVEL-based."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "mp2_train_amino": "Yes", "revel_train_amino": "No"},
        {"mavedb_variant_urn": "urn:mavedb:2", "mp2_train_amino": "No", "revel_train_amino": "Yes"},
    ])
    out = apply_notebook_exclusions(df, chek2_file, predictor="MP2")
    assert list(out["mavedb_variant_urn"]) == ["urn:mavedb:2"]


def test_no_train_exclusion_when_predictor_is_am(chek2_file):
    """AM has no training-circularity column in the checkpoint -- neither
    REVEL's nor MP2's training flag should exclude anything when
    predictor=AM."""
    df = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "revel_train_amino": "Yes", "mp2_train_amino": "Yes"},
        {"mavedb_variant_urn": "urn:mavedb:2", "revel_train_amino": "No", "mp2_train_amino": "No"},
    ])
    out = apply_notebook_exclusions(df, chek2_file, predictor="AM")
    assert set(out["mavedb_variant_urn"]) == {"urn:mavedb:1", "urn:mavedb:2"}


# --- add_points_columns -------------------------------------------------------------------


def test_excalibr_points_is_literal_value_even_when_op_overrides_functional():
    """TP53 uses OddsPath for Functional_points, but ExCALIBR_points should
    still hold the literal ExC_points_2025 value."""
    df = _checkpoint_frame([
        {"Gene": "TP53", "ExC_points_2025": 4, "OP_points": 8, "Fxn_points": 8},
    ])
    out = add_points_columns(df)
    assert out["ExCALIBR_points"].iloc[0] == 4
    assert out["OddsPath_points"].iloc[0] == 8
    assert out["Functional_points"].iloc[0] == 8


def test_excalibr_points_vintage_override_for_brca1_pten_msh2_tp53():
    df = _checkpoint_frame([
        {"Gene": "BRCA1", "ExC_points_2025": 2, "ExC_points_2018": 6},
        {"Gene": "TP53", "ExC_points_2025": 2, "ExC_points_2018": 6},
        {"Gene": "G1", "ExC_points_2025": 2, "ExC_points_2018": 6},
    ])
    out = add_points_columns(df)
    assert out.loc[out.Gene == "BRCA1", "ExCALIBR_points"].iloc[0] == 6
    assert out.loc[out.Gene == "TP53", "ExCALIBR_points"].iloc[0] == 6
    assert out.loc[out.Gene == "G1", "ExCALIBR_points"].iloc[0] == 2


def test_excalibr_points_vintage_override_falls_back_when_2018_missing():
    df = _checkpoint_frame([{"Gene": "PTEN", "ExC_points_2025": 3, "ExC_points_2018": None}])
    out = add_points_columns(df)
    assert out["ExCALIBR_points"].iloc[0] == 3


def test_combined_points_is_functional_plus_revel_gene_specific():
    df = _checkpoint_frame([{"Fxn_points": 5, "Points_REVEL_GeneSpecific_GenomeWide": 2}])
    out = add_points_columns(df)
    assert out["Combined_points"].iloc[0] == 7


def test_combined_points_treats_missing_revel_as_zero():
    df = _checkpoint_frame([{"Fxn_points": 5, "Points_REVEL_GeneSpecific_GenomeWide": None}])
    out = add_points_columns(df)
    assert out["Combined_points"].iloc[0] == 5


def test_revel_points_is_renamed_from_points_revel_gene_specific_genome_wide():
    df = _checkpoint_frame([{"Points_REVEL_GeneSpecific_GenomeWide": 3}])
    out = add_points_columns(df)
    assert out["REVEL_points"].iloc[0] == 3


def test_conflict_revel_gene_specific_is_renamed_from_conflicting_revel_gene_specific():
    df = _checkpoint_frame([{"Conflicting_REVEL_GeneSpecific": "Conflicting evidence"}])
    out = add_points_columns(df)
    assert out["Conflict_REVEL_GeneSpecific"].iloc[0] == "Conflicting evidence"


def test_points_output_columns_reduces_to_historical_names_for_default():
    assert points_output_columns("REVEL", "gene-specific") == [
        "ExCALIBR_points", "OddsPath_points", "Functional_points",
        "REVEL_points", "Conflict_REVEL_GeneSpecific", "Combined_points",
    ]


def test_add_points_columns_predictor_am_uses_am_gene_specific_columns():
    df = _checkpoint_frame([{
        "Fxn_points": 5,
        "Points_AM_GeneSpecific_GenomeWide": 2,
        "Conflicting_AM_GeneSpecific": "Conflicting evidence",
    }])
    out = add_points_columns(df, predictor="AM")
    assert out["AM_points"].iloc[0] == 2
    assert out["Conflict_AM_GeneSpecific"].iloc[0] == "Conflicting evidence"
    assert out["Combined_points"].iloc[0] == 7


def test_add_points_columns_predictor_mp2_genome_wide_scope():
    df = _checkpoint_frame([{
        "Fxn_points": 5,
        "Points_MP2_GenomeWide": 3,
        "Conflicting_MP2_GenomeWide": "Conflicting evidence",
    }])
    out = add_points_columns(df, predictor="MP2", predictor_scope="genome-wide")
    assert out["MP2_points"].iloc[0] == 3
    assert out["Conflict_MP2_GenomeWide"].iloc[0] == "Conflicting evidence"
    assert out["Combined_points"].iloc[0] == 8


# --- build_reclassification_dataset (end to end) -------------------------------------------


def test_build_reclassification_dataset_keeps_one_row_per_measurement(tmp_path, chek2_file):
    """Unlike the pipeline's per-category exports, this one is not
    deduplicated to one row per DNA variant: more than one dataset scoring
    the same variant should survive as separate rows."""
    checkpoint = _checkpoint_frame([
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_low", "Fxn_points": 2,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_high", "Fxn_points": 9,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
    ])
    checkpoint_path = tmp_path / "checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")

    out = build_reclassification_dataset(checkpoint_path, chek2_file)
    assert len(out) == 2
    assert set(out["Dataset"]) == {"Dataset_low", "Dataset_high"}
    assert list(out.columns) == OUTPUT_COLUMNS


def test_build_reclassification_dataset_dedups_by_dna_variant_when_enabled(tmp_path, chek2_file):
    checkpoint = _checkpoint_frame([
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_low", "Fxn_points": 2,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_high", "Fxn_points": 9,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
    ])
    checkpoint_path = tmp_path / "checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")

    out = build_reclassification_dataset(checkpoint_path, chek2_file, dedup=True)
    assert len(out) == 1
    assert out["Dataset"].iloc[0] == "Dataset_high"
    assert list(out.columns) == OUTPUT_COLUMNS


def test_build_reclassification_dataset_predictor_can_change_dedup_winner(tmp_path, chek2_file):
    """Combined_points (the dedup sort key) is predictor-dependent, so the
    surviving measurement for a variant scored by more than one dataset can
    differ by --predictor even when Functional_points is tied."""
    checkpoint = _checkpoint_frame([
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_low_revel_high_am", "Fxn_points": 2,
         "Points_REVEL_GeneSpecific_GenomeWide": 0, "Points_AM_GeneSpecific_GenomeWide": 8},
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_high_revel_low_am", "Fxn_points": 2,
         "Points_REVEL_GeneSpecific_GenomeWide": 5, "Points_AM_GeneSpecific_GenomeWide": 0},
    ])
    checkpoint_path = tmp_path / "checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")

    revel_out = build_reclassification_dataset(checkpoint_path, chek2_file, dedup=True, predictor="REVEL")
    am_out = build_reclassification_dataset(checkpoint_path, chek2_file, dedup=True, predictor="AM")

    assert revel_out["Dataset"].iloc[0] == "Dataset_high_revel_low_am"
    assert am_out["Dataset"].iloc[0] == "Dataset_low_revel_high_am"
    assert list(am_out.columns) == BASE_OUTPUT_COLUMNS + points_output_columns("AM", "gene-specific")


def test_build_reclassification_dataset_dedup_excludes_conflicting_fxn_data(tmp_path, chek2_file):
    checkpoint = _checkpoint_frame([
        {"mavedb_variant_urn": "urn:mavedb:1", "VariantNotes": "conflicting_fxn_data"},
        {"mavedb_variant_urn": "urn:mavedb:2", "VariantNotes": "max_fxn_pts"},
    ])
    checkpoint_path = tmp_path / "checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")

    out = build_reclassification_dataset(checkpoint_path, chek2_file, dedup=True)
    assert list(out["mavedb_variant_urn"]) == ["urn:mavedb:2"]
