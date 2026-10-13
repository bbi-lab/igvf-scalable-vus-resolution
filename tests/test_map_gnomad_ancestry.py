import shutil

import pandas as pd
import pysam
import pytest
from click.testing import CliRunner

from src.map_gnomad_ancestry import (
    clean_allele,
    left_align,
    main,
    map_gnomad_ancestry,
    match_record,
    normalize_row,
    parse_int,
    ref_check,
    revcomp,
    to_vcf,
    variant_type,
)

CHROM = "chr1"
REF_SEQ = "N" * 99 + "A" + "CGTACGTACG" + "N" * 890  # position 100 (1-based) == "A"


def _seq(fasta_path):
    fasta = pysam.FastaFile(fasta_path)

    def seq(chrom, start, end):
        return fasta.fetch(chrom, start - 1, end).upper()

    return seq


def _write_fasta(tmp_path):
    path = tmp_path / "ref.fa"
    with open(path, "w") as f:
        f.write(f">{CHROM}\n{REF_SEQ}\n")
    pysam.faidx(str(path))
    return path


def _write_exomes_vcf(tmp_path, records):
    """records: list of (pos, ref, alt, info_dict) tuples. Returns the exomes dir."""
    exomes_dir = tmp_path / "exomes"
    exomes_dir.mkdir()
    vcf_path = exomes_dir / f"gnomad.exomes.v4.1.sites.{CHROM}.vcf"
    info_ids = sorted({k for _, _, _, info in records for k in info})
    header_lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length={len(REF_SEQ)}>",
    ]
    for info_id in info_ids:
        number = "A" if info_id.startswith(("AC", "AF")) else "1"
        type_ = "Float" if info_id.startswith("AF") else "Integer"
        header_lines.append(f'##INFO=<ID={info_id},Number={number},Type={type_},Description="">')
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    with open(vcf_path, "w") as f:
        f.write("\n".join(header_lines) + "\n")
        for pos, ref, alt, info in records:
            info_str = ";".join(f"{k}={v}" for k, v in info.items())
            f.write(f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info_str}\n")
    bgz_path = pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    # gnomAD's own naming convention ends in .vcf.bgz, not .vcf.gz
    final_path = exomes_dir / f"gnomad.exomes.v4.1.sites.{CHROM}.vcf.bgz"
    shutil.move(bgz_path, final_path)
    shutil.move(bgz_path + ".tbi", str(final_path) + ".tbi")
    return exomes_dir


# ---------------- pure allele-normalization helpers ----------------
def test_clean_allele_normalizes_empty_markers():
    assert clean_allele("-") == ""
    assert clean_allele(None) == ""
    assert clean_allele(float("nan")) == ""
    assert clean_allele(" a ") == "A"


def test_revcomp():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AATT") == "AATT"
    assert revcomp("GATTACA") == "TGTAATC"


def test_parse_int():
    assert parse_int("100") == 100
    assert parse_int("100.0") == 100
    assert parse_int("not a number") is None
    assert parse_int(None) is None


def test_ref_check_match_revcomp_and_mismatch(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    assert ref_check(CHROM, 100, "A", seq) == "match"
    assert ref_check(CHROM, 100, "T", seq) == "revcomp_match"
    assert ref_check(CHROM, 100, "G", seq) == "mismatch"
    assert ref_check(CHROM, 100, "", seq) == "no_ref_allele"


def test_to_vcf_anchors_empty_alleles(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    # deletion: ref given, alt empty -> anchor to the preceding base
    assert to_vcf(CHROM, 101, None, "C", "", seq) == (100, "AC", "A")
    # insertion: alt given, ref empty, end given -> anchor to min(start, end)
    assert to_vcf(CHROM, 101, 101, "", "G", seq) == (101, "C", "CG")
    # both present: passed through unchanged
    assert to_vcf(CHROM, 100, 100, "A", "T", seq) == (100, "A", "T")
    # both empty: invalid
    assert to_vcf(CHROM, 100, 100, "", "", seq) is None


def test_left_align_trims_matching_right_base(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    # "AC" -> "A" at pos 100 is a C-deletion; already left-aligned
    assert left_align(CHROM, 100, "AC", "A", seq) == (100, "AC", "A")


def test_variant_type_classification():
    assert variant_type("A", "T") == "SNV"
    assert variant_type("AC", "TG") == "MNV"
    assert variant_type("AC", "A") == "deletion"
    assert variant_type("A", "AC") == "insertion"
    assert variant_type("AC", "GTA") == "complex_indel"


def test_normalize_row_rejects_unparseable_and_identical_alleles(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    assert normalize_row(CHROM, None, None, "A", "T", seq) == (None, "unparseable", "n/a")
    assert normalize_row(CHROM, 100, 100, "A", "A", seq) == (None, "invalid_alleles", "n/a")


def test_normalize_row_snv(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    key, vt, rc = normalize_row(CHROM, 100, 100, "A", "T", seq)
    assert key == (CHROM, 100, "A", "T")
    assert vt == "SNV"
    assert rc == "match"


# ---------------- gnomAD VCF matching ----------------
def test_match_record_finds_exact_allele(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    exomes_dir = _write_exomes_vcf(
        tmp_path, [(100, "A", "T", {"AC": 5, "AN": 1000, "AF": 0.005, "AC_nfe": 3, "AN_nfe": 600, "AF_nfe": 0.005})]
    )
    vcf = pysam.VariantFile(str(exomes_dir / f"gnomad.exomes.v4.1.sites.{CHROM}.vcf.bgz"))
    candidates = list(vcf.fetch(CHROM, 90, 110))
    rec, alt_i = match_record(candidates, (CHROM, 100, "A", "T"), seq)
    assert rec is not None
    assert alt_i == 0
    assert rec.info["AF_nfe"][0] == pytest.approx(0.005, abs=1e-6)


def test_match_record_no_match_for_different_allele(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    seq = _seq(str(fasta_path))
    exomes_dir = _write_exomes_vcf(tmp_path, [(100, "A", "T", {"AC": 5, "AN": 1000, "AF": 0.005})])
    vcf = pysam.VariantFile(str(exomes_dir / f"gnomad.exomes.v4.1.sites.{CHROM}.vcf.bgz"))
    candidates = list(vcf.fetch(CHROM, 90, 110))
    rec, alt_i = match_record(candidates, (CHROM, 100, "A", "G"), seq)
    assert rec is None
    assert alt_i is None


# ---------------- end-to-end ----------------
def _write_input_csv(path):
    pd.DataFrame(
        {
            "Chrom": ["1", "1"],
            "hg38_start": [100, 500],
            "ref_allele": ["A", "C"],
            "alt_allele": ["T", "G"],
        }
    ).to_csv(path, index=False)


def test_map_gnomad_ancestry_end_to_end(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    exomes_dir = _write_exomes_vcf(
        tmp_path,
        [(100, "A", "T", {"AC": 5, "AN": 1000, "AF": 0.005, "AC_nfe": 3, "AN_nfe": 600, "AF_nfe": 0.005})],
    )
    input_path = tmp_path / "variants.csv"
    _write_input_csv(input_path)

    df = map_gnomad_ancestry(input_path, fasta_path, exomes_dir, None, workers=1)
    df["hg38_start"] = df["hg38_start"].astype(int)

    assert df.loc[df["hg38_start"] == 100, "exomes_found"].iloc[0]
    assert df.loc[df["hg38_start"] == 100, "AF_nfe"].iloc[0] == pytest.approx(0.005, abs=1e-6)
    assert not df.loc[df["hg38_start"] == 500, "exomes_found"].iloc[0]
    assert not df.loc[df["hg38_start"] == 500, "gnomad_found_any"].iloc[0]


def test_map_gnomad_ancestry_rejects_missing_columns(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    exomes_dir = _write_exomes_vcf(tmp_path, [(100, "A", "T", {"AC": 1, "AN": 2, "AF": 0.5})])
    input_path = tmp_path / "variants.csv"
    pd.DataFrame({"Chrom": ["1"]}).to_csv(input_path, index=False)

    with pytest.raises(ValueError, match="missing required column"):
        map_gnomad_ancestry(input_path, fasta_path, exomes_dir, None, workers=1)


def test_cli_end_to_end(tmp_path):
    fasta_path = _write_fasta(tmp_path)
    exomes_dir = _write_exomes_vcf(
        tmp_path,
        [(100, "A", "T", {"AC": 5, "AN": 1000, "AF": 0.005, "AC_nfe": 3, "AN_nfe": 600, "AF_nfe": 0.005})],
    )
    input_path = tmp_path / "variants.csv"
    _write_input_csv(input_path)
    output_path = tmp_path / "out.csv"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ref-fasta",
            str(fasta_path),
            "--exomes-dir",
            str(exomes_dir),
            "--workers",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.exists()
    out = pd.read_csv(output_path)
    assert out.loc[out["hg38_start"] == 100, "exomes_found"].iloc[0]


def test_cli_requires_fasta_index(tmp_path):
    fasta_path = tmp_path / "ref.fa"
    with open(fasta_path, "w") as f:
        f.write(f">{CHROM}\n{REF_SEQ}\n")
    # deliberately not indexed
    exomes_dir = tmp_path / "exomes"
    exomes_dir.mkdir()
    input_path = tmp_path / "variants.csv"
    _write_input_csv(input_path)

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--input",
            str(input_path),
            "--output",
            str(tmp_path / "out.csv"),
            "--ref-fasta",
            str(fasta_path),
            "--exomes-dir",
            str(exomes_dir),
        ],
    )

    assert result.exit_code != 0
    assert "samtools faidx" in result.output
