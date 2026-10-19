import csv
import gzip

from src.build_training_variant_files import build_revel_training_file

_HEADER = ["gene_symbol", "unqualified_hgvs_p", "chromosome", "hg38_start", "hg38_end", "ref_allele", "alt_allele"]


def test_build_revel_training_file_writes_rows_from_source(tmp_path, capsys):
    src_path = tmp_path / "revel_training_overlap.csv.gz"
    dest_path = tmp_path / "revel_training_variants.tsv"

    with gzip.open(src_path, "wt", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["#CHROM", "POS", "REF", "ALT", "protein_variant", "short_id", "gene_symbol", "ensembl_prot_id"]
        )
        writer.writerow(["chr1", "32780173", "G", "T", "Q416K", "ENST00000373477", "YARS1", "ENSP00000362576.4"])

    build_revel_training_file(src_path, dest_path)

    with dest_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_symbol"] == "YARS1"
    assert rows[0]["unqualified_hgvs_p"] == "p.Gln416Lys"
    assert rows[0]["chromosome"] == "1"


def test_build_revel_training_file_skips_missing_source_without_raising(tmp_path, capsys):
    src_path = tmp_path / "does_not_exist.csv.gz"
    dest_path = tmp_path / "revel_training_variants.tsv"

    build_revel_training_file(src_path, dest_path)

    assert dest_path.exists()
    with dest_path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    assert rows == [_HEADER]
    assert "not found" in capsys.readouterr().out


def test_build_revel_training_file_missing_source_output_loads_as_empty_training_set(tmp_path):
    """The downstream loader (annotate_predictors._load_revel_training_variants) only requires
    the header row be present with the expected columns -- confirm the skip path produces that,
    not a truly empty (unparseable) file."""
    src_path = tmp_path / "does_not_exist.csv.gz"
    dest_path = tmp_path / "revel_training_variants.tsv"

    build_revel_training_file(src_path, dest_path)

    with dest_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert reader.fieldnames is not None
        required = ["chromosome", "hg38_start", "hg38_end", "ref_allele", "alt_allele"]
        assert all(col in reader.fieldnames for col in required)
        assert list(reader) == []
