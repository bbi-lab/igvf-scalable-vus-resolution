import base64
import hashlib

import pandas as pd

from src.submit_supplementary_data_1_to_igvf import (
    build_attachment,
    build_document_payload,
    build_tabular_file_payload,
    md5sum,
    xlsx_schema_to_tsv,
)


def test_md5sum_matches_hashlib(tmp_path):
    path = tmp_path / "data.bin"
    path.write_bytes(b"some file content" * 1000)
    assert md5sum(path) == hashlib.md5(path.read_bytes()).hexdigest()


def test_xlsx_schema_to_tsv_strips_header_whitespace(tmp_path):
    xlsx_path = tmp_path / "schema.xlsx"
    tsv_path = tmp_path / "schema.tsv"
    pd.DataFrame(
        {
            "Column name ": ["Gene", "Chrom"],
            "Definition": ["HGNC gene symbol", "Chromosome"],
            "DataType": ["string", "string"],
        }
    ).to_excel(xlsx_path, index=False)

    xlsx_schema_to_tsv(xlsx_path, tsv_path)

    result = pd.read_csv(tsv_path, sep="\t")
    assert list(result.columns) == ["Column name", "Definition", "DataType"]
    assert result["Column name"].tolist() == ["Gene", "Chrom"]


def test_build_attachment_base64_encodes_content(tmp_path):
    tsv_path = tmp_path / "format_spec.tsv"
    tsv_path.write_text("Column name\tDefinition\tDataType\nGene\tHGNC gene symbol\tstring\n")

    attachment = build_attachment(tsv_path, mime_type="text/tab-separated-values")

    assert attachment["download"] == "format_spec.tsv"
    assert attachment["type"] == "text/tab-separated-values"
    prefix = "data:text/tab-separated-values;base64,"
    assert attachment["href"].startswith(prefix)
    decoded = base64.b64decode(attachment["href"][len(prefix) :])
    assert decoded == tsv_path.read_bytes()


def test_build_document_payload_includes_supersedes_and_attachment():
    attachment = {"download": "f.tsv", "type": "text/tab-separated-values", "href": "data:..."}

    payload = build_document_payload(attachment, alias="lea-starita:test-document")

    assert payload["aliases"] == ["lea-starita:test-document"]
    assert payload["attachment"] == attachment
    assert payload["supersedes"] == ["/documents/2bfe67c5-937a-4fe0-b343-5008bdbf7241/"]
    assert payload["document_type"] == "file format specification"
    assert "standardized_file_format" not in payload  # admin-only; submitting it is a 422


def test_build_tabular_file_payload_includes_supersedes_and_md5():
    payload = build_tabular_file_payload(
        md5="d41d8cd98f00b204e9800998ecf8427e",
        file_size=12345,
        aliases=["lea-starita:test-file"],
        submitted_file_name="/path/to/file.tsv.gz",
        file_format_specification_id="/documents/new-doc-uuid/",
    )

    assert payload["md5sum"] == "d41d8cd98f00b204e9800998ecf8427e"
    assert payload["file_size"] == 12345
    assert payload["file_format_specifications"] == ["/documents/new-doc-uuid/"]
    assert payload["supersedes"] == ["/tabular-files/IGVFFI2521UGYG/"]
    assert payload["controlled_access"] is False
    assert payload["reference_files"] == ["/reference-files/IGVFFI0653VCGH/"]


def test_build_tabular_file_payload_omits_reference_files_when_none():
    payload = build_tabular_file_payload(
        md5="d41d8cd98f00b204e9800998ecf8427e",
        file_size=12345,
        aliases=["lea-starita:test-file"],
        submitted_file_name="file.tsv.gz",
        file_format_specification_id="/documents/new-doc-uuid/",
        reference_files=None,
    )

    assert "reference_files" not in payload
