import zipfile

import pytest

from src.lib.zenodo_client import render_manifest_summary
from src.submit_data_prerequisites_to_zenodo import (
    _EXCALIBR_JSON_ZIP_KEY,
    _GLOB_DIRS,
    _STATIC_FILES,
    FALLBACK_GITHUB_URL,
    MP2_PUBLIC_SOURCE_KEY,
    PREPRINT_DOI,
    build_file_manifest,
    build_related_identifiers,
)


def _make_complete_repo(root):
    """Create every file/dir `build_file_manifest` expects, each with distinct content."""
    for local_rel, _key, _restricted in _STATIC_FILES:
        path = root / local_rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(local_rel.encode())

    for dir_rel, _prefix in _GLOB_DIRS:
        dir_path = root / dir_rel
        dir_path.mkdir(parents=True, exist_ok=True)
        (dir_path / ".DS_Store").write_bytes(b"junk")
        (dir_path / "a.xlsx").write_bytes(b"a")
        (dir_path / "b.pdf").write_bytes(b"b")

    json_dir = root / "data/input/mave_calibration/excalibr/json"
    json_dir.mkdir(parents=True, exist_ok=True)
    (json_dir / "DatasetA.json").write_text('{"prior": 1}')
    (json_dir / "DatasetB.json").write_text('{"prior": 2}')

    mp2_path = root / "data/input/predictors/mp2_training_data.public.txt.gz"
    mp2_path.parent.mkdir(parents=True, exist_ok=True)
    mp2_path.write_bytes(b"mp2 public source")

    return "data/input/predictors/mp2_training_data.public.txt.gz"


def test_build_file_manifest_includes_every_static_file(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    keys = {e.key for e in entries}
    for _local_rel, key, _restricted in _STATIC_FILES:
        assert key in keys


def test_build_file_manifest_skips_dotfiles_in_glob_dirs(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    keys = {e.key for e in entries}
    assert "external_data_supplements__a.xlsx" in keys
    assert "external_data_supplements__b.pdf" in keys
    assert not any(k.endswith("DS_Store") for k in keys)


def test_build_file_manifest_zips_excalibr_json_directory(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    zip_entry = next(e for e in entries if e.key == _EXCALIBR_JSON_ZIP_KEY)
    with zipfile.ZipFile(zip_entry.local_path) as zf:
        assert sorted(zf.namelist()) == ["DatasetA.json", "DatasetB.json"]


def test_build_file_manifest_includes_mp2_public_source(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    assert any(e.key == MP2_PUBLIC_SOURCE_KEY for e in entries)


def test_build_file_manifest_marks_restricted_files(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    restricted_keys = {e.key for e in entries if e.restricted}
    assert "predictors__revel_hg38.tsv.gz" in restricted_keys
    assert "gnomAD__gnomad_with_ancestry_af.csv" in restricted_keys
    assert "maves__score_sets.tsv" not in restricted_keys


def test_build_file_manifest_keys_never_contain_a_slash(tmp_path):
    """Zenodo's legacy bucket upload 404s on a key containing "/" -- regression guard."""
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    assert not any("/" in entry.key for entry in entries)


def test_build_file_manifest_raises_on_missing_file(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    (tmp_path / "data/input/maves/score_sets.tsv").unlink()
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    with pytest.raises(ValueError, match="score_sets.tsv"):
        build_file_manifest(tmp_path, mp2_source, work_dir)


def test_build_file_manifest_raises_on_missing_mp2_source(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    (tmp_path / mp2_source).unlink()
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    with pytest.raises(ValueError, match="mp2_training_data"):
        build_file_manifest(tmp_path, mp2_source, work_dir)


def test_render_manifest_summary_flags_restricted_files(tmp_path):
    mp2_source = _make_complete_repo(tmp_path)
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    entries = build_file_manifest(tmp_path, mp2_source, work_dir)

    summary = render_manifest_summary(entries)

    assert "predictors__revel_hg38.tsv.gz" in summary
    revel_line = next(line for line in summary.splitlines() if "predictors__revel_hg38.tsv.gz" in line)
    assert "NON-COMMERCIAL RESTRICTED" in revel_line
    score_sets_line = next(line for line in summary.splitlines() if "maves__score_sets.tsv" in line)
    assert "NON-COMMERCIAL RESTRICTED" not in score_sets_line


def test_build_related_identifiers_includes_preprint_and_github():
    related = build_related_identifiers("https://github.com/example/repo")

    identifiers = {r["identifier"] for r in related}
    assert PREPRINT_DOI in identifiers
    assert "https://github.com/example/repo" in identifiers
    assert all(r["relation"] == "isSupplementTo" for r in related)


def test_build_related_identifiers_falls_back_to_known_github_url():
    related = build_related_identifiers(None)

    identifiers = {r["identifier"] for r in related}
    assert FALLBACK_GITHUB_URL in identifiers
