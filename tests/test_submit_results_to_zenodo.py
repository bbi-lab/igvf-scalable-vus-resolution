import pytest

from src.submit_results_to_zenodo import _STATIC_FILES, build_file_manifest


def _make_complete_repo(root):
    """Create every file `build_file_manifest` expects, each with distinct content."""
    for local_rel, _key in _STATIC_FILES:
        path = root / local_rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(local_rel.encode())


def test_build_file_manifest_includes_every_static_file(tmp_path):
    _make_complete_repo(tmp_path)

    entries = build_file_manifest(tmp_path)

    keys = {e.key for e in entries}
    for _local_rel, key in _STATIC_FILES:
        assert key in keys


def test_build_file_manifest_keys_are_flat_filenames(tmp_path):
    _make_complete_repo(tmp_path)

    entries = build_file_manifest(tmp_path)

    assert not any("/" in entry.key for entry in entries)


def test_build_file_manifest_raises_on_missing_file(tmp_path):
    _make_complete_repo(tmp_path)
    (tmp_path / "data/output/supplementary_data/Supplementary_Data_6.xlsx").unlink()

    with pytest.raises(ValueError, match="Supplementary_Data_6.xlsx"):
        build_file_manifest(tmp_path)
