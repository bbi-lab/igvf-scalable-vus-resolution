import subprocess

import pytest

from src.lib.zenodo_client import (
    FileEntry,
    ZenodoClient,
    check_no_slashed_keys,
    github_https_url,
    render_manifest_summary,
)


class _FakeResponse:
    def __init__(self, json_body, ok=True, status_code=200):
        self._json_body = json_body
        self.ok = ok
        self.status_code = status_code
        self.text = str(json_body)

    def json(self):
        return self._json_body


def test_github_https_url_converts_ssh_remote(tmp_path):
    subprocess.run(["git", "init"], cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(
        ["git", "remote", "add", "origin", "git@github.com:bbi-lab/igvf-scalable-vus-resolution.git"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )

    assert github_https_url(tmp_path) == "https://github.com/bbi-lab/igvf-scalable-vus-resolution"


def test_github_https_url_returns_none_without_git_repo(tmp_path):
    assert github_https_url(tmp_path) is None


def test_check_no_slashed_keys_passes_when_none_contain_a_slash():
    entries = [FileEntry("a", "maves__score_sets.tsv"), FileEntry("b", "LICENSE.md")]
    check_no_slashed_keys(entries)  # should not raise


def test_check_no_slashed_keys_raises_when_a_key_contains_a_slash():
    entries = [FileEntry("a", "maves/score_sets.tsv")]
    with pytest.raises(ValueError, match="maves/score_sets.tsv"):
        check_no_slashed_keys(entries)


def test_render_manifest_summary_flags_restricted_files(tmp_path):
    f1 = tmp_path / "a.txt"
    f1.write_bytes(b"x" * 10)
    f2 = tmp_path / "b.txt"
    f2.write_bytes(b"y" * 20)
    entries = [FileEntry(f1, "a.txt", restricted=True), FileEntry(f2, "b.txt", restricted=False)]

    summary = render_manifest_summary(entries)

    a_line = next(line for line in summary.splitlines() if "a.txt" in line)
    assert "NON-COMMERCIAL RESTRICTED" in a_line
    b_line = next(line for line in summary.splitlines() if " b.txt" in line)
    assert "NON-COMMERCIAL RESTRICTED" not in b_line


def test_resolve_latest_record_id_returns_the_latest_versions_own_id(monkeypatch):
    # A concept recid's public record GET resolves to the latest version's own id --
    # confirmed against production.
    captured_url = {}

    def fake_get(url):
        captured_url["url"] = url
        return _FakeResponse({"id": 18637474, "conceptrecid": "18637473"})

    monkeypatch.setattr("requests.get", fake_get)

    client = ZenodoClient("https://zenodo.org/api/", token="unused")
    latest_id = client.resolve_latest_record_id("18637473")

    assert latest_id == 18637474
    assert captured_url["url"] == "https://zenodo.org/api/records/18637473"


def test_resolve_latest_record_id_raises_on_error_response(monkeypatch):
    monkeypatch.setattr("requests.get", lambda url: _FakeResponse({"message": "not found"}, ok=False, status_code=404))

    client = ZenodoClient("https://zenodo.org/api/", token="unused")

    with pytest.raises(Exception, match="404"):
        client.resolve_latest_record_id("18637473")


def test_get_public_metadata_returns_just_the_metadata_field(monkeypatch):
    monkeypatch.setattr(
        "requests.get",
        lambda url: _FakeResponse({"id": 18637474, "metadata": {"title": "Example"}}),
    )

    client = ZenodoClient("https://zenodo.org/api/", token="unused")

    assert client.get_public_metadata("18637473") == {"title": "Example"}


def test_create_new_version_uses_the_versions_endpoint_not_the_legacy_action(monkeypatch):
    # Regression test: the legacy POST .../actions/newversion 400s with
    # "Please remove all files first" whenever an unpublished draft already
    # exists for the record -- create_new_version must use
    # POST /records/<id>/versions instead (see its docstring).
    post_urls = []
    get_urls = []

    def fake_post(url, json=None):
        post_urls.append(url)
        assert url == "https://zenodo.org/api/records/18637474/versions"
        return _FakeResponse({"id": 99999999})

    def fake_get(url):
        get_urls.append(url)
        assert url == "https://zenodo.org/api/deposit/depositions/99999999"
        return _FakeResponse({"id": 99999999, "submitted": False, "links": {"bucket": "https://example/bucket"}})

    client = ZenodoClient("https://zenodo.org/api/", token="unused")
    monkeypatch.setattr(client.session, "post", fake_post)
    monkeypatch.setattr(client.session, "get", fake_get)

    draft = client.create_new_version("18637474")

    assert draft["id"] == 99999999
    assert post_urls == ["https://zenodo.org/api/records/18637474/versions"]
    assert get_urls == ["https://zenodo.org/api/deposit/depositions/99999999"]


def test_create_new_version_raises_if_draft_is_not_a_distinct_unpublished_version(monkeypatch):
    client = ZenodoClient("https://zenodo.org/api/", token="unused")
    monkeypatch.setattr(client.session, "post", lambda url, json=None: _FakeResponse({"id": 18637474}))
    # Zenodo handed back the *original* record, not a new draft -- something went wrong.
    monkeypatch.setattr(client.session, "get", lambda url: _FakeResponse({"id": 18637474, "submitted": True}))

    with pytest.raises(Exception, match="unpublished new-version draft"):
        client.create_new_version("18637474")
