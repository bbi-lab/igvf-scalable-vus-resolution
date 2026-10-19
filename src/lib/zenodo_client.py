"""Shared Zenodo REST API client, file-entry type, and small helpers.

Used by `src/submit_data_prerequisites_to_zenodo.py` and
`src/submit_results_to_zenodo.py`. Factored out here rather than duplicated
because the non-obvious parts -- how to create/resume a new version
without hitting known Zenodo API bugs, and how file keys have to be shaped
-- apply identically to both scripts and should not be allowed to drift.

Three Zenodo API quirks this module works around, confirmed against
production rather than assumed:

1. A deposit file key containing a literal "/" 404s on upload (the legacy
   bucket-upload route doesn't accept a multi-segment path) -- callers must
   keep keys flat (e.g. join subdirectory components with "__" instead).
2. `POST /deposit/depositions/<id>/actions/newversion` requires the
   *specific latest version's own* id, not a concept recid, and 400s with
   `{"field": "files.enabled", "messages": ["Please remove all files
   first."]}` whenever an unpublished new-version draft already exists for
   the record (e.g. left over from a previous run that failed partway
   through, or that simply timed out on Zenodo's end after creating it) --
   a known Zenodo issue, not specific to any one project.
3. Both problems above are avoided by: (a) resolving a stable *concept*
   recid to the current latest version's id via the public records API
   before versioning (`resolve_latest_record_id`), and (b) creating/resuming
   the new-version draft via `POST /records/<id>/versions` instead of the
   legacy `actions/newversion` action (`ZenodoClient.create_new_version`).
"""

import subprocess
from pathlib import Path

import click
import requests


class FileEntry:
    def __init__(self, local_path, key, restricted=False):
        self.local_path = Path(local_path)
        self.key = key
        self.restricted = restricted


class ZenodoClient:
    """Thin wrapper over the Zenodo REST API (legacy `/api/deposit/depositions`)."""

    def __init__(self, host_url, token):
        self.host_url = host_url.rstrip("/") + "/"
        self.session = requests.Session()
        self.session.headers.update({"Authorization": f"Bearer {token}"})

    def _check(self, response, what):
        if not response.ok:
            raise click.ClickException(f"{what} failed with {response.status_code}: {response.text}")
        return response

    def get_public_record(self, record_id, public_host_url="https://zenodo.org/api/"):
        """GET the *published*, public record for `record_id` -- no auth needed or sent.

        `record_id` may be a concept recid or a specific version's own id --
        Zenodo's public records API resolves either to the record representing
        the *current latest* published version (confirmed against production:
        GET /records/<concept recid> and GET /records/<latest version's own id>
        both return that latest version's own `id`).
        """
        response = requests.get(public_host_url.rstrip("/") + f"/records/{record_id}")
        self._check(response, f"GET public record {record_id}")
        return response.json()

    def get_public_metadata(self, record_id, public_host_url="https://zenodo.org/api/"):
        """GET just the *published*, public metadata for `record_id` (see `get_public_record`).

        Used to seed sandbox test depositions with realistic metadata
        (creators, affiliations, etc.) without fabricating any of it.
        """
        return self.get_public_record(record_id, public_host_url)["metadata"]

    def resolve_latest_record_id(self, record_id, public_host_url="https://zenodo.org/api/"):
        """Resolve `record_id` (concept recid or any version's own id) to the
        current latest published version's own id -- the only id
        `create_new_version` below will accept (Zenodo's `newversion` action
        requires the specific latest version's id, not the concept id)."""
        return self.get_public_record(record_id, public_host_url)["id"]

    def create_new_version(self, deposition_id):
        """Create (or resume) an unpublished new-version draft of `deposition_id`.

        Deliberately uses the newer `POST /records/<id>/versions` endpoint
        rather than the legacy `POST /deposit/depositions/<id>/actions/newversion`
        action -- see this module's docstring, quirk 2. The newer endpoint
        resumes an existing draft instead of erroring. Its response id is in
        the same numbering space as the legacy deposit API, so it's
        re-fetched through `deposit/depositions/<id>` to get the familiar
        representation (`links.bucket`, `metadata`, etc.) callers expect.
        """
        response = self.session.post(self.host_url + f"records/{deposition_id}/versions", json={})
        self._check(response, f"POST versions for {deposition_id}")
        new_id = response.json()["id"]
        draft = self._check(
            self.session.get(self.host_url + f"deposit/depositions/{new_id}"), f"GET new draft {new_id}"
        ).json()
        if draft.get("submitted") is not False or str(draft["id"]) == str(deposition_id):
            raise click.ClickException(
                f"Expected an unpublished new-version draft distinct from {deposition_id}, got: {draft}"
            )
        return draft

    def create_new_deposition(self):
        response = self.session.post(self.host_url + "deposit/depositions", json={})
        return self._check(response, "POST new deposition").json()

    def update_metadata(self, deposition_id, metadata):
        response = self.session.put(self.host_url + f"deposit/depositions/{deposition_id}", json={"metadata": metadata})
        return self._check(response, f"PUT metadata for {deposition_id}").json()

    def list_files(self, deposition_id):
        response = self.session.get(self.host_url + f"deposit/depositions/{deposition_id}/files")
        return self._check(response, f"GET files for {deposition_id}").json()

    def delete_file(self, deposition_id, file_id):
        response = self.session.delete(self.host_url + f"deposit/depositions/{deposition_id}/files/{file_id}")
        self._check(response, f"DELETE file {file_id} from {deposition_id}")

    def upload_file(self, bucket_url, entry):
        url = bucket_url.rstrip("/") + "/" + entry.key
        with open(entry.local_path, "rb") as f:
            response = self.session.put(url, data=f)
        self._check(response, f"PUT {url}")

    def publish(self, deposition_id):
        response = self.session.post(self.host_url + f"deposit/depositions/{deposition_id}/actions/publish")
        return self._check(response, f"POST publish for {deposition_id}").json()


def render_manifest_summary(entries):
    lines = []
    total_bytes = 0
    for entry in entries:
        size = entry.local_path.stat().st_size
        total_bytes += size
        flag = "  [NON-COMMERCIAL RESTRICTED]" if entry.restricted else ""
        lines.append(f"  {entry.key:<60} {size:>14,} bytes{flag}")
    lines.append(f"\n  Total: {total_bytes:,} bytes ({total_bytes / 1e9:.2f} GB) across {len(entries)} files")
    return "\n".join(lines)


def check_no_slashed_keys(entries):
    """Raise `ValueError` if any entry's key contains "/" -- see this module's
    docstring, quirk 1. Call this at the end of a manifest-building function."""
    slashed = [entry.key for entry in entries if "/" in entry.key]
    if slashed:
        raise ValueError(f"Manifest key(s) contain '/', which Zenodo's bucket upload rejects: {slashed}")


def github_https_url(repo_root):
    """Best-effort `https://github.com/...` URL for this checkout's `origin` remote.

    Returns None (rather than guessing) if `git remote get-url origin` fails
    or isn't a recognizable GitHub remote.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "remote", "get-url", "origin"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    url = result.stdout.strip()
    if url.startswith("git@github.com:"):
        return "https://github.com/" + url[len("git@github.com:") :].removesuffix(".git")
    if url.startswith("https://github.com/"):
        return url.removesuffix(".git")
    return None
