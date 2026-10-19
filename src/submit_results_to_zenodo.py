#!/usr/bin/env python3
"""Upload the paper's result data files as a new version of their Zenodo record.

Creates a new version of the existing "Data from 'A scalable approach to
resolving variants of uncertain significance'" Zenodo record (identified by
its concept recid, DEFAULT_DEPOSITION_ID below -- see
src/lib/zenodo_client.py's docstring for why a concept recid rather than a
specific version id) and uploads the integrated variant effect dataset
(Supplementary Data 1) plus Supplementary Data 3, 5, and 6.

Unlike `submit_data_prerequisites_to_zenodo.py`, this script does not touch
the record's metadata (title, description, license, related identifiers) --
it's already correctly set from when the record was first created, and
there's nothing about it that needs to change here. This script's only job
is getting the files uploaded.

## Usage

Dry run (default -- no network calls; lists every file the manifest would
upload, with its size):

    poetry run python -m src.submit_results_to_zenodo

Real run against the sandbox first (a brand-new sandbox deposition, since
sandbox.zenodo.org doesn't have a copy of the production record):

    export ZENODO_TOKEN=...   # a sandbox.zenodo.org personal access token
    poetry run python -m src.submit_results_to_zenodo --host sandbox --new --execute

Once the sandbox draft looks right (check it in the sandbox web UI), create
a new version of the real production record:

    export ZENODO_TOKEN=...   # a zenodo.org personal access token (different from sandbox's)
    poetry run python -m src.submit_results_to_zenodo --host prod --execute

That leaves a draft for you to review at the printed URL. Only once you're
satisfied, publish it:

    poetry run python -m src.submit_results_to_zenodo --host prod --execute --publish
"""

import os
from pathlib import Path

import click

from src.lib.zenodo_client import FileEntry, ZenodoClient, check_no_slashed_keys, render_manifest_summary

ZENODO_HOSTS = {
    "sandbox": "https://sandbox.zenodo.org/api/",
    "prod": "https://zenodo.org/api/",
}

# The existing "Data from ..." record this script creates a new version of.
# This is its *concept* recid (10.5281/zenodo.23308563), confirmed via
# GET https://zenodo.org/api/records/23308563 -> {"id": 23308564, ...}
# (23308564 being the latest version as of this writing). See
# src/lib/zenodo_client.py's docstring for why a concept recid, not a
# specific version id. Only meaningful on --host prod -- sandbox has no copy
# of this record, so --new is required there.
DEFAULT_DEPOSITION_ID = "23308563"

# (local path relative to repo root, destination key in the deposit).
# All flat filenames -- no "/" to worry about, unlike the data-prerequisites
# manifest. Keep in sync with the file list printed in the PR/commit that
# introduces this script, and with docs/submit_results_to_zenodo.md.
_STATIC_FILES = [
    ("data/output/supplementary_data/Supplementary_Data_1.tsv.gz", "Supplementary_Data_1.tsv.gz"),
    ("data/output/supplementary_data/Supplementary_Data_1_schema.xlsx", "Supplementary_Data_1_schema.xlsx"),
    ("data/output/supplementary_data/Supplementary_Data_3.xlsx", "Supplementary_Data_3.xlsx"),
    ("data/output/supplementary_data/Supplementary_Data_5.xlsx", "Supplementary_Data_5.xlsx"),
    (
        "data/output/supplementary_data/Supplementary_Data_5.with_secondary_variants.xlsx",
        "Supplementary_Data_5.with_secondary_variants.xlsx",
    ),
    ("data/output/supplementary_data/Supplementary_Data_6.xlsx", "Supplementary_Data_6.xlsx"),
    (
        "data/output/supplementary_data/Supplementary_Data_6.with_secondary_variants.xlsx",
        "Supplementary_Data_6.with_secondary_variants.xlsx",
    ),
]


def build_file_manifest(repo_root):
    """Return the list of `FileEntry` this deposit should contain.

    Raises `ValueError` (wrapped as `click.ClickException` in `main`) if any
    expected file is missing -- better to fail loudly before touching the
    network than to silently publish a partial archive (the exact problem
    that prompted writing this script).
    """
    repo_root = Path(repo_root)
    entries = []
    missing = []
    for local_rel, key in _STATIC_FILES:
        path = repo_root / local_rel
        if not path.exists():
            missing.append(str(path))
            continue
        entries.append(FileEntry(path, key))

    if missing:
        raise ValueError("Missing expected input file(s):\n  " + "\n  ".join(missing))

    check_no_slashed_keys(entries)

    return entries


@click.command(help=__doc__)
@click.option(
    "--repo-root",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=Path(__file__).resolve().parent.parent,
    show_default=True,
    help="Repo checkout to read data/output/supplementary_data/ from.",
)
@click.option(
    "--host",
    type=click.Choice(sorted(ZENODO_HOSTS)),
    default="sandbox",
    show_default=True,
    help="Which Zenodo environment to talk to.",
)
@click.option(
    "--deposition-id",
    default=DEFAULT_DEPOSITION_ID,
    show_default=True,
    help=(
        "Concept recid (or any specific version's own id) of the existing record to create a new "
        "version of -- resolved to the current latest published version automatically. Ignored with --new."
    ),
)
@click.option(
    "--new",
    "create_new",
    is_flag=True,
    help=(
        "Create a brand-new deposition instead of versioning --deposition-id. Required on --host "
        "sandbox (sandbox has no copy of the production record) -- metadata is seeded from the "
        "production record's public GET, not fabricated."
    ),
)
@click.option(
    "--execute",
    is_flag=True,
    help="Actually create/update the draft and upload files. Without this, only prints the manifest.",
)
@click.option(
    "--publish",
    is_flag=True,
    help="After updating the draft, also publish it. IRREVERSIBLE. Requires --execute.",
)
@click.option("--yes", "skip_confirm", is_flag=True, help="Don't ask for interactive confirmation on --host prod.")
def main(repo_root, host, deposition_id, create_new, execute, publish, skip_confirm):
    if execute and host == "sandbox" and not create_new:
        raise click.ClickException(
            "--host sandbox has no copy of the production deposition -- pass --new to create a fresh "
            "sandbox test deposition instead of trying to version --deposition-id."
        )
    if publish and not execute:
        raise click.ClickException("--publish requires --execute.")

    try:
        entries = build_file_manifest(repo_root)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    click.echo(f"Manifest ({host}" + (", new deposition" if create_new else f", new version of {deposition_id}") + "):")
    click.echo(render_manifest_summary(entries))

    if not execute:
        click.echo("\nDry run only (pass --execute to create/update the draft).")
        return

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        raise click.ClickException("ZENODO_TOKEN must be set in the environment to --execute.")

    if host == "prod" and not skip_confirm:
        action = "publish a new version of" if publish else "create/update a draft new version of"
        click.confirm(
            f"This will {action} PRODUCTION record {deposition_id} (resolved to its current latest version). Continue?",
            abort=True,
        )

    client = ZenodoClient(ZENODO_HOSTS[host], token)

    if create_new:
        click.echo("\nSeeding metadata from the production record's public listing...")
        base_metadata = client.get_public_metadata(DEFAULT_DEPOSITION_ID)
        draft = client.create_new_deposition()
        draft_id = draft["id"]
        click.echo(f"Updating metadata on new sandbox draft {draft_id}...")
        client.update_metadata(draft_id, base_metadata)
    else:
        click.echo(f"\nResolving {deposition_id} to its current latest published version...")
        latest_id = client.resolve_latest_record_id(deposition_id)
        click.echo(f"Creating new version of {latest_id}...")
        draft = client.create_new_version(latest_id)
        draft_id = draft["id"]

    bucket_url = draft["links"]["bucket"]

    existing_files = client.list_files(draft_id)
    manifest_keys = {entry.key for entry in entries}
    for f in existing_files:
        if f["filename"] in manifest_keys:
            click.echo(f"Removing existing {f['filename']} before re-upload...")
            client.delete_file(draft_id, f["id"])

    for entry in entries:
        click.echo(f"Uploading {entry.key} ({entry.local_path.stat().st_size:,} bytes)...")
        client.upload_file(bucket_url, entry)

    if not publish:
        click.echo(f"\nDraft ready for review: {draft['links']['html']}")
        click.echo("Not publishing (pass --publish to go live -- irreversible).")
        return

    click.echo("\nPublishing...")
    published = client.publish(draft_id)
    click.echo(f"Published: {published['links']['record_html']}")


if __name__ == "__main__":
    main()
