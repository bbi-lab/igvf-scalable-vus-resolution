#!/usr/bin/env python3
"""Submit an updated Supplementary Data 1 to the IGVF data portal.

Posts a new `TabularFile` object (the updated
`data/output/supplementary_data/Supplementary_Data_1.tsv.gz`) and a new
`Document` object (its column-description/format-specification file) to
the IGVF portal, each with a `supersedes` link back to the object it
replaces:

- TabularFile `IGVFFI2521UGYG` (https://data.igvf.org/tabular-files/IGVFFI2521UGYG/)
- Document `2bfe67c5-937a-4fe0-b343-5008bdbf7241`
  (https://data.igvf.org/documents/2bfe67c5-937a-4fe0-b343-5008bdbf7241/)

The portal schema makes `status` (and the `replaced`/`revoked`/`archived`
transitions that implement "replacement") an admin-only, read-only field --
submitters cannot flip it themselves. This script can only create the new
objects and link them via `supersedes`; an IGVF DACC curator must still
transition the two objects above to `revoked`/`archived` status and release
the new ones. A summary of the created objects is printed at the end so
that request can be made.

The format-specification document must be a flat TSV (the portal's
`Document.attachment` MIME-type allowlist has no spreadsheet type), so if
`--format-spec` is an `.xlsx` file, it's first converted to TSV with pandas.

`content_type: "variants"` also triggers a schema rule requiring one of
`reference_files`, `submitted_assembly`, or
`submitted_transcriptome_annotation` -- the latter two are admin-only, so
this script links the GRCh38 `ReferenceFile`
(https://data.igvf.org/reference-files/IGVFFI0653VCGH/) via
`reference_files` by default, which is a plain submittable field.

## Usage

Dry run (default -- no network calls; builds both payloads, converts the
format spec, computes the new file's md5sum/size, and prints everything for
review):

    poetry run python -m src.submit_supplementary_data_1_to_igvf

Real submission, against staging first:

    export IGVF_API_KEY=...
    export IGVF_SECRET_KEY=...
    poetry run python -m src.submit_supplementary_data_1_to_igvf --host staging --execute

Then, once the staging result looks right, against production:

    poetry run python -m src.submit_supplementary_data_1_to_igvf --host prod --execute

`--execute --host prod` asks for an interactive y/N confirmation before
doing anything, unless `--yes` is also passed. Requires the `aws` CLI on
PATH (used to upload the file to the credentials the portal issues).
"""

import base64
import hashlib
import mimetypes
import os
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import click
import pandas as pd
import requests

IGVF_HOSTS = {
    # "sandbox" (api.sandbox.igvf.org) is deprecated by the portal in favor of "staging".
    "staging": "https://api.staging.igvf.org/",
    "prod": "https://api.data.igvf.org/",
}

EXISTING_TABULAR_FILE_ID = "/tabular-files/IGVFFI2521UGYG/"
EXISTING_DOCUMENT_ID = "/documents/2bfe67c5-937a-4fe0-b343-5008bdbf7241/"

LAB = "/labs/lea-starita/"
AWARD = "/awards/HG011969/"
FILE_SET = "/curated-sets/IGVFDS4011IBUI/"
CONTENT_TYPE = "variants"
FILE_FORMAT = "tsv"
DOCUMENT_TYPE = "file format specification"
DOCUMENT_DESCRIPTION = (
    "This document provides column definitions for the IGVF Coding Focus Group Pillar Project dataframe."
)
ATTACHMENT_MIME_TYPE = "text/tab-separated-values"

DEFAULT_DATA_FILE = Path("data/output/supplementary_data/Supplementary_Data_1.tsv.gz")
DEFAULT_FORMAT_SPEC = Path("data/output/supplementary_data/Supplementary_Data_1_schema.xlsx")
DEFAULT_SUBMITTED_FILE_NAME = "integrated_variant_effect_dataset.tsv.gz"
# GRCh38 ReferenceFile: https://data.igvf.org/reference-files/IGVFFI0653VCGH/
DEFAULT_REFERENCE_FILES = ("/reference-files/IGVFFI0653VCGH/",)


def md5sum(path):
    """Return the hex md5 digest of the file at `path`, read in chunks."""
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def xlsx_schema_to_tsv(xlsx_path, tsv_path):
    """Convert a single-sheet column-description `xlsx_path` to a TSV at `tsv_path`.

    The portal's `Document.attachment` MIME-type allowlist doesn't include a
    spreadsheet type, so the format-specification document must be
    submitted as plain text. Column headers are stripped of surrounding
    whitespace (the source sheet has a stray trailing space on one header).
    """
    df = pd.read_excel(xlsx_path)
    df.columns = [str(c).strip() for c in df.columns]
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)


def build_attachment(tsv_path, mime_type=ATTACHMENT_MIME_TYPE):
    """Return a `Document.attachment` dict (base64 data URI) for the file at `tsv_path`."""
    content = Path(tsv_path).read_bytes()
    encoded = base64.b64encode(content).decode("ascii")
    return {
        "download": Path(tsv_path).name,
        "type": mime_type,
        "href": f"data:{mime_type};base64,{encoded}",
    }


def build_document_payload(attachment, alias, supersedes=(EXISTING_DOCUMENT_ID,)):
    """Return the POST payload for the new format-specification Document.

    `standardized_file_format` is deliberately omitted: despite appearing on
    the existing document, it's an admin-only property (the portal rejects
    it from a submitter with a 422, even when set to the same value an
    admin previously set) and isn't in the schema's `required` list.
    """
    return {
        "lab": LAB,
        "award": AWARD,
        "document_type": DOCUMENT_TYPE,
        "description": DOCUMENT_DESCRIPTION,
        "aliases": [alias],
        "attachment": attachment,
        "supersedes": list(supersedes),
    }


def build_tabular_file_payload(
    md5,
    file_size,
    aliases,
    submitted_file_name,
    file_format_specification_id,
    supersedes=(EXISTING_TABULAR_FILE_ID,),
    reference_files=DEFAULT_REFERENCE_FILES,
):
    """Return the POST payload for the new Supplementary Data 1 TabularFile.

    `content_type: "variants"` triggers a schema rule requiring one of
    `reference_files`, `submitted_assembly`, or
    `submitted_transcriptome_annotation` to be present. The latter two are
    admin-only, so `reference_files` -- linking the GRCh38 `ReferenceFile`
    by default -- is the field a lab submitter can actually set; pass
    `reference_files=None` to omit it.
    """
    payload = {
        "lab": LAB,
        "award": AWARD,
        "file_set": FILE_SET,
        "content_type": CONTENT_TYPE,
        "file_format": FILE_FORMAT,
        "controlled_access": False,
        "md5sum": md5,
        "file_size": file_size,
        "aliases": list(aliases),
        "submitted_file_name": submitted_file_name,
        "file_format_specifications": [file_format_specification_id],
        "supersedes": list(supersedes),
    }
    if reference_files:
        payload["reference_files"] = list(reference_files)
    return payload


class IgvfClient:
    """Thin wrapper over the IGVF REST API: HTTP Basic auth, POST, and S3 upload credentials."""

    def __init__(self, host_url, api_key, secret_key):
        self.host_url = host_url.rstrip("/") + "/"
        self.session = requests.Session()
        self.session.auth = (api_key, secret_key)
        self.session.headers.update({"Content-Type": "application/json", "Accept": "application/json"})

    def post(self, collection, payload):
        response = self.session.post(self.host_url + collection + "/", json=payload)
        if not response.ok:
            raise click.ClickException(f"POST {collection}/ failed with {response.status_code}: {response.text}")
        return response.json()["@graph"][0]

    def get_upload_credentials(self, accession):
        response = self.session.post(self.host_url + "files/" + accession + "/@@upload", json={})
        if not response.ok:
            raise click.ClickException(
                f"Requesting upload credentials for {accession} failed with {response.status_code}: {response.text}"
            )
        return response.json()["@graph"][0]["upload_credentials"]


def upload_file_to_s3(local_path, upload_credentials):
    """Upload `local_path` to the S3 destination in `upload_credentials` via `aws s3 cp`."""
    if shutil.which("aws") is None:
        raise click.ClickException("The `aws` CLI is required to upload the file but was not found on PATH.")
    env = dict(os.environ)
    env["AWS_ACCESS_KEY_ID"] = upload_credentials["access_key"]
    env["AWS_SECRET_ACCESS_KEY"] = upload_credentials["secret_key"]
    env["AWS_SESSION_TOKEN"] = upload_credentials["session_token"]
    upload_url = upload_credentials["upload_url"]
    result = subprocess.run(["aws", "s3", "cp", str(local_path), upload_url], env=env, check=False)
    if result.returncode != 0:
        raise click.ClickException(f"`aws s3 cp` exited with status {result.returncode}.")


@click.command(help=__doc__)
@click.option(
    "--data-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_DATA_FILE,
    show_default=True,
    help="Updated Supplementary Data 1 TSV(.gz) to submit as the new TabularFile.",
)
@click.option(
    "--format-spec",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_FORMAT_SPEC,
    show_default=True,
    help="Column-description file for --data-file (.xlsx is converted to TSV automatically).",
)
@click.option(
    "--tabular-file-alias",
    default=lambda: f"lea-starita:integrated-variant-effect-dataset-{datetime.now(UTC):%Y%m%d}",
    show_default=True,
    help="Unique portal alias for the new TabularFile.",
)
@click.option(
    "--document-alias",
    default=lambda: f"lea-starita:variant-effect-dataset-specifications-{datetime.now(UTC):%Y%m%d}",
    show_default=True,
    help="Unique portal alias for the new Document.",
)
@click.option(
    "--attachment-filename",
    default=lambda: f"integrated_variant_effect_dataset_ffs_{datetime.now(UTC):%Y%m%d}.tsv",
    show_default=True,
    help="Filename recorded for the Document's attachment (no path -- only the basename is used).",
)
@click.option(
    "--submitted-file-name",
    default=DEFAULT_SUBMITTED_FILE_NAME,
    show_default=True,
    help="Filename recorded as the TabularFile's submitted_file_name (no path -- only the basename is used).",
)
@click.option(
    "--reference-file",
    "reference_files",
    multiple=True,
    default=DEFAULT_REFERENCE_FILES,
    show_default=True,
    help=(
        "ReferenceFile @id/accession to link via the TabularFile's reference_files, required by the "
        "portal schema for content_type='variants'. Repeat to link more than one."
    ),
)
@click.option(
    "--existing-document-id",
    default=None,
    help=(
        "Skip creating a new Document and use this @id (e.g. /documents/<uuid>/) as the TabularFile's "
        "file_format_specifications instead. Use this to retry after the Document was already created "
        "in a previous run -- re-running without it would try to create a second Document."
    ),
)
@click.option(
    "--upload-only-accession",
    default=None,
    help=(
        "Skip creating a Document and TabularFile entirely and just upload --data-file to this "
        "already-created TabularFile accession. Use this to retry the upload step on its own, e.g. "
        "after the metadata POST succeeded but the upload itself failed."
    ),
)
@click.option(
    "--host",
    type=click.Choice(sorted(IGVF_HOSTS)),
    default="staging",
    show_default=True,
    help="Which IGVF portal environment to submit to.",
)
@click.option(
    "--execute",
    is_flag=True,
    help="Actually make the API calls and upload the file. Without this, only builds and prints the payloads.",
)
@click.option(
    "--yes",
    "skip_confirm",
    is_flag=True,
    help="Don't ask for interactive confirmation before submitting to --host prod.",
)
def main(
    data_file,
    format_spec,
    tabular_file_alias,
    document_alias,
    attachment_filename,
    submitted_file_name,
    reference_files,
    existing_document_id,
    upload_only_accession,
    host,
    execute,
    skip_confirm,
):
    # Basename only: the portal records these as metadata, not a path on this machine.
    attachment_filename = Path(attachment_filename).name
    submitted_file_name = Path(submitted_file_name).name

    click.echo(f"Data file:   {data_file} ({data_file.stat().st_size:,} bytes)")
    click.echo("Computing md5sum...")
    md5 = md5sum(data_file)
    file_size = data_file.stat().st_size

    if upload_only_accession:
        click.echo(f"md5sum: {md5}")
        if not execute:
            click.echo(
                f"\nDry run only (pass --execute to upload). Would upload {data_file} to the TabularFile "
                f"{upload_only_accession} on {host}."
            )
            return
        api_key = os.environ.get("IGVF_API_KEY")
        secret_key = os.environ.get("IGVF_SECRET_KEY")
        if not api_key or not secret_key:
            raise click.ClickException("IGVF_API_KEY and IGVF_SECRET_KEY must be set in the environment to upload.")
        if host == "prod" and not skip_confirm:
            click.confirm(
                f"This will upload {data_file} to PRODUCTION accession {upload_only_accession}. Continue?", abort=True
            )
        client = IgvfClient(IGVF_HOSTS[host], api_key, secret_key)
        click.echo(f"Requesting upload credentials for {upload_only_accession}...")
        upload_credentials = client.get_upload_credentials(upload_only_accession)
        click.echo(f"Uploading {data_file} to {upload_credentials['upload_url']}...")
        upload_file_to_s3(data_file, upload_credentials)
        click.echo(f"\nDone. Uploaded to {upload_only_accession}.")
        return

    document_payload = None
    attachment = None
    if existing_document_id:
        click.echo(f"\nReusing existing Document {existing_document_id}; not creating a new one.")
    else:
        if format_spec.suffix.lower() in (".xlsx", ".xls"):
            tsv_spec_path = format_spec.with_name(attachment_filename)
            click.echo(f"Converting {format_spec} -> {tsv_spec_path}")
            xlsx_schema_to_tsv(format_spec, tsv_spec_path)
        else:
            tsv_spec_path = format_spec

        attachment = build_attachment(
            tsv_spec_path, mime_type=mimetypes.guess_type(attachment_filename)[0] or ATTACHMENT_MIME_TYPE
        )
        document_payload = build_document_payload(attachment, alias=document_alias)
        click.echo("\nDocument payload (attachment content omitted):")
        preview = {**document_payload, "attachment": {k: v for k, v in attachment.items() if k != "href"}}
        click.echo(preview)

    # file_format_specifications is filled in with the real document @id once it's posted (or immediately,
    # if --existing-document-id was given).
    tabular_file_payload = build_tabular_file_payload(
        md5=md5,
        file_size=file_size,
        aliases=[tabular_file_alias],
        submitted_file_name=submitted_file_name,
        file_format_specification_id=existing_document_id or "<pending: new document @id>",
        reference_files=reference_files or None,
    )
    click.echo("\nTabularFile payload:")
    click.echo(tabular_file_payload)

    if not execute:
        click.echo(
            "\nDry run only (pass --execute to submit). An IGVF DACC curator must still transition "
            f"{EXISTING_DOCUMENT_ID} to 'archived' and {EXISTING_TABULAR_FILE_ID} to 'revoked', and "
            "release the new objects -- submitters can't set status themselves."
        )
        return

    api_key = os.environ.get("IGVF_API_KEY")
    secret_key = os.environ.get("IGVF_SECRET_KEY")
    if not api_key or not secret_key:
        raise click.ClickException("IGVF_API_KEY and IGVF_SECRET_KEY must be set in the environment to --execute.")

    if host == "prod" and not skip_confirm:
        click.confirm("This will create new objects on the PRODUCTION IGVF portal. Continue?", abort=True)

    client = IgvfClient(IGVF_HOSTS[host], api_key, secret_key)

    if existing_document_id:
        new_document_id = existing_document_id
    else:
        click.echo(f"\nPosting Document to {host}...")
        new_document = client.post("documents", document_payload)
        new_document_id = new_document["@id"]
        click.echo(f"Created {new_document_id}")

    tabular_file_payload["file_format_specifications"] = [new_document_id]
    click.echo(f"Posting TabularFile to {host}...")
    new_tabular_file = client.post("tabular-files", tabular_file_payload)
    click.echo(f"Created {new_tabular_file['@id']}")

    click.echo("Requesting upload credentials...")
    upload_credentials = client.get_upload_credentials(new_tabular_file["accession"])
    click.echo(f"Uploading {data_file} to {upload_credentials['upload_url']}...")
    upload_file_to_s3(data_file, upload_credentials)

    click.echo("\nDone. Summary:")
    click.echo(f"  Document:         {new_document_id} (supersedes {EXISTING_DOCUMENT_ID})")
    click.echo(f"  New TabularFile:  {new_tabular_file['@id']} (supersedes {EXISTING_TABULAR_FILE_ID})")
    click.echo(
        "\nA DACC curator must still transition the old objects to 'archived'/'revoked' status and "
        "release the new ones -- those fields are admin-only."
    )


if __name__ == "__main__":
    main()
