# Submit Supplementary Data 1 to IGVF

`src/submit_supplementary_data_1_to_igvf.py` submits the updated
Supplementary Data 1 (`data/output/supplementary_data/Supplementary_Data_1.tsv.gz`)
and its column-description document to the IGVF data portal as new, versioned
objects, replacing:

- TabularFile [`IGVFFI2521UGYG`](https://data.igvf.org/tabular-files/IGVFFI2521UGYG/)
- Document [`2bfe67c5-937a-4fe0-b343-5008bdbf7241`](https://data.igvf.org/documents/2bfe67c5-937a-4fe0-b343-5008bdbf7241/)

## How "replacement" works on the portal

IGVF's `File` and `Document` schemas don't let a submitter overwrite an
already-released object in place. Instead, the pattern is:

1. POST a new object with a `supersedes` link back to the one it replaces.
2. A DACC curator transitions the old object's `status` to `revoked`
   (File) or `archived` (Document), and releases the new one.

`status` is an admin-only, read-only field in both schemas, so step 2 can't
be done by this script or by any submitter credential -- **after running
this script, ask IGVF DACC curators to complete the status transition** for
the two objects above.

The format-specification document must be plain text: the portal's
`Document.attachment` MIME-type allowlist has no spreadsheet type. If
`--format-spec` is an `.xlsx` file (the default,
`data/output/supplementary_data/Supplementary_Data_1_schema.xlsx`), the
script converts it to TSV first.

The filenames recorded on the portal (the TabularFile's
`submitted_file_name` and the Document attachment's filename) match the
original submission's naming -- `integrated_variant_effect_dataset.tsv.gz`
and `integrated_variant_effect_dataset_ffs_<date>.tsv` -- rather than this
repo's local `Supplementary_Data_1*` filenames, and are recorded as bare
filenames with no local path.

### `reference_files`

The portal schema requires a `content_type: "variants"` TabularFile to have
one of `reference_files`, `submitted_assembly`, or
`submitted_transcriptome_annotation`. The latter two are admin-only (a
lab-submitter account gets a 422/`permission 'admin_only' required` trying
to set either), but `reference_files` is a plain submittable field, so the
script links the GRCh38 `ReferenceFile`
([`IGVFFI0653VCGH`](https://data.igvf.org/reference-files/IGVFFI0653VCGH/))
by default. Pass `--reference-file` (repeatable) to link a different or
additional `ReferenceFile`.

## Prerequisites

This script is an internal submission utility for the project team -- it
isn't needed to reproduce any published result. Its Python dependencies
(`click`, `pandas`, `requests`) are already part of the project's Poetry
environment, so a regular `poetry install` is all that's needed on that
front.

- `IGVF_API_KEY` / `IGVF_SECRET_KEY` environment variables, for a portal
  account with submit permission for the `lea-starita` lab / `HG011969`
  award. Required only for `--execute`.
- The `aws` CLI on `PATH` -- the portal issues temporary S3 credentials and
  the upload itself is a plain `aws s3 cp`. This is a standalone system
  tool, not a Python package, so it's intentionally not declared in
  `pyproject.toml`; install it separately (e.g. `brew install awscli` or
  the official AWS installer). No local AWS profile/credentials are
  needed -- the script passes the portal's temporary upload credentials to
  `aws s3 cp` via environment variables for each call.

## Usage

Dry run (default): builds both payloads, converts the format spec, computes
the new file's md5sum, and prints everything for review. No network calls.

```bash
poetry run python -m src.submit_supplementary_data_1_to_igvf
```

Submit to staging first:

```bash
export IGVF_API_KEY=...
export IGVF_SECRET_KEY=...
poetry run python -m src.submit_supplementary_data_1_to_igvf --host staging --execute
```

Once the staging result looks right, submit to production (asks for
interactive confirmation unless `--yes` is also passed):

```bash
poetry run python -m src.submit_supplementary_data_1_to_igvf --host prod --execute
```

### Options

| Option | Default | Purpose |
|---|---|---|
| `--data-file` | `data/output/supplementary_data/Supplementary_Data_1.tsv.gz` | Local file to read and upload as the new TabularFile. |
| `--format-spec` | `data/output/supplementary_data/Supplementary_Data_1_schema.xlsx` | Column-description file; `.xlsx` is auto-converted to TSV. |
| `--tabular-file-alias` | `lea-starita:integrated-variant-effect-dataset-<today>` | Must be unique on the portal. |
| `--document-alias` | `lea-starita:variant-effect-dataset-specifications-<today>` | Must be unique on the portal. |
| `--attachment-filename` | `integrated_variant_effect_dataset_ffs_<today>.tsv` | Basename recorded on the Document's attachment -- matches the original submission's naming. |
| `--submitted-file-name` | `integrated_variant_effect_dataset.tsv.gz` | Basename recorded as the TabularFile's `submitted_file_name` -- matches the original submission's naming, not the local `--data-file` path. |
| `--reference-file` | `/reference-files/IGVFFI0653VCGH/` (GRCh38) | TabularFile's `reference_files`; required by the schema for `content_type: "variants"` (see above). Repeatable. |
| `--existing-document-id` | none | Skip creating a new Document and reuse this `@id` instead -- use when retrying after the Document already succeeded in a prior run (see below). |
| `--upload-only-accession` | none | Skip creating a Document/TabularFile entirely; just upload `--data-file` to this already-created TabularFile accession -- useful for retrying just the upload step. |
| `--host` | `staging` | `staging` or `prod` (the portal deprecated the old `sandbox` host in favor of `staging`). |
| `--execute` | off | Without it, dry run only. |
| `--yes` | off | Skip the confirmation prompt when submitting to `prod`. |

### Retrying after a partial failure

The script always POSTs the Document before the TabularFile. If the
Document succeeds but the TabularFile POST then fails (e.g. a schema
validation error), re-running the script as-is would try to create a
*second* Document. Instead, pass the `@id` the first run printed:

```bash
poetry run python -m src.submit_supplementary_data_1_to_igvf \
  --host prod --execute --existing-document-id /documents/<uuid-from-first-run>/
```

## After submission

The script prints the new Document and TabularFile `@id`s. Pass those along
with the two accessions/UUIDs being replaced to IGVF DACC curators so they
can revoke/archive the old objects and release the new ones.
