# Submit Result Data Files to Zenodo

`src/submit_results_to_zenodo.py` uploads the paper's result data files as a
new version of the existing *"Data from 'A scalable approach to resolving
variants of uncertain significance'"* Zenodo record
([zenodo.org/records/23308564](https://zenodo.org/records/23308564) at the
time of writing). Earlier attempts to publish Supplementary Data 5 and 6 to
this record didn't go through; this script replaces whatever ad hoc process
was used (there's no prior script or doc for it in this repo -- it appears
the files were uploaded by hand through Zenodo's web UI) with the same
API-based approach already proven out by
[`submit_data_prerequisites_to_zenodo.py`](submit_data_prerequisites_to_zenodo.md)
for the sibling "data prerequisites" record, sharing its client code via
`src/lib/zenodo_client.py`.

Unlike that script, this one does not touch the record's metadata (title,
description, license, related identifiers) -- it's already correct from
when the record was first created. This script's only job is getting the
files uploaded.

## What gets uploaded

| File | Size (approx.) |
|---|---|
| `Supplementary_Data_1.tsv.gz` (the integrated variant effect dataset) | 52 MB |
| `Supplementary_Data_1_schema.xlsx` (its column documentation) | 19 KB |
| `Supplementary_Data_3.xlsx` | 660 KB |
| `Supplementary_Data_5.xlsx` | 211 MB |
| `Supplementary_Data_5.with_secondary_variants.xlsx` | 225 MB |
| `Supplementary_Data_6.xlsx` | 181 MB |
| `Supplementary_Data_6.with_secondary_variants.xlsx` | 192 MB |

~861 MB total. All read from `data/output/supplementary_data/` --
regenerate them per [`docs/figures.md`](figures.md) /
[`README.md`](../README.md)'s Section 2 if they're missing or stale before
running this script.

**Open question, worth confirming before a real run:** Supplementary Data 5
and 6 each also have a `.xlsx.gz` form on disk (gzip-compressed, ~20%
smaller, identical content). This manifest uses the plain `.xlsx` files, for
consistency with `Supplementary_Data_3.xlsx` (not gzipped) and because
that's the form someone would open directly in Excel. Switch `_STATIC_FILES`
in the script to the `.xlsx.gz` paths instead if you'd rather minimize
upload/download size.

## Usage

Dry run (default -- no network calls; lists every file the manifest would
upload with its size):

```bash
poetry run python -m src.submit_results_to_zenodo
```

Real run against the sandbox first. Sandbox is a separate database from
production with no copy of this record, so `--new` creates a fresh test
deposition instead of versioning it -- its metadata is seeded from
production's own public listing via an unauthenticated GET, not fabricated:

```bash
export ZENODO_TOKEN=...   # a sandbox.zenodo.org personal access token
poetry run python -m src.submit_results_to_zenodo --host sandbox --new --execute
```

Check the result in the sandbox web UI (the script prints the draft's URL).
Once it looks right, create a new version of the real production record
(asks for interactive confirmation unless `--yes` is also passed):

```bash
export ZENODO_TOKEN=...   # a zenodo.org personal access token -- different from sandbox's
poetry run python -m src.submit_results_to_zenodo --host prod --execute
```

This only creates/updates a **draft** -- review it at the printed URL.
Publishing a Zenodo version is irreversible (the DOI is permanent; a
published file can't be edited or removed, only superseded by a further new
version), so it needs a second, explicit flag once you're satisfied:

```bash
poetry run python -m src.submit_results_to_zenodo --host prod --execute --publish
```

### Options

| Option | Default | Purpose |
|---|---|---|
| `--repo-root` | this checkout | Where to read `data/output/supplementary_data/` from. |
| `--host` | `sandbox` | `sandbox` or `prod`. |
| `--deposition-id` | `23308563` (concept recid) | Record to create a new version of, resolved to its current latest published version automatically. Ignored with `--new`. |
| `--new` | off | Create a brand-new deposition instead of versioning `--deposition-id`. Required on `--host sandbox`. |
| `--execute` | off | Without it, dry run only. |
| `--publish` | off | Publish the draft after updating it. Requires `--execute`. Irreversible. |
| `--yes` | off | Skip the confirmation prompt on `--host prod`. |

Large files (up to 225 MB here) over a flaky connection can still fail
mid-upload even via the API -- if a run fails partway through, just re-run
it; `create_new_version` resumes the existing draft rather than erroring
(see `src/lib/zenodo_client.py`), and already-uploaded files are simply
replaced on re-upload.
