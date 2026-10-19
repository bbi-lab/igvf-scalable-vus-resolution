# Submit Data Prerequisites to Zenodo

`src/submit_data_prerequisites_to_zenodo.py` creates a new version of the
existing *"A scalable approach to resolving variants of uncertain
significance (supplementary data and supporting files)"* Zenodo record
(18637474), retitled *"Data prerequisites for ..."*, and uploads the
pipeline's not-checked-into-GitHub input data plus the project's own
figure/table source files. See
[`docs/data_prerequisites_zenodo_manifest.md`](data_prerequisites_zenodo_manifest.md)
for the full per-file license/citation rationale behind what is and isn't
included -- that document is itself uploaded into the deposit as
`LICENSE_AND_CITATIONS.md`, so downstream users have it without needing this
repo.

## Before running

Fetch the public MutPred2 training-variant source (not done automatically --
dry runs make no network calls, on principle):

```bash
curl -s https://mutpred.mutdb.org/wo_exclusive_hgmd_mp2_training_data.txt \
  | gzip > data/input/predictors/mp2_training_data.public.txt.gz
```

Every other file the manifest needs should already be in place per
[`docs/data.md`](data.md).

## Usage

Dry run (default -- no network calls; lists every file the manifest would
upload with its size, and flags the two non-commercial-restricted groups):

```bash
poetry run python -m src.submit_data_prerequisites_to_zenodo
```

Real run against the sandbox first. Sandbox is a separate database from
production with no copy of record 18637474, so `--new` creates a fresh test
deposition instead of versioning it -- its metadata (creators, affiliations,
etc.) is seeded from production's own public listing via an unauthenticated
GET, not fabricated:

```bash
export ZENODO_TOKEN=...   # a sandbox.zenodo.org personal access token
poetry run python -m src.submit_data_prerequisites_to_zenodo \
  --host sandbox --new --execute
```

Check the result in the sandbox web UI (the script prints the draft's URL).
Once it looks right, create a new version of the real production record
(asks for interactive confirmation unless `--yes` is also passed):

```bash
export ZENODO_TOKEN=...   # a zenodo.org personal access token -- different from sandbox's
poetry run python -m src.submit_data_prerequisites_to_zenodo --host prod --execute
```

This only creates/updates a **draft** -- review it at the printed URL.
Publishing a Zenodo version is irreversible (the DOI is permanent; a
published file can't be edited or removed, only superseded by a further new
version), so it needs a second, explicit flag once you're satisfied:

```bash
poetry run python -m src.submit_data_prerequisites_to_zenodo --host prod --execute --publish
```

### Options

| Option | Default | Purpose |
|---|---|---|
| `--repo-root` | this checkout | Where to read `data/`/`docs/` from. |
| `--mp2-public-source` | `data/input/predictors/mp2_training_data.public.txt.gz` | Local cache of the public MutPred2 training source (see "Before running"). |
| `--host` | `sandbox` | `sandbox` or `prod`. |
| `--deposition-id` | `18637473` (the record's concept recid) | Record to create a new version of, resolved to its current latest published version automatically -- see "Why a concept recid, not a version id" below. Ignored with `--new`. |
| `--new` | off | Create a brand-new deposition instead of versioning `--deposition-id`. Required on `--host sandbox`. |
| `--execute` | off | Without it, dry run only. |
| `--publish` | off | Publish the draft after updating it. Requires `--execute`. Irreversible. |
| `--yes` | off | Skip the confirmation prompt on `--host prod`. |

### Why a concept recid, not a version id

Zenodo's `newversion` action requires the specific latest version's own
record id -- it will not accept the concept id directly, and that id
changes every time a new version is published. Hardcoding a version id
would need updating by hand after every release. Instead, `--deposition-id`
defaults to the record's *concept* recid (`18637473`, stable forever --
it's the numeric suffix of the concept DOI, `10.5281/zenodo.18637473`), and
the script resolves it to the actual latest version's id itself via an
unauthenticated `GET /records/<id>` against the public records API (Zenodo
resolves any id in a version chain -- concept or a specific version -- to
the record representing the current latest published version). You can
still pass a specific version's id instead if you ever need to, since the
same resolution call is a no-op for an id that's already the latest.

## What gets uploaded

See `build_file_manifest()` in the script and
[`docs/data_prerequisites_zenodo_manifest.md`](data_prerequisites_zenodo_manifest.md)
for the authoritative, up-to-date list. In brief: this project's own curated
inputs and figure/table sources, permissively-licensed third-party reference
data, and two groups of files flagged non-commercial-restricted (REVEL
predictor scores; the gnomAD-ancestry files, which also embed REVEL and
SpliceAI scores). The 87 per-dataset ExCALIBR calibration JSON files are
zipped into a single `mave_calibration__excalibr__json.zip` rather than
uploaded individually. Keys use `__` rather than `/` as the path separator
(Zenodo's legacy upload API rejects a literal `/` in a filename) -- see
"Where each file goes" in `docs/data_prerequisites_zenodo_manifest.md`.

## What this script does not decide for you

It does not choose the file list -- that's `build_file_manifest()`, which
encodes the decisions in `docs/data_prerequisites_zenodo_manifest.md`.
Update both together if the file list changes. It also does not invent
deposit metadata (creators, affiliations): new-version mode inherits that
automatically from the record being versioned; `--new` mode reuses
production's own public metadata instead of fabricating it.

## Retrying after a partial failure

Re-running the same command is safe: `--execute` always re-derives the
manifest and re-uploads every file (overwriting same-keyed files already in
the draft), and metadata is a plain overwrite rather than an append (aside
from `related_identifiers`, which accumulates -- harmless duplicate entries
if you re-run, but worth checking the draft before `--publish`). If a run
fails partway through the upload loop, just re-run it; already-uploaded
files are simply replaced.

### A known Zenodo quirk: "Please remove all files first"

If an earlier run left an unpublished new-version draft sitting on the
record (a previous failure partway through, or even just a `newversion`
call that timed out on Zenodo's end after actually creating the draft),
Zenodo's legacy `actions/newversion` endpoint 400s on every subsequent
attempt with:

```
{"field": "files.enabled", "messages": ["Please remove all files first."]}
```

This is a known Zenodo API issue, not specific to this project --
`create_new_version` works around it by using the newer
`POST /records/<id>/versions` endpoint instead, which resumes the existing
draft rather than erroring. You shouldn't need to do anything about this
yourself; it's mentioned here so the error (if you ever see it from some
other tool, or from a version of this script predating the fix) doesn't
look like a sign your draft is broken -- it isn't, and nothing needs to be
deleted or recreated by hand.
