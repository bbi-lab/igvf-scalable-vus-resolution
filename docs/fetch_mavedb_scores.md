# Fetch MaveDB Scores (via the API)

`src/fetch_mavedb_scores.py` is an API-based equivalent of
`src/mavedb_scores.sql`. Both produce `cvfg_variants.0.tsv`-shaped output -- one row
per variant, with score, mapped-HGVS, and target-gene columns -- but the SQL
version queries a local mirror of MaveDB's Postgres database, while this
script hits the public API (`https://api.mavedb.org`) directly. Use it when a
DB mirror isn't available.

## Output columns

Same header as `cvfg_variants.0.tsv`: `gene_symbol`, `score_set_title`,
`variant_urn`, `variant_id`, `raw_hgvs_nt`, `raw_hgvs_pro`, `score`,
`rna_score`, `rna_score_d6`, `rna_score_d20`, `score_data`, `reference_based`,
`mavedb_mapped_hgvs_g`, `mavedb_mapped_hgvs_c`, `mavedb_mapped_hgvs_p`,
`mavedb_mapping_error`, `target_sequence_type`, `target_sequence`,
`preferred_transcript`.

## Known differences from the SQL/DB version

Validated by diffing this script's output against a `src/mavedb_scores.sql`
export covering the same ~90 score sets: total row counts, per-score-set row
counts, `variant_urn` sets, and a field-by-field comparison of every shared
row. None of the differences found affect the variant-annotation pipeline's
results -- see the "doesn't matter because" note under each one.

- **No internal integer `variant_id`.** MaveDB's public API only exposes it
  via a per-variant GET (`/api/v1/variants/{urn}`), which isn't practical to
  call once per variant across the tens of thousands of variants in this
  dataset. The column is kept, for schema parity, but left blank.

  *Doesn't matter because* nothing downstream reads `variant_id`: it rides
  through the pipeline unused and is dropped by Step 21's final column
  selection before either output file is written.

- **`mavedb_mapped_hgvs_c`/`_p` are blank for many rows**, and
  `mavedb_mapped_hgvs_g` is occasionally expressed in a non-normalized form.
  The DB's `mapped_variants` table stores `hgvs_g`/`_c`/`_p` simultaneously
  (computed together by MaveDB's internal mapping pipeline), but the public
  API's mapped-variant model exposes only one HGVS expression per mapping --
  genomic for DNA-level score sets, protein for protein-level ones. In
  practice this leaves `_c`/`_p` blank in 32 of the 90 shared score sets (6 of
  those -- `urn:mavedb:00001224-a-1`, `00001226-a-1`/`-b-1`, `00001197-a-5`,
  `00001198-a-1`, `00001254-a-1` -- entirely blank; the other 26 only
  partially), and leaves `mavedb_mapped_hgvs_g` itself non-normalized
  (`delins` where the DB gives the canonical `del`/`dup` for the same change)
  on ~4,500 rows across ~15 score sets.

  *Doesn't matter because* `mavedb_mapped_hgvs_g`/`_c`/`_p` are produced but
  never read anywhere in the pipeline: they pass through Steps 1-20 untouched
  and are dropped at Step 21's final column selection. (The pipeline's actual
  HGVS mapping -- used for ClinGen allele lookups and the final output's own
  `hgvs_g`/`_c`/`_p` columns -- comes from a different, unprefixed
  `mapped_hgvs_g/c/p` column family, computed fresh in Step 1 by
  `vendor/variant-annotation/src/map_variants.py` from `raw_hgvs_nt`/
  `raw_hgvs_pro`.)

- **CALM1 (`urn:mavedb:00000001-c-1`): 497 of 2,800 rows pick a different
  duplicate-group representative.** This score set was uploaded at DNA level
  with the same protein-level score repeated across several synonymous
  nucleotide variants. The SQL dedupes with
  `distinct on (hgvs_pro) order by hgvs_pro, v.id`, keeping the row with the
  lowest internal DB id per group. The public API never exposes `v.id`, so
  this script substitutes the variant URN's numeric suffix as the tie-break
  key -- and that ordering doesn't track `v.id` at all (confirmed directly
  against the DB export's `variant_id` column: URN suffix and `v.id` are
  unrelated orderings).

  *Doesn't matter because* every member of a duplicate group shares the same
  `raw_hgvs_pro` and `score` by construction of the original upload, and the
  pipeline ignores `raw_hgvs_nt` here (the whole reason this score set needs
  deduplication). Both versions end up with the same 2,800 deduplicated
  protein-level variants and the same scores; only `variant_urn` and
  `mavedb_mapped_hgvs_g` (itself unused, see above) differ for the 497 rows
  that picked a different group member.

- **Floating-point noise in `score` and `score_data`.** Re-serializing
  through the API's JSON/CSV layer round-trips floats to within ~1e-16 of
  the DB's stored value (e.g. `0.0012318956870124` vs
  `0.00123189568701249`) -- last-bit rounding, not a different measurement.
  `score_data` also differs cosmetically in JSON key order; in whether a
  null-valued key is present as explicit `null` (DB) or omitted entirely
  (API, via `row.dropna()` before serializing); and, for a handful of score
  sets, in value type (e.g. integer `4` vs float `4.0`, or a number vs its
  string form). Checked directly: every non-null, non-formatting value
  matched.

  *Doesn't matter because* Step 11's score-range classification
  (`vendor/variant-annotation/src/annotate_mavedb.py`) does ordinary
  `>=`/`<=` bound comparisons against calibration thresholds, far coarser
  than 1e-16, and `score_data` itself is never parsed anywhere downstream --
  it rides through as an opaque column. The one place an *exact* float match
  would matter -- `src/build_variant_reclassification_dataset.py`'s CHEK2
  join against an Excel-sourced score column -- uses
  `urn:mavedb:00001205-a-1`, which has zero `score` differences between the
  two sources in this comparison.

- **Some score sets have no mapped-variants data in MaveDB**, and the API
  returns `404` for them rather than an empty list. The script treats this as
  "no mapped variants" and continues with that score set's scores; its
  `mavedb_mapped_hgvs_*`, `clingen_allele_id`, and `mavedb_mapping_error`
  columns are blank, same as the partial-coverage case above.

- **Three score sets this script couldn't fetch at all** (as of this
  writing) are simply absent from its output, rather than present with blank
  columns: `urn:mavedb:00001267-0-3`, `urn:mavedb:00001278-a-1`,
  `urn:mavedb:00001278-b-1`. Confirmed these are the *only* rows missing
  outright -- every other shared score set's `variant_urn` set matches the DB
  export exactly (aside from the CALM1 tie-break above).

- The live API returns transient `504`s fairly often; the script retries
  automatically (`--retries`, default 8, exponential backoff), but a full run
  across all ~90 score sets can still take a while and may need a retry.

## Source of truth: `src/mavedb_scores.sql`

The score-set URN list and the five dataset-specific manual overrides at the
end of `src/mavedb_scores.sql` (the CHEK2 preferred transcript, the BRCA1_Findlay
MANE Select version lift, the NDUFAF6/PTEN `raw_hgvs_nt` blanking, and the
JAG1/TARDBP(x2)/SCN5A/KCNH2 target-sequence corrections) are **parsed
directly out of that SQL file at run time** rather than duplicated here as
Python literals -- several of them are thousand-character DNA/protein
sequences, so re-typing them would risk a silent transcription error.
Parsing the SQL keeps this script from drifting out of sync with it if the
overrides are ever revised. Point `--sql-source` at a different file if you
need to override this.

## Usage

```bash
poetry run python -m src.fetch_mavedb_scores
```

This writes directly to `data/input/maves/cvfg_variants.0.tsv` by default --
see [Where `cvfg_variants.0.tsv` comes from](data.md#where-cvfg_variants0tsv-comes-from)
in `docs/data.md`. Pass `--output` to write elsewhere instead, e.g. when
comparing against the SQL/DB version's output.

Fetch a subset instead of the full ~90-dataset list (useful for testing
against the flaky live API), by repeating `--score-set-urn`:

```bash
poetry run python -m src.fetch_mavedb_scores \
  --output /tmp/smoke.tsv \
  --score-set-urn urn:mavedb:00000060-a-1 \
  --score-set-urn urn:mavedb:00001205-a-1
```

## CLI options

| Option | Default | Description |
|---|---|---|
| `--output` | `data/input/maves/cvfg_variants.0.tsv` | Where to write the resulting TSV |
| `--score-set-urn` | (full list from `src/mavedb_scores.sql`) | Fetch only this score-set URN; repeatable |
| `--sql-source` | `src/mavedb_scores.sql` | Source of the score-set list and manual overrides |
| `--retries` | `8` | HTTP retries per request on 5xx responses |
| `--timeout` | `120` | Per-request timeout, in seconds |
| `--sleep-between-requests` | `0.0` | Seconds to sleep between score sets, to go easier on the API |
