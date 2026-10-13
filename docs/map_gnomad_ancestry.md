# Map gnomAD Ancestry

`src/map_gnomad_ancestry.py` maps a variant set to gnomAD v4.1 exome
ancestry-specific allele frequencies (and genome frequencies for reference).
It is the second step of the Extended Data Figure 10 data-prep workflow,
run on `src/filter_pillar_data.py`'s output -- see
`docs/gnomad_ancestry_resolution.md` for the full walkthrough from
obtaining gnomAD data through generating the figure.

## What it does

For each variant:

1. Checks `ref_allele` against the GRCh38 reference (`match` /
   `revcomp_match` / `mismatch`).
2. Gives empty alleles (`-`, `''`, NaN) a VCF anchor base from the
   reference.
3. Left-aligns and trims against the reference (the `vt`/`bcftools norm`
   algorithm), so indels written at a different position within a repeat
   still match gnomAD's own representation.
4. Looks the normalized variant up in the gnomAD exome (and, optionally,
   genome) sites VCFs, normalizing each gnomAD record the same way; only an
   exact allele match counts.

The PASS filter is **not** applied here -- `exomes_FILTER` is kept as-is and
filtered later by `src/make_extended_data_figure_10.py`'s default
PASS-only behavior (`--include-filtered` to disable it).

## Inputs

- `--input` -- a variant CSV with `Chrom`, `hg38_start` (and optionally
  `hg38_end`), `ref_allele`, `alt_allele` columns -- typically one of
  `src/filter_pillar_data.py`'s `replicated_<category>_category_REVEL.csv`
  outputs.
- `--ref-fasta` -- an indexed GRCh38 FASTA (`.fai` alongside it; create with
  `samtools faidx`). Defaults to `data/input/gnomAD/reference/GRCh38.fa`.
- `--exomes-dir` -- a directory of gnomAD v4.1 exome sites VCFs, bgzipped
  and tabix-indexed, one per chromosome, named
  `gnomad.exomes.v4.1.sites.<chrom>.vcf.bgz` (gnomAD's own download naming
  convention, e.g. `gnomad.exomes.v4.1.sites.chr1.vcf.bgz`). Defaults to
  `data/input/gnomAD/vcf/exomes`.
- `--genomes-dir` -- the genome-sites equivalent
  (`gnomad.genomes.v4.1.sites.<chrom>.vcf.bgz`), optional; genome AF/AN/AC
  are recorded for reference but not used by the ancestry-group analysis. No
  default (conventionally `data/input/gnomAD/vcf/genomes` if provided).
- `--workers` (default 4) -- one lookup job per (source, chromosome) pair,
  run in a `multiprocessing.Pool`.

`--ref-fasta`/`--exomes-dir` are user-downloaded data, not checked into this
repo -- see `docs/gnomad_ancestry_resolution.md` for where to get
them and why they default to `data/input/gnomAD/` (gitignored).

## Output

`--output` gets the input columns plus:

- `gnomad_variant_type`, `ref_check`, `norm_pos`, `norm_ref`, `norm_alt`
- `exomes_found`, `exomes_AC`/`AN`/`AF`, `exomes_AC_<group>`/`AN_<group>`,
  `exomes_FILTER`, `exomes_record`
- `genomes_found`, `genomes_AC`/`AN`/`AF`, `genomes_FILTER`,
  `genomes_record` (if `--genomes-dir` given)
- `AF`, `AF_<group>` -- exome AFs (overall and per ancestry group), the
  columns `src/make_extended_data_figure_10.py` reads directly
- `exomes_MAF`, `genomes_MAF`, `gnomad_found_any`, `gnomad_sources`

Ancestry groups: `afr`, `amr`, `asj`, `eas`, `fin`, `mid`, `nfe`,
`remaining`, `sas` (`src/lib/gnomad.py`'s `GROUPS`, shared with
`src/make_extended_data_figure_10.py`).

Any row whose `gnomad_MAF` was already populated upstream but that didn't
match any gnomAD record is also written to
`<output>_unmatched_with_gnomad_MAF.csv`, for spot-checking normalization
misses.

## Usage

```bash
poetry run python -m src.map_gnomad_ancestry \
    --input data/intermediate/figures/extended_data_figure_10/replicated_VUS_category_REVEL.csv \
    --output data/intermediate/figures/extended_data_figure_10/vus_with_ancestry_af.csv
```

(assumes the reference FASTA and exome sites VCFs are already in place at
the `data/input/gnomAD/` defaults above; pass `--ref-fasta`/`--exomes-dir`/
`--genomes-dir` explicitly otherwise.)

Run once per category (`vus_with_ancestry_af.csv` / `gnomad_with_ancestry_af.csv`).
This writes to the gitignored `data/intermediate/figures/extended_data_figure_10/`
working directory, not directly to `data/input/gnomAD/` -- copy the result
there yourself when you want to update the checked-in snapshot
`src/make_extended_data_figure_10.py` reads by default; see
`docs/gnomad_ancestry_resolution.md`.

Requires `pysam` (an indexed-FASTA/VCF library wrapping `htslib`), already a
project dependency via `pyproject.toml`.
