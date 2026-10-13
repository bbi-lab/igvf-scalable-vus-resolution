#!/usr/bin/env python
"""Map a variant set to gnomAD v4.1 exome ancestry-specific allele frequencies
(and genome frequencies for reference).

First step of the Extended Data Figure 10 data-prep workflow -- see
`docs/gnomad_ancestry_resolution.md` for the full
obtain-gnomAD-files-through-generate-the-figure walkthrough. Takes a variant
CSV (e.g. one of `src/filter_pillar_data.py`'s `replicated_<category>_category
_REVEL.csv` outputs) and adds the gnomAD columns `src/make_extended_data_
figure_10.py` expects (`Class_REVEL`, `AF_<group>`, `exomes_FILTER`).

For each variant:
  1. Check ref_allele against the GRCh38 reference (match / revcomp_match / mismatch).
  2. Give empty alleles ('-', '', NaN) a VCF anchor base from the reference.
  3. Left-align and trim against the reference (vt / bcftools-norm algorithm), so
     indels written at a different position in a repeat still match gnomAD.
  4. Look the normalized variant up in the gnomAD exome (and genome) sites VCFs,
     normalizing each gnomAD record the same way; only an exact allele match counts.

Output = input columns plus:
  gnomad_variant_type, ref_check, norm_pos, norm_ref, norm_alt
  exomes_found, exomes_AC/AN/AF, exomes_AC_<group>, exomes_AN_<group>, exomes_FILTER, exomes_record
  genomes_found, genomes_AC/AN/AF, genomes_FILTER, genomes_record
  AF, AF_<group>       exome AFs (overall and per ancestry group) used by the analysis
  exomes_MAF, genomes_MAF, gnomad_found_any, gnomad_sources
Groups: afr, amr, asj, eas, fin, mid, nfe, remaining, sas.
The PASS filter is applied later, by `src/make_extended_data_figure_10.py`'s
`--include-filtered`/default PASS-only logic (exomes_FILTER is kept here).

Requires an indexed GRCh38 FASTA (.fai alongside it) and local gnomAD v4.1
sites VCFs (bgzipped + tabix-indexed, one per chromosome, named
`gnomad.exomes.v4.1.sites.<chrom>.vcf.bgz` / `gnomad.genomes.v4.1.sites.<chrom>.vcf.bgz`
-- gnomAD's own download naming convention). These are large, user-downloaded
files that don't belong in the repo -- by convention they live under
`data/input/gnomAD/` (gitignored, like `data/input/biobank/`), which is also
where this script's own output is meant to land. See
`docs/gnomad_ancestry_resolution.md` for where to get these.

Usage:
    poetry run python -m src.map_gnomad_ancestry \\
        --input variants.csv \\
        --output data/input/gnomAD/variants_with_ancestry_af.csv
    # --ref-fasta/--exomes-dir/--genomes-dir default to data/input/gnomAD/reference
    # and data/input/gnomAD/vcf/{exomes,genomes} -- override with explicit paths
    # if your copies live elsewhere.
"""

import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import click
import numpy as np
import pandas as pd

from src.lib.gnomad import GROUPS

DEFAULT_REF_FASTA = Path("data/input/gnomAD/reference/GRCh38.fa")
DEFAULT_EXOMES_DIR = Path("data/input/gnomAD/vcf/exomes")
DEFAULT_GENOMES_DIR = Path("data/input/gnomAD/vcf/genomes")

VCF_NAME = {
    "exomes": "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz",
    "genomes": "gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz",
}
FIELDS = {  # gnomAD INFO fields pulled from each source
    "exomes": ["AC", "AN", "AF"] + [f"{x}_{g}" for x in ("AC", "AN", "AF") for g in GROUPS],
    "genomes": ["AC", "AN", "AF"],
}
EMPTY = {"", "-", ".", "nan", "none", "na", "<na>"}
COMP = str.maketrans("ACGTN", "TGCAN")


def log(msg=""):
    print(msg, file=sys.stderr, flush=True)


# ---------------- allele normalization (seq(chrom, start, end): 1-based inclusive) ----------------
def clean_allele(a):
    if a is None or (isinstance(a, float) and np.isnan(a)):
        return ""
    a = str(a).strip().upper()
    return "" if a.lower() in EMPTY else a


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse_int(x):
    try:
        return int(float(str(x).strip()))
    except (TypeError, ValueError):
        return None


def ref_check(chrom, start, ref, seq):
    if not ref:
        return "no_ref_allele"
    genome = seq(chrom, start, start + len(ref) - 1)
    if genome == ref:
        return "match"
    if genome == revcomp(ref):
        return "revcomp_match"
    return "mismatch"


def to_vcf(chrom, start, end, ref, alt, seq):
    """Anchor empty alleles so both are non-empty (VCF style)."""
    if ref and alt:
        return start, ref, alt
    if ref and not alt:  # deletion: start = first deleted base
        a = seq(chrom, start - 1, start - 1)
        return start - 1, a + ref, a
    if alt and not ref:  # insertion: anchor = left flanking base
        anchor = min(start, end) if end is not None else start - 1
        a = seq(chrom, anchor, anchor)
        return anchor, a, a + alt
    return None


def left_align(chrom, pos, ref, alt, seq):
    """Trim a matching right base, extend left when an allele empties, repeat;
    then trim matching left bases while both alleles are >= 2 long."""
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1] and not (len(ref) == 1 and len(alt) == 1):
            ref, alt, changed = ref[:-1], alt[:-1], True
        if (not ref or not alt) and pos > 1:
            b = seq(chrom, pos - 1, pos - 1)
            ref, alt, pos, changed = b + ref, b + alt, pos - 1, True
        if not changed:
            break
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"  # gnomAD stores these as separate SNVs, so they can't match one record
    if len(alt) == 1 and ref[0] == alt[0]:
        return "deletion"
    if len(ref) == 1 and alt[0] == ref[0]:
        return "insertion"
    return "complex_indel"


def normalize_row(chrom, start, end, ref, alt, seq):
    """-> (key or None, variant_type, ref_check)."""
    if start is None:
        return None, "unparseable", "n/a"
    ref, alt = clean_allele(ref), clean_allele(alt)
    if not ref and not alt or ref == alt:
        return None, "invalid_alleles", "n/a"
    rc = ref_check(chrom, start, ref, seq)
    v = to_vcf(chrom, start, end, ref, alt, seq)
    if v is None:
        return None, "invalid_alleles", rc
    pos, r, a = left_align(chrom, *v, seq)
    return (chrom, pos, r, a), variant_type(r, a), rc


def match_record(candidates, key, seq):
    """Return (record, alt_index) whose normalized allele equals `key`."""
    chrom, pos, ref, alt = key
    for rec in candidates:
        for i, ra in enumerate(rec.alts or []):
            if ra is None or ra.startswith("<") or ra == "*":
                continue
            if left_align(chrom, rec.pos, rec.ref.upper(), ra.upper(), seq) == (pos, ref, alt):
                return rec, i
    return None, None


# ---------------- gnomAD lookup ----------------
def make_seq(fasta):
    names = set(fasta.references)

    def contig(chrom):
        if chrom in names:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in names:
            return alt
        if chrom in ("chrM", "chrMT", "MT", "M"):
            for c in ("chrM", "MT", "M"):
                if c in names:
                    return c
        raise KeyError(f"{chrom} not in reference FASTA")

    def seq(chrom, start, end):
        return fasta.fetch(contig(chrom), start - 1, end).upper()

    return seq


def extract(rec, alt_i, fields):
    out = {}
    for f in fields:
        try:
            v = rec.info.get(f)
        except (KeyError, ValueError):
            v = None
        if isinstance(v, tuple):
            v = v[alt_i] if len(v) > alt_i else (v[0] if len(v) == 1 else None)
        if v is not None:
            out[f] = v
    keys = list(rec.filter.keys())
    out["FILTER"] = ";".join(keys) if keys else "PASS"
    out["record"] = f"{rec.contig}:{rec.pos} {rec.ref}>{rec.alts[alt_i]}"
    return out


def worker(job):
    """Look up one chromosome's variants in one gnomAD source."""
    import pysam

    source, path, chrom, keys, fasta_path = job
    t0 = time.time()
    try:
        vcf = pysam.VariantFile(path)
    except (OSError, ValueError) as e:
        log(f"WARNING [{source}] cannot open {path}: {e}")
        return source, {}
    seq = make_seq(pysam.FastaFile(fasta_path))
    in_header = {r.get("ID") for r in vcf.header.records if r.type == "INFO"}
    fields = [f for f in FIELDS[source] if f in in_header]
    contig = chrom if chrom in vcf.header.contigs else chrom.replace("chr", "")
    res = {}
    for key in keys:
        _, pos, ref, _alt = key
        try:
            cands = vcf.fetch(contig, max(pos - 6, 0), pos + len(ref) + 5)
        except ValueError:
            continue
        rec, i = match_record(cands, key, seq)
        if rec is not None:
            res[key] = extract(rec, i, fields)
    log(f"[{source}] {chrom}: {len(res)}/{len(keys)} matched ({time.time() - t0:.0f}s)")
    return source, res


# ---------------- main ----------------
def map_gnomad_ancestry(input_path, ref_fasta, exomes_dir, genomes_dir, workers):
    """Run the full mapping and return the annotated dataframe. Writes
    nothing; `main()` below handles output. Raises ValueError for bad
    input (missing required columns, unreadable CSV)."""
    import pysam

    try:
        df = pd.read_csv(input_path, dtype=str, low_memory=False)
    except pd.errors.ParserError as e:
        raise ValueError(
            f"Could not parse {input_path}: {e}\n'EOF inside string' usually means the file is "
            f"truncated (e.g. an interrupted copy) -- compare its size/md5 with the original."
        ) from e
    missing = [c for c in ["Chrom", "hg38_start", "ref_allele", "alt_allele"] if c not in df.columns]
    if missing:
        raise ValueError(f"{input_path} is missing required column(s): {', '.join(missing)}")
    log(f"Loaded {len(df)} rows from {input_path}")

    sources = {"exomes": exomes_dir}
    if genomes_dir:
        sources["genomes"] = genomes_dir

    # --- normalize every variant ---
    seq = make_seq(pysam.FastaFile(str(ref_fasta)))
    chroms = "chr" + df["Chrom"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True).str.replace(
        "^chr", "", regex=True
    )
    starts = df["hg38_start"].map(parse_int)
    ends = df["hg38_end"].map(parse_int) if "hg38_end" in df.columns else pd.Series([None] * len(df))
    keys, vtypes, rchecks = [], [], []
    for c, s, e, r, a in zip(chroms, starts, ends, df["ref_allele"], df["alt_allele"]):
        try:
            k, vt, rc = normalize_row(c, s, e, r, a, seq)
        except (KeyError, ValueError) as err:
            k, vt, rc = None, "reference_lookup_failed", str(err)[:40]
        keys.append(k)
        vtypes.append(vt)
        rchecks.append(rc)
    df["gnomad_variant_type"], df["ref_check"] = vtypes, rchecks
    df["norm_pos"] = [k[1] if k else None for k in keys]
    df["norm_ref"] = [k[2] if k else None for k in keys]
    df["norm_alt"] = [k[3] if k else None for k in keys]

    log("\nVariant types after normalization:\n" + df["gnomad_variant_type"].value_counts().to_string())

    # --- look up each unique variant, one job per source x chromosome ---
    by_chrom = {}
    for k in sorted({k for k in keys if k is not None}):
        by_chrom.setdefault(k[0], []).append(k)
    jobs = [
        (s, os.path.join(d, VCF_NAME[s].format(chrom=c)), c, ks, str(ref_fasta))
        for s, d in sources.items()
        for c, ks in by_chrom.items()
    ]
    log(
        f"\nMapping {sum(map(len, by_chrom.values()))} unique variants against {list(sources)} "
        f"({len(jobs)} jobs, {workers} workers)"
    )
    found = {s: {} for s in sources}
    with Pool(workers) as pool:
        for s, res in pool.imap_unordered(worker, jobs):
            found[s].update(res)

    # --- output columns ---
    for s in sources:
        recs = [found[s].get(k) if k else None for k in keys]
        df[f"{s}_found"] = [r is not None for r in recs]
        for f in [f for f in FIELDS[s] if not (s == "exomes" and f.startswith("AF_"))] + ["FILTER", "record"]:
            df[f"{s}_{f}"] = [r.get(f) if r else None for r in recs]
        if s == "exomes":  # exome AFs used by the analysis, unprefixed
            for f in ["AF"] + [f"AF_{g}" for g in GROUPS]:
                df[f] = pd.to_numeric([r.get(f) if r else None for r in recs], errors="coerce")
        af = pd.to_numeric(df[f"{s}_AF"], errors="coerce")
        df[f"{s}_MAF"] = np.minimum(af, 1 - af)
    found_cols = [f"{s}_found" for s in sources]
    df["gnomad_found_any"] = df[found_cols].any(axis=1)
    df["gnomad_sources"] = df[found_cols].apply(
        lambda r: "+".join(s for s, v in zip(sources, r) if v) or "not_found", axis=1
    )
    return df


def print_qc(df, sources, output_path):
    print(
        f"\n=== {len(df)} rows; found in gnomAD (any source): {df['gnomad_found_any'].sum()} "
        f"({df['gnomad_found_any'].mean() * 100:.1f}%) ==="
    )
    print(df["gnomad_sources"].value_counts().to_string())
    for s in sources:
        print(
            f"[{s}] matched records failing FILTER (non-PASS): "
            f"{(df[f'{s}_found'] & (df[f'{s}_FILTER'] != 'PASS')).sum()}"
        )
    print("\nMatch rate by variant type:")
    print(
        df.groupby("gnomad_variant_type")["gnomad_found_any"]
        .agg(["size", "sum", "mean"])
        .rename(columns={"size": "n", "sum": "found", "mean": "rate"})
        .to_string()
    )
    if "gnomad_MAF" in df.columns:
        has_maf = df["gnomad_MAF"].notna()
        print("\nAgreement with the dataset's own gnomad_MAF column:")
        print(
            pd.crosstab(
                has_maf.map({True: "gnomad_MAF present", False: "gnomad_MAF absent"}),
                df["gnomad_found_any"].map({True: "found", False: "not found"}),
            ).to_string()
        )
        miss = df[has_maf & ~df["gnomad_found_any"]]
        if len(miss):
            print(f"\n{len(miss)} rows have gnomad_MAF but no match -- by type / ref_check:")
            print(pd.crosstab(miss["gnomad_variant_type"], miss["ref_check"]).to_string())
            unmatched_path = output_path.with_name(output_path.stem + "_unmatched_with_gnomad_MAF.csv")
            miss.to_csv(unmatched_path, index=False)


@click.command(help=__doc__)
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="CSV of variants (Chrom, hg38_start, hg38_end, ref_allele, alt_allele)",
)
@click.option("--output", "output_path", required=True, type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--ref-fasta",
    default=DEFAULT_REF_FASTA,
    show_default=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="indexed GRCh38 FASTA",
)
@click.option(
    "--exomes-dir",
    default=DEFAULT_EXOMES_DIR,
    show_default=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="directory of gnomAD v4.1 exome sites VCFs",
)
@click.option(
    "--genomes-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help=f"directory of gnomAD v4.1 genome sites VCFs (optional; conventionally {DEFAULT_GENOMES_DIR})",
)
@click.option("--workers", type=int, default=4, show_default=True)
def main(input_path, output_path, ref_fasta, exomes_dir, genomes_dir, workers):
    if not Path(str(ref_fasta) + ".fai").exists():
        raise click.ClickException(f"No index {ref_fasta}.fai -- create it with: samtools faidx {ref_fasta}")

    try:
        df = map_gnomad_ancestry(input_path, ref_fasta, exomes_dir, genomes_dir, workers)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    sources = ["exomes"] + (["genomes"] if genomes_dir else [])
    print_qc(df, sources, output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    grp = [f"AF_{g}" for g in GROUPS]
    print(
        f"\nRows with exome ancestry AFs: {df[grp].notna().any(axis=1).sum()}; nonzero AF per group: "
        + ", ".join(f"{g}={(df[f'AF_{g}'] > 0).sum()}" for g in GROUPS)
    )
    click.echo(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
