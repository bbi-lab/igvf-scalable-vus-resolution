#!/usr/bin/env python
"""Build the VUS / gnomAD / Unobserved category datasets used by Extended Data
Figure 10, replicating `notebooks/analysis/Variant_Classification_analysis.ipynb`'s
steps in the same order:

First step of the Extended Data Figure 10 data-prep workflow -- see
`docs/gnomad_ancestry_resolution.md` for the full
obtain-gnomAD-files-through-generate-the-figure walkthrough. Its output feeds
`src/map_gnomad_ancestry.py` as `--input`.

   1. Functional points for every row: excalibr_points, with the Supp 4 ExCALIBR
      intervals filling rows ExCALIBR left unscored (flag / vep_splice /
      spliceai); OddsPath points for F9/TP53; missing points -> 0
   2. Keep only the meta-analysis dataset for TP53 and F9
   3. Tag (don't yet remove) splice variants their assay doesn't measure, and
      start-lost variants in amino-acid-level assays
   4. Tag variants whose functional points conflict in sign between assays,
      within nucleotide-level and within amino-acid-level assays
   5. Spread splice and REVEL-training flags to every nucleotide variant
      encoding the same amino-acid change
   6. Keep only nucleotide- and amino-acid-level rows
   7. Tag sign conflicts between nucleotide- and amino-acid-level assays, and
      extend each tag to every row sharing that measurement (mavedb_variant_urn)
   8. Remove SFPQ
   9. Apply the CHEK2 QC flag (redundant: already set upstream)
  10. Remove tagged conflicts, unmeasured splice variants and start-lost variants
  11. Remove rows with Flag == '*'
  12. Select the category (VUS / gnomAD / Unobserved), keep one row per variant
      (greatest absolute functional points, Dataset name as tie-break), and
      drop REVEL training variants
  13. Classify each variant: functional points + REVEL points (gene-specific
      threshold, genome-wide fallback, else 0) -> Class_REVEL (P >=10, LP 6-9,
      Uncertain 0-5, LB -6..-1, B <=-7); flag functional-vs-REVEL disagreement
      in Conflicting_REVEL

Usage:
  poetry run python -m src.filter_pillar_data \\
      --points-input data/input/maves/dataframe_with_points.csv.gz \\
      --supp4-input data/output/supplementary_data/Supplementary_Data_4.xlsx \\
      --supp3-input data/input/maves/Supplementary_Data_3.xlsx \\
      --category VUS
Ablations: --no-excalibr-backfill, --no-urn-propagation
"""

import re
import sys
import unicodedata
from pathlib import Path

import click
import numpy as np
import pandas as pd

DEFAULT_POINTS_INPUT = Path("data/input/maves/dataframe_with_points.csv.gz")
DEFAULT_OUTDIR = Path("data/intermediate/figures/extended_data_figure_10")

TP53_F9_EXCLUDED = [
    "TP53_Boettcher_2019",
    "TP53_Fortuno_2021",
    "TP53_Giacomelli_2018_combined_score",
    "TP53_Giacomelli_2018_p53WT_Nutlin3",
    "TP53_Giacomelli_2018_p53null_Nutlin3",
    "TP53_Giacomelli_2018_p53null_etoposide",
    "TP53_Kato_2003_AIP1nWT",
    "TP53_Kato_2003_BAXnWT",
    "TP53_Kato_2003_GADD45nWT",
    "TP53_Kato_2003_MDM2nWT",
    "TP53_Kato_2003_NOXAnWT",
    "TP53_Kato_2003_P53R2nWT",
    "TP53_Kato_2003_WAF1nWT",
    "TP53_Kato_2003_h1433snWT",
    "F9_Popp_2025_carboxy_F9_specific",
    "F9_Popp_2025_carboxy_gla_motif",
    "F9_Popp_2025_heavy_chain",
    "F9_Popp_2025_light_chain",
    "F9_Popp_2025_strep_2",
]
OP_GENES = ["F9", "TP53"]
PRIORITY_GENES_2018 = ["BRCA1", "PTEN", "MSH2"]  # genes using 2018-vintage ExCALIBR calibrations
UNSCORED_REASONS = ["flag", "vep_splice", "spliceai"]
DIS = [  # tags removed before category selection
    "conflicting_fxn_data",
    "splice_variant_not_measured",
    "splice_variant_not_measured;conflicting_fxn_data",
    "start_lost_variant_not_measured",
]
DISALLOWED = [  # tags the cross-resolution pass leaves as is
    "splice_variant_not_measured",
    "splice_variant_not_measured;conflicting_fxn_data",
    "conflicting_fxn_data",
]
FUNC_MAP = {"normal": "Normal", "abnormal": "Abnormal", "not_specified": "Indeterminate"}


def log(msg=""):
    print(msg, file=sys.stderr)


def nfc(s):
    return unicodedata.normalize("NFC", s) if isinstance(s, str) else s


def group_any(df, keys, values, required=()):
    """Per-group any(values), broadcast to rows, returned as a plain bool Series.
    Rows with NaN in a `required` column are left out of grouping and get False --
    mimics the notebook's default groupby(dropna=True). Groups on a single
    concatenated string key so behavior doesn't depend on the pandas version's
    groupby(dropna=...) handling."""
    out = np.zeros(len(df), dtype=bool)
    ok = np.ones(len(df), dtype=bool)
    for c in required:
        ok &= df[c].notna().to_numpy()
    if ok.any():
        sub = df.loc[ok, keys].astype(str)
        key = sub[keys[0]]
        for c in keys[1:]:
            key = key + "\x1f" + sub[c]
        v = pd.Series(np.asarray(values, dtype=bool)[ok].astype(np.int8), index=key.index)
        out[ok] = v.groupby(key.values).transform("max").to_numpy() > 0
    return pd.Series(out, index=df.index, dtype=bool)


def opposite_signs(df, keys, pts, required=()):
    """has_opposite_signs: group has a strictly positive AND strictly negative value."""
    return (group_any(df, keys, pts > 0, required) & group_any(df, keys, pts < 0, required)).astype(bool)


# ---------------- ExCALIBR functional points ----------------
def load_excalibr(path):
    ex = pd.read_excel(path, sheet_name="ExCALIBR_calibrations")
    ex["dataset"] = ex["dataset"].apply(nfc)
    is18 = ex["dataset"].str.endswith("_clinvar_2018")
    main_, e18 = ex[~is18].copy(), ex[is18].copy()
    main_["base"] = main_["dataset"]
    e18["base"] = e18["dataset"].str.replace("_clinvar_2018", "", regex=False)
    m = main_.add_suffix("_2025").merge(e18.add_suffix("_2018"), left_on="base_2025", right_on="base_2018", how="outer")
    m["base"] = m["base_2025"].fillna(m["base_2018"])
    m = m.drop(columns=["base_2025", "base_2018"])
    if m["base"].duplicated().any():
        log("WARNING: duplicate dataset names in ExCALIBR sheet -- keeping first")
        m = m.drop_duplicates("base")
    return m


def points_from_intervals(row, cols):
    score = row["_score"]
    for col in cols:
        cell = row[col]
        if not isinstance(cell, str):
            continue
        for chunk in cell.split(","):
            lo, hi = chunk.strip().split()
            lower = -np.inf if lo.lower() == "-inf" else float(lo)
            upper = np.inf if hi.lower() == "inf" else float(hi)
            if ((lower == -np.inf) or (score >= lower)) and ((upper == np.inf) or (score < upper)):
                return int(re.search(r"^range_(-?\d+)", col).group(1))
    return np.nan


def exc_fxn(sub, ex):
    """Fxn_points from ExCALIBR intervals, using 2018-vintage calibrations for BRCA1/PTEN/MSH2."""
    if len(sub) == 0:
        return pd.Series(dtype=float)
    m = sub[["Dataset", "Gene", "auth_reported_score"]].copy()
    m["_rowid"] = m.index
    m["_score"] = pd.to_numeric(m["auth_reported_score"].astype(str).str.strip(), errors="coerce")
    m = m.merge(ex, left_on="Dataset", right_on="base", how="left").set_index("_rowid")
    r25 = [c for c in m.columns if c.startswith("range_") and c.endswith("_2025")]
    r18 = [c for c in m.columns if c.startswith("range_") and c.endswith("_2018")]
    p25 = m.apply(lambda r: points_from_intervals(r, r25), axis=1).astype(float)
    p18 = m.apply(lambda r: points_from_intervals(r, r18), axis=1).astype(float)
    fxn = p25.copy()
    pri = m["Gene"].isin(PRIORITY_GENES_2018)
    fxn[pri] = p18[pri].fillna(p25[pri])
    return fxn.reindex(sub.index)


# ---------------- OddsPath functional points for F9/TP53 ----------------
def parse_bound(s):
    s = s.strip().lower()
    if s in ["-inf", "-infinity"]:
        return float("-inf")
    if s in ["inf", "infinity", "+inf"]:
        return float("inf")
    return float(s)


def find_matching_class(row, score):
    for i in range(1, 10):
        rc, cc = f"Interval {i} Range", f"Interval {i} Class"
        if rc in row.index and isinstance(row[rc], str):
            iv = row[rc]
            try:
                lo_s, hi_s = iv[1:-1].split(",")
                lo, hi = parse_bound(lo_s), parse_bound(hi_s)
                ok_lo = score >= lo if iv[0] == "[" else score > lo
                ok_hi = score <= hi if iv[-1] == "]" else score < hi
                if ok_lo and ok_hi:
                    return row[cc]
            except (ValueError, IndexError):
                continue
    return "Not specified"


def op_points(sub, supp4, supp3):
    op = pd.read_excel(supp4, sheet_name="OddsPath_calibrations").drop_duplicates("Dataset")
    m = sub.copy()
    m["_rowid"] = m.index
    m = m.merge(op[["Dataset", "OddsNormal", "OddsAbnormal"]], on="Dataset", how="left")
    cat = m["auth_reported_func_class_category"]
    cls = cat.astype(str).str.strip().str.lower().map(FUNC_MAP)
    cls[cat.isna()] = np.nan
    m["StandardizedClass"] = cls
    if supp3:
        mave = pd.read_excel(supp3, sheet_name="Curation", header=0)[
            ["Dataset Name", "Score Intervals Reported?", "Functional Classification Provided?"]
        ].drop_duplicates("Dataset Name")
        m = m.merge(mave, left_on="Dataset", right_on="Dataset Name", how="left")
        fcp = m["Functional Classification Provided?"]
        # as in the notebook, interval-based classes OVERWRITE the category-based class (no isna() check)
        mask = (fcp.isin(["No"]) | fcp.isna()) & (m["Score Intervals Reported?"] == "Reported")
        if mask.any():
            score = pd.to_numeric(m["auth_reported_score"], errors="coerce")
            m.loc[mask, "StandardizedClass"] = m.loc[mask].apply(
                lambda r: find_matching_class(r, score.loc[r.name]), axis=1
            )
        log(f"  Supp 3 interval reclassification applied to {mask.sum()} F9/TP53 rows")
    cls = m["StandardizedClass"].astype(str).str.upper()
    n = pd.to_numeric(m["OddsNormal"], errors="coerce")
    a = pd.to_numeric(m["OddsAbnormal"], errors="coerce")
    pts = pd.Series(np.nan, index=m.index)
    normal, abnormal = cls == "NORMAL", cls == "ABNORMAL"
    pts[normal] = np.select([n[normal] < 0.053, n[normal] < 0.23, n[normal] < 0.48], [-4, -2, -1], 0)
    pts[abnormal] = np.select(
        [a[abnormal] > 350, a[abnormal] > 18.7, a[abnormal] > 4.3, a[abnormal] > 2.1], [8, 4, 2, 1], 0
    )
    pts.index = m["_rowid"].values
    return pts.reindex(sub.index)


# ---------------- REVEL points + final classification ----------------
REVEL_GW_BINS = [  # genome-wide REVEL bins, applied to REVEL rounded to 3 decimals; 0.291-0.643 is indeterminate
    (None, 0.016, "BP4_Strong"),
    (0.017, 0.052, "BP4_Moderate+"),
    (0.053, 0.183, "BP4_Moderate"),
    (0.184, 0.290, "BP4_Supporting"),
    (0.644, 0.772, "PP3_Supporting"),
    (0.773, 0.878, "PP3_Moderate"),
    (0.879, 0.931, "PP3_Moderate+"),
    (0.932, None, "PP3_Strong"),
]
CODE_POINTS = {  # points for each evidence code
    "PP3_Supporting": 1,
    "PP3_Moderate": 2,
    "PP3_Moderate+": 3,
    "PP3_Strong": 4,
    "PP3_Very Strong": 8,
    "BP4_Supporting": -1,
    "BP4_Moderate": -2,
    "BP4_Moderate+": -3,
    "BP4_Strong": -4,
    "BP4_Very Strong": -8,
}
BP4_ORDER = ["BP4_Very Strong", "BP4_Strong", "BP4_Moderate+", "BP4_Moderate", "BP4_Supporting"]
PP3_ORDER = ["PP3_Very Strong", "PP3_Strong", "PP3_Moderate+", "PP3_Moderate", "PP3_Supporting"]
CLASS_ORDER = ["Pathogenic", "Likely Pathogenic", "Uncertain", "Likely Benign", "Benign"]


def add_revel_classification(d, supp4, keep_bp4_very_strong=False):
    """Total_Points_REVEL = Fxn_points + gene-specific REVEL points (genome-wide
    fallback, then 0), mapped to Class_REVEL, as in the notebook."""
    d = d.copy()
    revel = pd.to_numeric(d["REVEL"], errors="coerce").to_numpy()
    r3 = np.round(revel, 3)

    # genome-wide REVEL bins
    conds = []
    for lo, hi, _ in REVEL_GW_BINS:
        c = ~np.isnan(r3)
        if lo is not None:
            c &= r3 >= lo
        if hi is not None:
            c &= r3 <= hi
        conds.append(c)
    gw = np.select(conds, [b[2] for b in REVEL_GW_BINS], default="")
    d["Points_REVEL_GenomeWide"] = pd.Series(gw, index=d.index).map(CODE_POINTS)

    # gene-specific REVEL thresholds: BP4 checked strongest->weakest, then PP3
    tab = pd.read_excel(supp4, sheet_name="REVEL_gene_specific_calibration")
    tab = tab[tab["Gene"].isin(d["Gene"].unique())].drop_duplicates("Gene")
    if "BP4_Very Strong" in tab.columns and not keep_bp4_very_strong:
        log(
            "[note] Supp 4 still has a 'BP4_Very Strong' REVEL column -- the current pipeline "
            "removed it (benign predictor evidence capped at Strong), so it is ignored. "
            "Use --keep-bp4-very-strong to apply it."
        )
        tab = tab.drop(columns=["BP4_Very Strong"])
    thr_cols = [c for c in BP4_ORDER + PP3_ORDER if c in tab.columns]
    thr = d[["Gene"]].merge(tab[["Gene"] + thr_cols], on="Gene", how="left")
    conds, codes = [], []
    for c in BP4_ORDER:
        if c in thr_cols:
            conds.append(revel <= pd.to_numeric(thr[c], errors="coerce").to_numpy())
            codes.append(c)
    for c in PP3_ORDER:
        if c in thr_cols:
            conds.append(revel >= pd.to_numeric(thr[c], errors="coerce").to_numpy())
            codes.append(c)
    gs = np.select(conds, codes, default="") if conds else np.full(len(d), "")
    d["Points_REVEL_GeneSpecific"] = pd.Series(gs, index=d.index).map(CODE_POINTS)

    # REVEL points used: gene-specific, fall back to genome-wide, then 0
    d["Points_REVEL_GeneSpecific_GenomeWide"] = (
        d["Points_REVEL_GeneSpecific"].fillna(d["Points_REVEL_GenomeWide"]).fillna(0)
    )
    f = pd.to_numeric(d["points"], errors="coerce").fillna(0).astype(float)
    pr = d["Points_REVEL_GeneSpecific_GenomeWide"].astype(float)
    tot = f + pr
    d["Fxn_points"] = f
    d["Total_Points_REVEL"] = tot

    # total points -> classification
    t = tot.astype(int).to_numpy()
    d["Class_REVEL"] = np.select([t >= 10, t >= 6, t >= 0, t >= -6], CLASS_ORDER[:4], CLASS_ORDER[4])

    # flag (not exclude) functional vs predictor disagreement
    conflict = (((f > 0) & (pr < 0)) | ((f < 0) & (pr > 0))).to_numpy()
    d["Conflicting_REVEL"] = np.where(conflict, "Conflicting evidence", tot.astype(object).to_numpy())

    n_gs = d["Points_REVEL_GeneSpecific"].notna().sum()
    n_gw = (d["Points_REVEL_GeneSpecific"].isna() & d["Points_REVEL_GenomeWide"].notna()).sum()
    log(
        f"REVEL evidence: {n_gs} gene-specific, {n_gw} genome-wide fallback, "
        f"{len(d) - n_gs - n_gw} none (0 points); {int(conflict.sum())} flagged "
        f"functional-vs-REVEL conflicting"
    )
    return d


def print_class_summary(d, label):
    n = len(d)
    if n == 0:
        return
    vc = d["Class_REVEL"].value_counts()
    print(f"\n--- Final REVEL classification: {label} (n={n}) ---")
    for c in CLASS_ORDER:
        k = int(vc.get(c, 0))
        print(f"  {c:<18}{k:>8}  ({k / n * 100:5.1f}%)")
    blb = int(vc.get("Benign", 0) + vc.get("Likely Benign", 0))
    plp = int(vc.get("Pathogenic", 0) + vc.get("Likely Pathogenic", 0))
    print(f"  {'B/LB':<18}{blb:>8}  ({blb / n * 100:5.1f}%)")
    print(f"  {'P/LP':<18}{plp:>8}  ({plp / n * 100:5.1f}%)")
    print(f"  {'Resolved (not Uncertain)':<26}{blb + plp:>8}  ({(blb + plp) / n * 100:5.1f}%)")


def filter_pillar_data(
    points_input,
    category,
    supp4_input=None,
    supp3_input=None,
    chek2_input=None,
    no_excalibr_backfill=False,
    no_urn_propagation=False,
    keep_bp4_very_strong=False,
    validate_n=20000,
):
    """Run the full filtering pipeline and return (dedup, revel) dataframes --
    the deduplicated category selection, and the same with REVEL-training
    variants dropped. Raises ValueError for bad input."""
    df = pd.read_csv(points_input, dtype=str, low_memory=False)
    log(f"Loaded {len(df)} rows")
    required = [
        "Dataset",
        "Gene",
        "excalibr_points",
        "nucleotide_or_aa",
        "Chrom",
        "hg38_start",
        "ref_allele",
        "alt_allele",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{points_input} is missing required column(s): {', '.join(missing)}")
    df["Dataset"] = df["Dataset"].apply(nfc)
    df["Gene"] = df["Gene"].astype(str)

    # ---------------- functional points for every row ----------------
    df["points"] = pd.to_numeric(df["excalibr_points"], errors="coerce")
    is_op = df["Gene"].isin(OP_GENES)
    reason = (
        df["excalibr_filter_reason"] if "excalibr_filter_reason" in df.columns else pd.Series(np.nan, index=df.index)
    )

    if supp4_input and not no_excalibr_backfill:
        ex = load_excalibr(supp4_input)
        if validate_n:
            scored = df[df["points"].notna() & ~is_op]
            s = scored.sample(n=min(validate_n, len(scored)), random_state=42)
            agree = (exc_fxn(s, ex).fillna(0).values == s["points"].values).mean()
            log(
                f"[validation] Supp 4 intervals reproduce excalibr_points on {agree:.2%} "
                f"of {len(s)} sampled scored rows (low agreement = Supp 4 is stale; "
                f"rerun with --no-excalibr-backfill)"
            )
        bf = df["points"].isna() & reason.isin(UNSCORED_REASONS) & ~is_op
        df.loc[bf, "points"] = exc_fxn(df[bf], ex)
        log(
            f"Backfilled ExCALIBR points for {bf.sum()} rows excalibr left unscored, "
            f"{(df.loc[bf, 'points'].fillna(0) != 0).sum()} nonzero; by reason: "
            f"{reason[bf].value_counts().to_dict()}"
        )
    elif not no_excalibr_backfill:
        log("NOTE: no --supp4-input -- rows excalibr left unscored stay at 0")

    op_rows = is_op & ~df["Dataset"].isin(TP53_F9_EXCLUDED)
    if supp4_input and op_rows.any():
        df.loc[op_rows, "points"] = op_points(df[op_rows], supp4_input, supp3_input)
        log(f"F9/TP53 meta-dataset rows with OP_points: {df.loc[op_rows, 'points'].notna().sum()}/{op_rows.sum()}")
    df["points"] = df["points"].fillna(0)  # missing points -> 0

    # ---------------- TP53/F9 dataset restriction; keep nt/aa rows only ----------------
    df = df[~df["Dataset"].isin(TP53_F9_EXCLUDED)].copy()
    log(f"After TP53/F9 restriction: {len(df)} rows")
    nt_aa = df["nucleotide_or_aa"].isin(["nt", "aa"])
    if (~nt_aa).any():
        log(
            f"Dropping {(~nt_aa).sum()} rows with nucleotide_or_aa not nt/aa (neither nucleotide- nor amino-acid-level)"
        )
    df = df[nt_aa].copy()
    is_nt = df["nucleotide_or_aa"] == "nt"
    is_aa = ~is_nt

    for c in ["Chrom", "hg38_start", "ref_allele", "alt_allele"]:
        df["_s_" + c] = df[c].astype(str)  # compared as strings, as in the notebook
    df["_pos"] = pd.to_numeric(df["hg38_start"], errors="coerce")
    df["_aa_pos"] = pd.to_numeric(df["aa_pos"], errors="coerce")
    df["_aa_ref"] = df["aa_ref"].astype(str)
    df["_aa_alt"] = df["aa_alt"].astype(str)
    df["_tx"] = df["RefSeq Transcript ID"].str.replace(r"\.\d+$", "", regex=True).astype(str)
    pts = df["points"]
    nt_keys = ["Gene", "_s_Chrom", "_s_hg38_start", "_s_ref_allele", "_s_alt_allele"]
    aa_keys = ["Gene", "_aa_ref", "_aa_pos", "_aa_alt", "_tx"]
    g_keys = ["Gene", "_s_Chrom", "_pos", "_s_ref_allele", "_s_alt_allele"]

    # ---------------- tag splice / start-lost variants (removed later) ----------------
    sai = df[["spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"]].apply(
        pd.to_numeric, errors="coerce"
    )
    splice_variant = (sai >= 0.2).any(axis=1) | (df["simplified_consequence"] == "splice_site_variant")
    if "splice_variant" in df.columns:
        n_diff = (np.where(splice_variant, "Yes", "No") != df["splice_variant"].astype(str)).sum()
        log(
            f"[diagnostic] recomputed splice_variant differs from the file's column on {n_diff} "
            f"rows (the notebook recomputes it, so the recomputed value is used)"
        )
    unmeasured = df["splice_measure"] != "Yes"
    notes = pd.Series(np.where(splice_variant & unmeasured, "splice_variant_not_measured", ""), index=df.index)
    start_tag = is_aa & (df["simplified_consequence"] == "start_lost")
    st = start_tag.to_numpy(dtype=bool)
    nv = notes.to_numpy(dtype=object)
    nv = np.where(
        st,
        np.where(nv != "", nv + ";start_lost_variant_not_measured", "start_lost_variant_not_measured"),
        nv,
    )
    notes = pd.Series(nv, index=df.index, dtype=object)

    # ---------------- sign conflicts within nt-level and within aa-level assays ----------------
    c12 = pd.Series(False, index=df.index, dtype=bool)
    c12.loc[is_nt] = opposite_signs(df[is_nt], nt_keys, pts[is_nt]).to_numpy(dtype=bool)
    c12.loc[is_aa] = opposite_signs(df[is_aa], aa_keys, pts[is_aa], required=["_aa_pos"]).to_numpy(dtype=bool)
    c12 = c12.astype(bool)
    nv = notes.to_numpy(dtype=object)
    nv = np.where(c12.to_numpy(), np.where(nv != "", nv + ";conflicting_fxn_data", "conflicting_fxn_data"), nv)
    notes = pd.Series(nv, index=df.index, dtype=object)
    log(f"Within-resolution conflicts: {c12.sum()} rows")

    # ---------------- spread splice / REVEL-training flags across each amino-acid change ----------------
    svamino = splice_variant.copy().astype(bool)
    svamino.loc[is_aa] = group_any(df[is_aa], aa_keys, splice_variant[is_aa], required=["_aa_pos"]).to_numpy(dtype=bool)
    rt = df["REVEL_train"].astype(str).str.strip().str.lower() == "true"
    revel_train = rt.copy().astype(bool)
    revel_train.loc[is_aa] = group_any(df[is_aa], aa_keys, rt[is_aa], required=["_aa_pos"]).to_numpy(dtype=bool)

    # ---------------- sign conflicts between nt- and aa-level assays ----------------
    c3_direct = opposite_signs(df, g_keys, pts, required=["_pos"])
    if "mavedb_variant_urn" in df.columns and not no_urn_propagation:
        c3 = df["mavedb_variant_urn"].isin(df.loc[c3_direct, "mavedb_variant_urn"].dropna().tolist()).astype(bool)
        log(
            f"Cross-resolution conflicts: {c3_direct.sum()} rows directly, {c3.sum()} after URN "
            f"propagation (+{(c3 & ~c3_direct).sum()} sibling rows)"
        )
    else:
        c3 = c3_direct
        log(f"Cross-resolution conflicts: {c3.sum()} rows (no URN propagation)")
    hit = c3.to_numpy(dtype=bool) & ~notes.isin(DISALLOWED).to_numpy(dtype=bool)
    notes = pd.Series(np.where(hit, "conflicting_fxn_data", notes.to_numpy(dtype=object)), index=df.index, dtype=object)
    quirk = notes == "start_lost_variant_not_measured;conflicting_fxn_data"
    if quirk.any():
        log(
            f"[note] {quirk.sum()} rows tagged 'start_lost...;conflicting_fxn_data' -- not in "
            f"the removal list, so the notebook keeps them (replicated)"
        )

    df["_notes"], df["_svamino"], df["_revel_train"] = notes, svamino, revel_train

    # ---------------- remove SFPQ ----------------
    df = df[df["Gene"] != "SFPQ"].copy()

    # ---------------- CHEK2 QC flag ----------------
    if chek2_input:
        ch = pd.read_excel(chek2_input)
        key = ch[["hgvs_pro", "score", "Filter_CI"]].drop_duplicates(["hgvs_pro", "score"])
        tmp = (
            pd.DataFrame(
                {
                    "hgvs_pro": df["hgvs_p"].str.replace(r"^[^:]+:", "", regex=True).values,
                    "score": pd.to_numeric(df["auth_reported_score"], errors="coerce").values,
                    "_rowid": df.index,
                }
            )
            .merge(key, on=["hgvs_pro", "score"], how="left")
            .set_index("_rowid")
        )
        fail = tmp["Filter_CI"].reindex(df.index) == 1
        log(f"CHEK2 QC: {fail.sum()} rows fail, {(fail & (df['Flag'] != '*')).sum()} not already flagged upstream")
        df.loc[fail, "Flag"] = "*"

    # ---------------- remove tagged and flagged rows ----------------
    n0 = len(df)
    rm = df["_notes"].isin(DIS).to_numpy(dtype=bool) | (
        df["_svamino"].to_numpy(dtype=bool) & (df["splice_measure"] != "Yes").to_numpy(dtype=bool)
    )
    df = df[~rm]
    log(f"Removed conflicting / unmeasured splice / start-lost rows: {n0 - len(df)}")
    n0 = len(df)
    df = df[df["Flag"] != "*"]
    log(f"Removed Flag == '*' rows: {n0 - len(df)}")

    # ---------------- select category, one row per variant, drop REVEL training ----------------
    if category == "VUS":
        cat = df["clinvar_sig_2025"] == "Uncertain significance"
    elif category == "gnomAD":
        cat = df["gnomad_MAF"].notna()
    else:
        cat = df["clinvar_sig_2025"].isna() & df["gnomad_MAF"].isna()
    c = df[cat].copy()
    log(f"Rows matching {category}: {len(c)}")
    c["_abs"] = c["points"].abs()
    c = c.sort_values(["_abs", "Dataset"], ascending=[False, True], na_position="last")
    dedup = c.drop_duplicates(subset=["Gene", "_pos", "ref_allele", "alt_allele"], keep="first")
    if category == "Unobserved":
        dedup = dedup[(dedup["ref_allele"].str.len() == 1) & (dedup["alt_allele"].str.len() == 1)]
    if supp4_input and "REVEL" in dedup.columns:
        dedup = add_revel_classification(dedup, supp4_input, keep_bp4_very_strong)
    elif supp4_input:
        log("WARNING: no REVEL column -- final classification skipped")
    revel = dedup[~dedup["_revel_train"].to_numpy(dtype=bool)]

    # internal working columns (and `points`, duplicated as Fxn_points) are not written out
    drop = [col for col in dedup.columns if col.startswith("_")]
    if "Fxn_points" in dedup.columns:
        drop.append("points")
    return dedup.drop(columns=drop), revel.drop(columns=drop)


@click.command(help=__doc__)
@click.option(
    "--points-input",
    "points_input",
    default=DEFAULT_POINTS_INPUT,
    show_default=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="dataframe_with_points.csv.gz",
)
@click.option(
    "--supp4-input",
    "supp4_input",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Supplementary_Data_4.xlsx (ExCALIBR backfill + OddsPath); typically "
    "data/output/supplementary_data/Supplementary_Data_4.xlsx",
)
@click.option(
    "--supp3-input",
    "supp3_input",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Supplementary_Data_3.xlsx (interval-based F9/TP53 functional classes); typically "
    "data/input/maves/Supplementary_Data_3.xlsx",
)
@click.option(
    "--chek2-input",
    "chek2_input",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="CHEK2_Gebbia_2024.xlsx (optional; Flag already set upstream); typically "
    "data/input/maves/CHEK2_Gebbia_2024.xlsx",
)
@click.option("--category", type=click.Choice(["VUS", "gnomAD", "Unobserved"]), default="VUS", show_default=True)
@click.option(
    "--outdir",
    "outdir",
    default=DEFAULT_OUTDIR,
    show_default=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="folder for replicated_<category>_category[_REVEL].csv (created if needed)",
)
@click.option(
    "--no-excalibr-backfill",
    "no_excalibr_backfill",
    is_flag=True,
    help="ablation: leave rows excalibr didn't score at 0 (old behavior)",
)
@click.option(
    "--no-urn-propagation",
    "no_urn_propagation",
    is_flag=True,
    help="ablation: cross-resolution conflict flags only direct rows, not URN siblings",
)
@click.option(
    "--keep-bp4-very-strong",
    "keep_bp4_very_strong",
    is_flag=True,
    help="apply a BP4_Very Strong REVEL tier if Supp 4 has one (current pipeline drops it)",
)
@click.option(
    "--validate-n",
    "validate_n",
    type=int,
    default=20000,
    show_default=True,
    help="rows sampled to check Supp 4 intervals reproduce excalibr_points (0 = skip)",
)
def main(
    points_input,
    supp4_input,
    supp3_input,
    chek2_input,
    category,
    outdir,
    no_excalibr_backfill,
    no_urn_propagation,
    keep_bp4_very_strong,
    validate_n,
):
    try:
        dedup, revel = filter_pillar_data(
            points_input,
            category,
            supp4_input=supp4_input,
            supp3_input=supp3_input,
            chek2_input=chek2_input,
            no_excalibr_backfill=no_excalibr_backfill,
            no_urn_propagation=no_urn_propagation,
            keep_bp4_very_strong=keep_bp4_very_strong,
            validate_n=validate_n,
        )
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"replicated_{category}_category"
    dedup.to_csv(f"{out}.csv", index=False)
    revel.to_csv(f"{out}_REVEL.csv", index=False)
    print(f"\n=== {category}: {len(dedup)} deduplicated ===")
    print(f"=== FINAL (REVEL): {len(revel)} ({len(dedup) - len(revel)} REVEL training variants removed) ===")
    if "Class_REVEL" in revel.columns:
        print_class_summary(revel, f"{category} (REVEL set)")
    click.echo(f"Wrote {out}.csv and {out}_REVEL.csv")


if __name__ == "__main__":
    main()
