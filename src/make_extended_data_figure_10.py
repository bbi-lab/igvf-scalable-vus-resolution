#!/usr/bin/env python3
"""Build Extended Data Figure 10: classification ("resolution") rate by gnomAD
ancestry group, for both the VUS and gnomAD variant sets.

Reads `data/input/gnomAD/vus_with_ancestry_af.csv` and
`data/input/gnomAD/gnomad_with_ancestry_af.csv` -- gnomAD v4 ancestry-group
allele frequencies joined onto this project's REVEL classification
(`Class_REVEL`), one row per variant. These files are produced upstream of
this script (by `map_gnomad_ancestry_af_v2.py`, not part of this repo) and
are checked in directly rather than rebuilt by the pipeline.

For each variant set, two framings:
  present -- exome AF > 0 in that ancestry group (groups may overlap)
  unique  -- exome AF > 0 in that group and no other

Resolved = Class_REVEL is Benign, Likely Benign, Likely Pathogenic, or
Pathogenic (not "Uncertain"). Per group: n, % resolved with a 95% bootstrap
CI, % Benign/Likely Benign, % Pathogenic/Likely Pathogenic, the difference
vs. the nfe (Non-Finnish European) reference group, a two-proportion z-test
p-value vs. nfe, and its Bonferroni correction across the 8 non-reference
groups. CIs are 95% percentile bootstrap intervals (resampling variants with
replacement).

By default only exome records passing all gnomAD filters (FILTER == PASS)
are used; records failing a filter are treated as absent. --include-filtered
uses all records.

In the "present" framing a variant shared by several groups counts in each,
so groups are not independent and the z-test is approximate; "unique" groups
are fully independent.

Outputs (in --outdir, default data/output/figures/assets/extended_data_figure_10/):
  ancestry_resolution_results.csv   all numbers, one row per set/framing/group
  ancestry_resolution_<set>.png     one two-panel figure per variant set
"""

import re
from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import norm

DEFAULT_VUS_INPUT = Path("data/input/gnomAD/vus_with_ancestry_af.csv")
DEFAULT_GNOMAD_INPUT = Path("data/input/gnomAD/gnomad_with_ancestry_af.csv")
DEFAULT_OUTDIR = Path("data/output/figures/assets/extended_data_figure_10")

GROUPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]
REF_GROUP = "nfe"
N_COMP = len(GROUPS) - 1
NAMES = {
    "afr": "African/African American",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "mid": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "remaining": "Remaining individuals",
    "sas": "South Asian",
}

# ======================================================================
# statistics
# ======================================================================


def bootstrap_ci(values, n_boot, seed=42):
    arr = np.asarray(values, dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    rates = arr[rng.integers(0, len(arr), size=(n_boot, len(arr)))].mean(axis=1)
    lo, hi = np.percentile(rates, [2.5, 97.5])
    return float(lo), float(hi)


def ztest(x1, n1, x2, n2):
    if n1 == 0 or n2 == 0:
        return np.nan
    p_pool = (x1 + x2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
    if se == 0:
        return np.nan
    return 2 * (1 - norm.cdf(abs((x1 / n1 - x2 / n2) / se)))


# ======================================================================
# data
# ======================================================================


def load(path, pass_only=False):
    """Read one ancestry-mapped CSV, coercing AF_<group> columns to numeric
    (missing/unparseable -> 0). Raises ValueError if Class_REVEL or any
    AF_<group> column is absent.

    If `pass_only`, exome records that fail a gnomAD quality filter (per
    `exomes_FILTER`) have all their AF_<group> values zeroed out, i.e. are
    treated as absent from gnomAD exomes."""
    d = pd.read_csv(path, low_memory=False)
    missing = [c for c in ["Class_REVEL"] + [f"AF_{g}" for g in GROUPS] if c not in d.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    for g in GROUPS:
        d[f"AF_{g}"] = pd.to_numeric(d[f"AF_{g}"], errors="coerce").fillna(0)
    if pass_only:
        if "exomes_FILTER" not in d.columns:
            raise ValueError(f"--pass-only needs an exomes_FILTER column, not found in {path}")
        fail = d["exomes_FILTER"].notna() & (d["exomes_FILTER"] != "PASS")
        had_af = (d[[f"AF_{g}" for g in GROUPS]] > 0).any(axis=1)
        for g in GROUPS:
            d.loc[fail, f"AF_{g}"] = 0
        print(
            f"{path}: PASS only -- {fail.sum():,} variants have a non-PASS exome record; "
            f"{(fail & had_af).sum():,} of them had exome AF > 0 and are now treated as absent"
        )
    return d


def analyze(d, label, n_bootstrap, results):
    """Append one results row per (framing, group) -- plus "ALL" -- to
    `results`, and print the same summary table shown in the manuscript
    analysis. `d` is one `load()`-ed dataframe."""
    n_in = len(d)
    d = d[d["Class_REVEL"].notna()]
    present = pd.DataFrame({g: d[f"AF_{g}"] > 0 for g in GROUPS}, index=d.index)
    observed = present.any(axis=1)
    d, present = d[observed], present[observed]
    resolved = (d["Class_REVEL"] != "Uncertain").to_numpy()
    blb = d["Class_REVEL"].isin(["Benign", "Likely Benign"]).to_numpy()
    plp = d["Class_REVEL"].isin(["Pathogenic", "Likely Pathogenic"]).to_numpy()
    n_groups = present.sum(axis=1).to_numpy()

    print(
        f"\n{'#' * 110}\n# {label}: {len(d):,} of {n_in:,} variants have exome ancestry data "
        f"({(n_groups == 1).sum():,} unique to one group)\n{'#' * 110}"
    )

    for framing in ["present", "unique"]:
        masks = {g: present[g].to_numpy() & ((n_groups == 1) if framing == "unique" else True) for g in GROUPS}
        overall = (n_groups == 1) if framing == "unique" else np.ones(len(d), dtype=bool)
        ref_m = masks[REF_GROUP]
        k_ref, n_ref = int(resolved[ref_m].sum()), int(ref_m.sum())

        print(
            "\n=== "
            + label
            + " -- "
            + (
                "present in population (may overlap other populations)"
                if framing == "present"
                else "unique to population"
            )
            + " ==="
        )
        print(
            f"{'group':<11}{'n':>7}{'n_res':>7}{'%res':>7}{'95% CI':>17}{'%B/LB':>7}{'%P/LP':>7}"
            f"{'diff_vs_nfe':>13}{'p':>10}{'p_corr':>9}"
        )

        for g in GROUPS + ["ALL"]:
            m = overall if g == "ALL" else masks[g]
            n, k = int(m.sum()), int(resolved[m].sum())
            if n == 0:
                print(f"{g:<11}{0:>7}")
                continue
            rate = k / n
            lo, hi = bootstrap_ci(resolved[m], n_bootstrap)
            diff = p = p_corr = np.nan
            diff_s = "ref" if g == REF_GROUP else ""
            p_s = pc_s = ""
            if g not in (REF_GROUP, "ALL") and n_ref:
                diff = rate - k_ref / n_ref
                p = ztest(k, n, k_ref, n_ref)
                p_corr = min(p * N_COMP, 1.0) if not np.isnan(p) else np.nan
                diff_s = f"{diff * 100:+.2f}%"
                p_s, pc_s = f"{p:.3g}", f"{p_corr:.3g}"
            print(
                f"{g:<11}{n:>7}{k:>7}{rate * 100:>6.1f}%  [{lo * 100:5.1f}%, {hi * 100:5.1f}%]"
                f"{blb[m].mean() * 100:>6.1f}%{plp[m].mean() * 100:>6.1f}%"
                f"{diff_s:>13}{p_s:>10}{pc_s:>9}"
            )
            results.append(
                {
                    "variant_set": label,
                    "framing": framing,
                    "group": g,
                    "n": n,
                    "n_resolved": k,
                    "pct_resolved": rate * 100,
                    "ci_low": lo * 100,
                    "ci_high": hi * 100,
                    "pct_BLB": blb[m].mean() * 100,
                    "pct_PLP": plp[m].mean() * 100,
                    "diff_vs_nfe": diff * 100,
                    "p": p,
                    "p_corr": p_corr,
                }
            )
    print(
        f"\n(CIs: 95% bootstrap, {n_bootstrap:,} resamples; p: two-proportion z-test vs {REF_GROUP}; "
        f"p_corr: Bonferroni across {N_COMP} comparisons)"
    )


# ======================================================================
# figures
# ======================================================================
# pastel_46: base colors in GROUPS order, each blended 46% toward white
BASE_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#17becf",
    "#7f7f7f",
    "#bcbd22",
]
BLEND = 0.46
YLABEL = {"VUS": "VUS resolution rate (%)"}
DEFAULT_YLABEL = "Variants classified B/LB/LP/P (%)"
NOUN = {"VUS": "VUS", "gnomAD": "gnomAD variants"}
EDGE, SPINE, TEXT = "#4a4a4a", "#E0E4E0", "#24302A"


COLORS = {g: tuple(c * (1 - BLEND) + BLEND for c in mcolors.to_rgb(h)) for g, h in zip(GROUPS, BASE_COLORS)}

# Two `panel()` presets: DIAGNOSTIC_STYLE for `make_figures`'s large exploratory
# PNGs (sized for on-screen inspection, DejaVu Sans), CALIBRATED_STYLE for
# `build_calibrated_figures`'s exact-print-size PDFs (6.5in wide -- this
# project's standard full-text-width, matching `save_calibrated_ablation_
# figure` in ablation_variant_reclassification.py -- Arial, pdf.fonttype=42,
# same reasoning as that figure and make_extended_data_figure_6.py).
DIAGNOSTIC_STYLE = {
    "tick_fontsize": 22,
    "ref_tick_fontsize": 23,
    "ylabel_fontsize": 23,
    "title_fontsize": 28,
    "title_pad": 20,
    "star_fontsize": 28,
    "bar_linewidth": 1.2,
    "capsize": 5,
    "error_linewidth": 1.5,
    "error_capthick": 1.5,
    "legend_fontsize": 21,
    "legend_handlelength": 1.8,
    "legend_handleheight": 1.6,
    "legend_labelspacing": 0.7,
}
CALIBRATED_STYLE = {
    "tick_fontsize": 6,
    "ref_tick_fontsize": 6.5,
    "ylabel_fontsize": 6.5,
    "title_fontsize": 7,
    "title_pad": 4,
    "star_fontsize": 8,
    "bar_linewidth": 0.6,
    "capsize": 2,
    "error_linewidth": 0.7,
    "error_capthick": 0.7,
}


def tick_label(g, lowercase):
    if lowercase:
        return g
    return "Remaining" if g == "remaining" else g.upper()


def panel_title(vset, framing):
    noun = NOUN.get(vset, "variants")
    if framing == "present":
        return f"All {noun} present in population\n(may overlap other populations)"
    return f"{noun} unique to population"


def panel(ax, sub, title, ylabel, stars, colors, style, lowercase=False, bar_width=0.6, show_legend=True):
    sub = sub.set_index("group").reindex(GROUPS)
    x = np.arange(len(GROUPS))
    rate = sub["pct_resolved"].to_numpy()
    yerr = [rate - sub["ci_low"].to_numpy(), sub["ci_high"].to_numpy() - rate]
    ax.bar(
        x,
        rate,
        yerr=yerr,
        capsize=style["capsize"],
        color=[colors[g] for g in GROUPS],
        edgecolor=EDGE,
        linewidth=style["bar_linewidth"],
        width=bar_width,
        error_kw={"elinewidth": style["error_linewidth"], "ecolor": EDGE, "capthick": style["error_capthick"]},
    )
    if pd.notna(sub.loc[REF_GROUP, "pct_resolved"]):
        ax.axhline(sub.loc[REF_GROUP, "pct_resolved"], color=EDGE, linestyle="--", linewidth=1.5, alpha=0.7, zorder=0)
    if stars:
        for xi, g in zip(x, GROUPS):
            pc = sub.loc[g, "p_corr"]
            if pd.notna(pc) and pc < 0.05:
                ax.text(
                    xi,
                    sub.loc[g, "ci_high"] + 1.5,
                    "*",
                    ha="center",
                    va="bottom",
                    fontsize=style["star_fontsize"],
                    fontweight="bold",
                    color=TEXT,
                )

    ax.set_xticks(x)
    ax.set_xticklabels([tick_label(g, lowercase) for g in GROUPS], fontsize=style["tick_fontsize"], color=TEXT)
    for t in ax.get_xticklabels():
        if t.get_text() == tick_label(REF_GROUP, lowercase):
            t.set_fontweight("bold")
            t.set_fontsize(style["ref_tick_fontsize"])
    ax.tick_params(axis="y", colors=TEXT, labelsize=style["tick_fontsize"] - 1)
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel(ylabel, fontsize=style["ylabel_fontsize"], color=TEXT)
    ax.set_title(title, fontsize=style["title_fontsize"], fontweight="bold", color=TEXT, pad=style["title_pad"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(SPINE)
    ax.grid(False)

    if not show_legend:
        return

    handles = [
        Patch(
            facecolor=colors[g],
            edgecolor=EDGE,
            label=f"{NAMES[g]} (n={int(sub.loc[g, 'n']) if pd.notna(sub.loc[g, 'n']) else 0:,})",
        )
        for g in GROUPS
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=1,
        fontsize=style["legend_fontsize"],
        frameon=False,
        labelcolor=TEXT,
        handlelength=style["legend_handlelength"],
        handleheight=style["legend_handleheight"],
        labelspacing=style["legend_labelspacing"],
    )


def make_figures(results, outdir, width, height, no_stars, ci_label, lowercase_labels=False):
    plt.rcParams["font.family"] = "DejaVu Sans"
    colors = COLORS
    outdir.mkdir(parents=True, exist_ok=True)
    for vset in results["variant_set"].unique():
        d = results[results["variant_set"] == vset]
        ylab = YLABEL.get(vset, DEFAULT_YLABEL) + (f"\n{ci_label}" if ci_label else "")
        fig, axes = plt.subplots(1, 2, figsize=(width, height), facecolor="white")
        for ax, framing in zip(axes, ["present", "unique"]):
            panel(
                ax,
                d[d["framing"] == framing],
                panel_title(vset, framing),
                ylab,
                not no_stars,
                colors,
                DIAGNOSTIC_STYLE,
                lowercase=lowercase_labels,
            )
        plt.tight_layout()
        stem = "ancestry_resolution_" + re.sub(r"[^A-Za-z0-9]+", "_", vset).strip("_")
        path = outdir / f"{stem}.png"
        plt.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {path}")


# Standard full-text-width for this project's calibrated (exact-print-size)
# figures -- same width as save_calibrated_ablation_figure's Extended Data
# Figure 10.
CALIBRATED_WIDTH_IN = 6.5
CALIBRATED_HEIGHT_IN = 3.0


def build_calibrated_figures(results, outdir, no_stars, ci_label, lowercase_labels=False):
    """Write one 6.5in-wide, exact-print-size two-panel PDF per variant set
    (ancestry_resolution_<set>_calibrated.pdf) -- the same present/unique
    panels as `make_figures`'s diagnostic PNGs, at print font sizes (Arial,
    pdf.fonttype=42, same reasoning as `save_calibrated_ablation_figure` in
    ablation_variant_reclassification.py and make_extended_data_figure_6.py).

    Unlike the diagnostic panels, these omit the per-bar-color legend: at
    this width the x-axis group-code ticks already identify each bar, so a
    nine-entry "<full name> (n=...)" legend would cost more vertical space
    than it adds -- the exact n per group/framing is already in
    ancestry_resolution_results.csv.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"font.family": "Arial", "font.size": 7, "pdf.fonttype": 42}):
        for vset in results["variant_set"].unique():
            d = results[results["variant_set"] == vset]
            ylab = YLABEL.get(vset, DEFAULT_YLABEL) + (f"\n{ci_label}" if ci_label else "")
            fig, axes = plt.subplots(
                1, 2, figsize=(CALIBRATED_WIDTH_IN, CALIBRATED_HEIGHT_IN), facecolor="white", dpi=300
            )
            for ax, framing in zip(axes, ["present", "unique"]):
                panel(
                    ax,
                    d[d["framing"] == framing],
                    panel_title(vset, framing),
                    ylab,
                    not no_stars,
                    COLORS,
                    CALIBRATED_STYLE,
                    lowercase=lowercase_labels,
                    bar_width=0.7,
                    show_legend=False,
                )
            fig.tight_layout(w_pad=1.5)
            stem = "ancestry_resolution_" + re.sub(r"[^A-Za-z0-9]+", "_", vset).strip("_")
            path = outdir / f"{stem}_calibrated.pdf"
            fig.savefig(path, bbox_inches="tight")
            plt.close(fig)
            print(f"Wrote {path}")


# ======================================================================
@click.command(help=__doc__)
@click.option(
    "--vus-input", "vus_input", default=DEFAULT_VUS_INPUT, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--gnomad-input",
    "gnomad_input",
    default=DEFAULT_GNOMAD_INPUT,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--outdir",
    "outdir",
    default=DEFAULT_OUTDIR,
    type=click.Path(file_okay=False, path_type=Path),
    show_default=True,
    help="folder for all outputs (created if needed)",
)
@click.option(
    "--output",
    "output",
    default="ancestry_resolution_results.csv",
    help="results CSV file name, written inside --outdir",
)
@click.option("--n-bootstrap", "n_bootstrap", type=int, default=1000, show_default=True)
@click.option("--no-stars", "no_stars", is_flag=True, help="omit the * marking groups with Bonferroni p < 0.05 vs nfe")
@click.option(
    "--ci-label", "ci_label", default="", help='optional second y-axis line, e.g. "error bars: 95%% bootstrap CI"'
)
@click.option("--width", "width", type=float, default=28, show_default=True, help="figure width in inches")
@click.option("--height", "height", type=float, default=14.5, show_default=True, help="figure height in inches")
@click.option("--no-figures", "no_figures", is_flag=True, help="statistics only")
@click.option(
    "--lowercase-labels",
    "lowercase_labels",
    is_flag=True,
    help="x-axis labels as gnomAD codes (afr, nfe, ...) instead of AFR, NFE, ...",
)
@click.option(
    "--include-filtered",
    "include_filtered",
    is_flag=True,
    help="also use exome records that fail gnomAD quality filters (default: PASS sites only)",
)
def main(
    vus_input,
    gnomad_input,
    outdir,
    output,
    n_bootstrap,
    no_stars,
    ci_label,
    width,
    height,
    no_figures,
    lowercase_labels,
    include_filtered,
):
    pass_only = not include_filtered
    click.echo(
        "gnomAD exome records used: "
        + (
            "PASS sites only (default)"
            if pass_only
            else "all records, including those failing gnomAD filters (--include-filtered)"
        )
    )
    try:
        results = []
        analyze(load(vus_input, pass_only), "VUS", n_bootstrap, results)
        analyze(load(gnomad_input, pass_only), "gnomAD", n_bootstrap, results)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    outdir.mkdir(parents=True, exist_ok=True)
    res = pd.DataFrame(results)
    csv_path = outdir / output
    res.to_csv(csv_path, index=False)
    click.echo(f"\nWrote {csv_path}")
    if not no_figures:
        make_figures(res, outdir, width, height, no_stars, ci_label, lowercase_labels)
        build_calibrated_figures(res, outdir, no_stars, ci_label, lowercase_labels)


if __name__ == "__main__":
    main()
