#!/usr/bin/env python3
"""Generate Supplementary Figure S4: GTDB-R207 genome quality vs MAE.

Four panels show calibrated Syn2bANI v5 MAE binned by:
  (a) minimum pair completeness,
  (b) maximum pair contamination,
  (c) maximum pair contig count,
  (d) minimum mean contig length.
Error bars are 95% bootstrap confidence intervals.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

set_publication_style()

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results" / "gtdb50k"
OUT = ROOT / "paper" / "figures" / "supplementary" / "fig_s4_gtdb_quality_vs_mae"

EST_COL = "ani_cal"
RNG_SEED = 42
N_BOOT = 1000


def load_metadata():
    cols = ["accession", "checkm_completeness", "checkm_contamination",
            "contig_count", "mean_contig_length"]
    bac = pd.read_csv(DATA / "gtdb_metadata" / "bac120_metadata_r207.tsv", sep="\t", usecols=cols)
    ar = pd.read_csv(DATA / "gtdb_metadata" / "ar53_metadata_r207.tsv", sep="\t", usecols=cols)
    meta = pd.concat([bac, ar], ignore_index=True)
    meta["accession"] = meta["accession"].astype(str).str.strip()
    # Strip GTDB GB_/RS_ prefixes to match NCBI-style accessions in pair IDs.
    meta["accession"] = meta["accession"].str.replace(r"^(GB_|RS_)", "", regex=True)
    return meta


def parse_pairid(pairid):
    q, r = pairid.split("__", 1)
    return q, r


def bootstrap_mae(errors, rng, n_boot=N_BOOT):
    """Return (low, high) 95% CI for MAE by resampling signed errors."""
    if len(errors) == 0:
        return np.nan, np.nan
    errors = np.asarray(errors)
    boot = []
    for _ in range(n_boot):
        sample = rng.choice(errors, size=len(errors), replace=True)
        boot.append(np.abs(sample).mean())
    boot = np.sort(boot)
    return boot[int(0.025 * n_boot)], boot[int(0.975 * n_boot)]


def summarize_bin(df, col, rng):
    sub = df.dropna(subset=[col, "anim_ani"])
    if len(sub) == 0:
        return {"n": 0, "mae": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    errors = sub[col].values - sub["anim_ani"].values
    mae = np.abs(errors).mean()
    ci_low, ci_high = bootstrap_mae(errors, rng)
    return {"n": len(sub), "mae": mae, "ci_low": ci_low, "ci_high": ci_high}


def plot_panel(ax, summary, title, xlabel, bin_labels, ylim=(0, 1.0)):
    x = np.arange(len(summary))
    mae_vals = np.array([s["mae"] for s in summary])
    lows = np.array([s["ci_low"] for s in summary])
    highs = np.array([s["ci_high"] for s in summary])
    err_low = mae_vals - lows
    err_high = highs - mae_vals

    ax.bar(x, mae_vals, color=COLORS["orange"], width=0.65, edgecolor="white", linewidth=0.5)
    ax.errorbar(x, mae_vals, yerr=[err_low, err_high], fmt="none",
                ecolor="black", capsize=2.5, capthick=0.8, elinewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("MAE (ANI points)")
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=8, pad=6)
    ax.set_ylim(ylim)
    for i, s in enumerate(summary):
        if s["n"] > 0:
            ax.text(i, mae_vals[i] + (ylim[1] - ylim[0]) * 0.02,
                    f"n={s['n']}", ha="center", va="bottom", fontsize=5.5)


def main():
    s2b = pd.read_csv(RES / "s2b_50k.tsv", sep="\t")
    truth = pd.read_csv(RES / "truth_50k.tsv", sep="\t")
    meta = load_metadata()

    df = s2b.merge(truth, on="pairid", how="inner")
    q_acc, r_acc = zip(*df["pairid"].apply(parse_pairid))
    df["q_acc"] = q_acc
    df["r_acc"] = r_acc

    df = df.merge(meta.rename(columns=lambda c: f"q_{c}" if c != "accession" else "q_acc"),
                  on="q_acc", how="left")
    df = df.merge(meta.rename(columns=lambda c: f"r_{c}" if c != "accession" else "r_acc"),
                  on="r_acc", how="left")

    df["min_completeness"] = df[["q_checkm_completeness", "r_checkm_completeness"]].min(axis=1)
    df["max_contamination"] = df[["q_checkm_contamination", "r_checkm_contamination"]].max(axis=1)
    df["max_contig_count"] = df[["q_contig_count", "r_contig_count"]].max(axis=1)
    df["min_mean_contig_len"] = df[["q_mean_contig_length", "r_mean_contig_length"]].min(axis=1)

    df = df[df[EST_COL].notna()].copy()

    rng = np.random.default_rng(RNG_SEED)

    # Completeness bins
    comp_bins = [0, 70, 80, 90, 95, 100]
    df["comp_bin"] = pd.cut(df["min_completeness"], comp_bins, right=False)
    comp_summary = []
    comp_labels = []
    for interval in df["comp_bin"].cat.categories:
        sub = df[df["comp_bin"] == interval]
        comp_summary.append(summarize_bin(sub, EST_COL, rng))
        comp_labels.append(f"[{int(interval.left)},{int(interval.right)})")

    # Contamination bins
    cont_bins = [0, 1, 2, 5, 10, 100]
    df["cont_bin"] = pd.cut(df["max_contamination"], cont_bins, right=False)
    cont_summary = []
    cont_labels = []
    for interval in df["cont_bin"].cat.categories:
        sub = df[df["cont_bin"] == interval]
        cont_summary.append(summarize_bin(sub, EST_COL, rng))
        cont_labels.append(f"[{interval.left:g},{interval.right:g})")

    # Contig count bins
    cc_bins = [0, 50, 100, 200, 500, 10000]
    df["cc_bin"] = pd.cut(df["max_contig_count"], cc_bins, right=False)
    cc_summary = []
    cc_labels = []
    for interval in df["cc_bin"].cat.categories:
        sub = df[df["cc_bin"] == interval]
        cc_summary.append(summarize_bin(sub, EST_COL, rng))
        cc_labels.append(f"[{int(interval.left)},{int(interval.right)})")

    # Mean contig length bins
    cl_bins = [0, 5e3, 1e4, 5e4, 1e5, 1e9]
    df["cl_bin"] = pd.cut(df["min_mean_contig_len"], cl_bins, right=False)
    cl_summary = []
    cl_labels = []
    for interval in df["cl_bin"].cat.categories:
        sub = df[df["cl_bin"] == interval]
        cl_summary.append(summarize_bin(sub, EST_COL, rng))
        if interval.right >= 1e8:
            cl_labels.append(f"[{int(interval.left/1000)},∞)")
        else:
            cl_labels.append(f"[{int(interval.left/1000)},{int(interval.right/1000)})")

    overall_err = df[EST_COL] - df["anim_ani"]
    overall_mae = overall_err.abs().mean()
    overall_r = df[EST_COL].corr(df["anim_ani"])
    print(f"Overall: n = {len(df)}, MAE = {overall_mae:.4f}, r = {overall_r:.4f}")

    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.85))
    plot_panel(axes[0, 0], comp_summary, "CheckM completeness (min of pair)",
               "Completeness (%)", comp_labels)
    plot_panel(axes[0, 1], cont_summary, "CheckM contamination (max of pair)",
               "Contamination (%)", cont_labels)
    plot_panel(axes[1, 0], cc_summary, "Contig count (max of pair)",
               "Contigs", cc_labels)
    plot_panel(axes[1, 1], cl_summary, "Mean contig length (min of pair)",
               "Mean contig length (kb)", cl_labels)

    label_panel(axes[0, 0], "a")
    label_panel(axes[0, 1], "b")
    label_panel(axes[1, 0], "c")
    label_panel(axes[1, 1], "d")

    fig.suptitle(
        f"Calibrated Syn2bANI v5 accuracy is stable across GTDB-R207 genome quality "
        f"(n = {len(df)}, MAE = {overall_mae:.3f}, r = {overall_r:.3f})",
        fontsize=9, y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
