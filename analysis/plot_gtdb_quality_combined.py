#!/usr/bin/env python3
"""Combined figure: Syn2bANI ANI accuracy vs GTDB-R207 genome quality metrics."""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results" / "gtdb50k"
OUT_FIG = ROOT / "figures" / "report"
OUT_FIG.mkdir(parents=True, exist_ok=True)


def load_metadata():
    bac = pd.read_csv(DATA / "gtdb_metadata" / "bac120_metadata_r207.tsv", sep="\t",
                      usecols=["accession", "checkm_completeness", "checkm_contamination",
                               "contig_count", "mean_contig_length"])
    ar = pd.read_csv(DATA / "gtdb_metadata" / "ar53_metadata_r207.tsv", sep="\t",
                     usecols=["accession", "checkm_completeness", "checkm_contamination",
                              "contig_count", "mean_contig_length"])
    meta = pd.concat([bac, ar], ignore_index=True)
    meta["accession"] = meta["accession"].astype(str).str.strip()
    meta["accession"] = meta["accession"].str.replace(r"^(GB_|RS_)", "", regex=True)
    return meta


def parse_pairid(pairid):
    parts = pairid.split("__")
    return parts[0], parts[1]


def summarize(df, col):
    df = df.dropna(subset=[col, "anim_ani"])
    if len(df) == 0:
        return pd.Series({"n": 0, "MAE": np.nan, "bias": np.nan, "r": np.nan})
    err = df[col] - df["anim_ani"]
    return pd.Series({
        "n": len(df),
        "MAE": err.abs().mean(),
        "bias": err.mean(),
        "r": df[col].corr(df["anim_ani"]),
    })


def plot_axis(ax, summary, title, xlabel, ylim=None):
    x = np.arange(len(summary))
    ax.bar(x, summary["MAE"].values, color="steelblue", width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(summary.index, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("MAE (ANI points)")
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    if ylim is not None:
        ax.set_ylim(ylim)
    for i, (v, n) in enumerate(zip(summary["MAE"].values, summary["n"].values)):
        if not np.isnan(v):
            ax.text(i, v, f"n={int(n)}", ha="center", va="bottom", fontsize=6)


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

    est_col = "ani_cal"
    df = df[df[est_col].notna()].copy()

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    comp_bins = [0, 70, 80, 90, 95, 100]
    df["comp_bin"] = pd.cut(df["min_completeness"], comp_bins, right=False)
    comp_summary = df.groupby("comp_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    comp_summary.index = [f"[{int(interval.left)},{int(interval.right)})" for interval in comp_summary.index]
    plot_axis(axes[0, 0], comp_summary, "(a) CheckM completeness (min of pair)", "Completeness (%)", ylim=(0, 1.0))

    cont_bins = [0, 1, 2, 5, 10, 100]
    df["cont_bin"] = pd.cut(df["max_contamination"], cont_bins, right=False)
    cont_summary = df.groupby("cont_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cont_summary.index = [f"[{interval.left},{interval.right})" for interval in cont_summary.index]
    plot_axis(axes[0, 1], cont_summary, "(b) CheckM contamination (max of pair)", "Contamination (%)", ylim=(0, 1.0))

    cc_bins = [0, 50, 100, 200, 500, 10000]
    df["cc_bin"] = pd.cut(df["max_contig_count"], cc_bins, right=False)
    cc_summary = df.groupby("cc_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cc_summary.index = [f"[{int(interval.left)},{int(interval.right)})" for interval in cc_summary.index]
    plot_axis(axes[1, 0], cc_summary, "(c) Contig count (max of pair)", "Contigs", ylim=(0, 1.0))

    cl_bins = [0, 5e3, 1e4, 5e4, 1e5, 1e9]
    df["cl_bin"] = pd.cut(df["min_mean_contig_len"], cl_bins, right=False)
    cl_summary = df.groupby("cl_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cl_summary.index = [f"[{int(interval.left/1000)},{int(interval.right/1000)})" for interval in cl_summary.index]
    plot_axis(axes[1, 1], cl_summary, "(d) Mean contig length (min of pair)", "Mean contig length (kb)", ylim=(0, 1.0))

    plt.suptitle("Syn2bANI calibrated ANI error is stable across GTDB-R207 genome quality", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / "gtdb_quality_vs_mae_combined.png", dpi=300)
    fig.savefig(OUT_FIG / "gtdb_quality_vs_mae_combined.pdf")
    plt.close(fig)

    print(f"Wrote {OUT_FIG / 'gtdb_quality_vs_mae_combined.png'}")


if __name__ == "__main__":
    main()
