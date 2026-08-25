#!/usr/bin/env python3
"""Analyze whether GTDB genome quality affects Syn2bANI ANI accuracy vs ANIm.

Joins the GTDB-R207 held-out benchmark with GTDB metadata (completeness,
contamination, contig count, mean contig length) and reports MAE/bias/r by
quality bins.
"""

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
                               "contig_count", "mean_contig_length", "genome_size"])
    ar = pd.read_csv(DATA / "gtdb_metadata" / "ar53_metadata_r207.tsv", sep="\t",
                     usecols=["accession", "checkm_completeness", "checkm_contamination",
                              "contig_count", "mean_contig_length", "genome_size"])
    meta = pd.concat([bac, ar], ignore_index=True)
    meta["accession"] = meta["accession"].astype(str).str.strip()
    # GTDB metadata prefixes accessions with GB_/RS_; strip to match NCBI-style IDs.
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
        "median_abs_err": err.abs().median(),
        "bias": err.mean(),
        "r": df[col].corr(df["anim_ani"]),
    })


def plot_quality_bins(summary_list, ylabel, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    titles = ["MAE", "Bias", "Pearson r"]
    keys = ["MAE", "bias", "r"]
    for ax, title, key in zip(axes, titles, keys):
        x = np.arange(len(summary_list))
        vals = summary_list[key].values
        ax.bar(x, vals, color="steelblue")
        ax.set_xticks(x)
        ax.set_xticklabels(summary_list.index, rotation=30, ha="right", fontsize=8)
        ax.set_title(title)
        ax.axhline(0, color="black", lw=0.5)
        if key == "MAE":
            ax.set_ylabel("ANI error (percentage points)")
        # annotate n on top
        for i, (v, n) in enumerate(zip(vals, summary_list["n"].values)):
            if not np.isnan(v):
                ax.text(i, v, f"n={n}", ha="center", va="bottom", fontsize=6)
    plt.suptitle(ylabel)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    s2b = pd.read_csv(RES / "s2b_50k.tsv", sep="\t")
    truth = pd.read_csv(RES / "truth_50k.tsv", sep="\t")
    meta = load_metadata()

    df = s2b.merge(truth, on="pairid", how="inner")
    q_acc, r_acc = zip(*df["pairid"].apply(parse_pairid))
    df["q_acc"] = q_acc
    df["r_acc"] = r_acc

    # Use minimum completeness / maximum contamination across the pair
    df = df.merge(meta.rename(columns=lambda c: f"q_{c}" if c != "accession" else "q_acc"),
                  on="q_acc", how="left")
    df = df.merge(meta.rename(columns=lambda c: f"r_{c}" if c != "accession" else "r_acc"),
                  on="r_acc", how="left")

    df["min_completeness"] = df[["q_checkm_completeness", "r_checkm_completeness"]].min(axis=1)
    df["max_contamination"] = df[["q_checkm_contamination", "r_checkm_contamination"]].max(axis=1)
    df["max_contig_count"] = df[["q_contig_count", "r_contig_count"]].max(axis=1)
    df["min_mean_contig_len"] = df[["q_mean_contig_length", "r_mean_contig_length"]].min(axis=1)
    df["min_genome_size"] = df[["q_genome_size", "r_genome_size"]].min(axis=1)

    # Use calibrated estimate as primary
    est_col = "ani_cal"
    df = df[df[est_col].notna()].copy()

    reports = []

    # Completeness bins
    comp_bins = [0, 70, 80, 90, 95, 100]
    df["comp_bin"] = pd.cut(df["min_completeness"], comp_bins, right=False)
    comp_summary = df.groupby("comp_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    comp_summary.index = [f"[{int(interval.left)},{int(interval.right)})" for interval in comp_summary.index]
    reports.append(("Completeness (min of pair)", comp_summary))
    plot_quality_bins(comp_summary, "ANI error by CheckM completeness", OUT_FIG / "gtdb_mae_by_completeness.png")

    # Contamination bins
    cont_bins = [0, 1, 2, 5, 10, 100]
    df["cont_bin"] = pd.cut(df["max_contamination"], cont_bins, right=False)
    cont_summary = df.groupby("cont_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cont_summary.index = [f"[{interval.left},{interval.right})" for interval in cont_summary.index]
    reports.append(("Contamination (max of pair)", cont_summary))
    plot_quality_bins(cont_summary, "ANI error by CheckM contamination", OUT_FIG / "gtdb_mae_by_contamination.png")

    # Contig count bins (log)
    cc_bins = [0, 50, 100, 200, 500, 10000]
    df["cc_bin"] = pd.cut(df["max_contig_count"], cc_bins, right=False)
    cc_summary = df.groupby("cc_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cc_summary.index = [f"[{int(interval.left)},{int(interval.right)})" for interval in cc_summary.index]
    reports.append(("Contig count (max of pair)", cc_summary))
    plot_quality_bins(cc_summary, "ANI error by contig count", OUT_FIG / "gtdb_mae_by_contig_count.png")

    # Mean contig length bins
    cl_bins = [0, 5e3, 1e4, 5e4, 1e5, 1e9]
    df["cl_bin"] = pd.cut(df["min_mean_contig_len"], cl_bins, right=False)
    cl_summary = df.groupby("cl_bin", observed=False).apply(lambda g: summarize(g, est_col), include_groups=False)
    cl_summary.index = [f"[{int(interval.left/1000)},{int(interval.right/1000)}) kb" for interval in cl_summary.index]
    reports.append(("Mean contig length (min of pair)", cl_summary))
    plot_quality_bins(cl_summary, "ANI error by mean contig length", OUT_FIG / "gtdb_mae_by_contig_length.png")

    # Overall
    overall = summarize(df, est_col)

    # Write report
    report_path = RES / "QUALITY_VS_MAE_REPORT.md"
    with open(report_path, "w") as fh:
        fh.write("# GTDB-R207 genome quality vs Syn2bANI ANI accuracy\n\n")
        fh.write(f"Estimator: `{est_col}`\n\n")
        fh.write(f"Overall: n = {overall['n']:.0f}, MAE = {overall['MAE']:.4f}, "
                 f"bias = {overall['bias']:.4f}, r = {overall['r']:.4f}\n\n")
        for title, summ in reports:
            fh.write(f"## {title}\n\n")
            fh.write(summ.round(4).to_string())
            fh.write("\n\n")
        fh.write("Figures:\n")
        fh.write("- `figures/report/gtdb_mae_by_completeness.png`\n")
        fh.write("- `figures/report/gtdb_mae_by_contamination.png`\n")
        fh.write("- `figures/report/gtdb_mae_by_contig_count.png`\n")
        fh.write("- `figures/report/gtdb_mae_by_contig_length.png`\n")

    print(f"Wrote {report_path}")
    print("\nOverall:", overall.to_dict())
    for title, summ in reports:
        print(f"\n{title}")
        print(summ.round(4).to_string())


if __name__ == "__main__":
    main()
