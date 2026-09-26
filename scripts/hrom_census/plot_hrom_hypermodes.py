#!/usr/bin/env python3
"""Two-panel SynTracker-style figure for HROM hypermode enrichment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
QC_DIR = SCRIPT_DIR.parent / "figure_qc"
if str(QC_DIR) not in sys.path:
    sys.path.insert(0, str(QC_DIR))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

COLORS = {
    "background": "#B9B9B9",
    "ani_top": "#2166AC",
    "structural_top": "#B2182B",
    "both_top": "#762A83",
}
LABELS = {
    "background": "All qualifying pairs",
    "ani_top": "ANI top 5%",
    "structural_top": "structural top 5%",
    "both_top": "both top 5%",
}
MODE_COLORS = {
    "structural_enriched": "#B2182B",
    "ani_enriched": "#2166AC",
    "both_enriched": "#762A83",
    "not_enriched": "#8C8C8C",
}


def save_figure(fig, stem, dpi=600):
    require_matplotlib_panel_alignment(
        fig,
        json_out=f"{stem}.alignment.json",
        overlay_svg=f"{stem}.alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, bbox_inches="tight")


def label_species(ax, table, max_labels=12):
    x = pd.to_numeric(table["log2_struct_over_ani"], errors="coerce")
    y = pd.to_numeric(table["neg_log10_min_p"], errors="coerce")
    significant = table[
        (table["struct_q"] < 0.05) | (table["ani_q"] < 0.05)
    ].copy()
    significant = significant.sort_values(
        ["neg_log10_min_p", "total_pairs"], ascending=[False, False]
    ).head(max_labels)
    x_span = np.ptp(x) if len(x) > 1 else 1.0
    y_top = y.max() if len(y) else 1.0
    for i, row in enumerate(significant.itertuples()):
        x_val = float(row.log2_struct_over_ani)
        y_val = float(row.neg_log10_min_p)
        dx = 3.0 if x_val >= 0 else -3.0
        ha = "left" if dx > 0 else "right"
        dy = 0.8 + 0.45 * (i % 2)
        ax.annotate(
            row.species,
            (x_val, y_val),
            xytext=(x_val + dx, min(y_val + dy, y_top * 1.03)),
            fontsize=5.5,
            ha=ha,
            va="bottom",
            arrowprops=dict(arrowstyle="-", lw=0.4, color="0.35",
                            shrinkA=0, shrinkB=1),
            zorder=5,
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--points", required=True,
                    help="hypermode_figure_points.tsv.gz")
    ap.add_argument("--enrichment", required=True,
                    help="hypermode_species_enrichment.tsv")
    ap.add_argument("--report", required=True,
                    help="hypermode_analysis_report.json")
    ap.add_argument("--outstem", required=True)
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
    })

    points = pd.read_csv(args.points, sep="\t")
    table = pd.read_csv(args.enrichment, sep="\t")
    with open(args.report) as fh:
        report = json.load(fh)

    fig = plt.figure(figsize=(180 / 25.4, 72 / 25.4))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0],
                             wspace=0.28, left=0.065, right=0.985,
                             top=0.90, bottom=0.14)

    # Panel a, including marginal densities as in SynTracker Fig. 4a.
    nested = outer[0].subgridspec(2, 2, width_ratios=[4, 1],
                                  height_ratios=[1, 4], hspace=0.04,
                                  wspace=0.04)
    ax_top = fig.add_subplot(nested[0, 0])
    ax_main = fig.add_subplot(nested[1, 0])
    ax_right = fig.add_subplot(nested[1, 1])

    order = ["background", "ani_top", "structural_top", "both_top"]
    for group in order:
        sub = points[points["group"] == group]
        if sub.empty:
            continue
        color = COLORS[group]
        size = 1.0 if group == "background" else 2.2
        alpha = 0.13 if group == "background" else 0.42
        ax_main.scatter(sub["structural"], sub["ani"], s=size, c=color,
                        alpha=alpha, linewidths=0, rasterized=True,
                        label=LABELS[group], zorder=2 if group == "background" else 4)
        x0, x1 = sub["structural"].min(), sub["structural"].max()
        y0, y1 = sub["ani"].min(), sub["ani"].max()
        bins_x = np.linspace(x0, x1 if x1 > x0 else x0 + 1, 45)
        bins_y = np.linspace(y0, y1 if y1 > y0 else y0 + 1, 45)
        ax_top.hist(sub["structural"], bins=bins_x, density=True,
                    histtype="step", lw=0.8, color=color, alpha=0.85,
                    zorder=3)
        ax_right.hist(sub["ani"], bins=bins_y, density=True,
                      histtype="step", lw=0.8, color=color, alpha=0.85,
                      orientation="horizontal", zorder=3)
    ax_top.set_xticklabels([])
    ax_right.set_yticklabels([])
    ax_top.spines.right.set_visible(False)
    ax_top.spines.top.set_visible(False)
    ax_right.spines.top.set_visible(False)
    ax_right.spines.right.set_visible(False)
    ax_top.set_ylabel("Density", fontsize=6)
    ax_right.set_xlabel("Density", fontsize=6)
    ax_main.set_xlabel("Syn2b structural similarity")
    ax_main.set_ylabel("skani ANI (%)")
    ax_main.legend(loc="lower right", handletextpad=0.1, borderaxespad=0.2)
    ax_main.text(0.02, 0.97, "ANI-enriched", transform=ax_main.transAxes,
                 ha="left", va="top", color=COLORS["ani_top"], fontsize=6)
    ax_main.text(0.98, 0.03, "structural-enriched",
                 transform=ax_main.transAxes, ha="right", va="bottom",
                 color=COLORS["structural_top"], fontsize=6)

    # Panel b: species enrichment in the two top-5% sets.
    ax_b = fig.add_subplot(outer[1])
    for mode, sub in table.groupby("mode"):
        ax_b.scatter(sub["log2_struct_over_ani"], sub["neg_log10_min_p"],
                     s=13 if mode == "not_enriched" else 18,
                     c=MODE_COLORS[mode], alpha=0.65 if mode == "not_enriched" else 0.88,
                     linewidths=0, label=mode.replace("_", " "), zorder=3)
    ax_b.axvline(0, color="0.25", lw=0.7, ls="--", zorder=1)
    ax_b.axhline(-np.log10(0.05), color="0.55", lw=0.6, ls=":", zorder=1)
    ax_b.text(0.02, 0.96, "structural-enriched", transform=ax_b.transAxes,
              ha="left", va="top", color=MODE_COLORS["structural_enriched"],
              fontsize=6)
    ax_b.text(0.98, 0.96, "ANI-enriched", transform=ax_b.transAxes,
              ha="right", va="top", color=MODE_COLORS["ani_enriched"],
              fontsize=6)
    ax_b.set_xlabel(r"$\log_2$ enrichment-ratio ratio (structural / ANI)")
    ax_b.set_ylabel(r"$-\log_{10}$ minimum enrichment $P$")
    ax_b.legend(loc="lower right", ncol=1, handletextpad=0.15)
    label_species(ax_b, table)

    # Panel labels in Nature style.
    fig.text(0.008, 0.955, "a", ha="left", va="top", fontsize=8, weight="bold")
    fig.text(0.625, 0.955, "b", ha="left", va="top", fontsize=8, weight="bold")

    counts = report["counts"]
    ax_main.text(0.02, 0.03,
                 f"n = {counts['kept_pairs']:,} pairs; top 5% = "
                 f"{counts['ani_top_pairs']:,} per axis",
                 transform=ax_main.transAxes, ha="left", va="bottom",
                 fontsize=5.5, color="0.25")

    # Draw the complete layout before the alignment gate.
    fig.canvas.draw()
    save_figure(fig, args.outstem, dpi=args.dpi)


if __name__ == "__main__":
    main()
