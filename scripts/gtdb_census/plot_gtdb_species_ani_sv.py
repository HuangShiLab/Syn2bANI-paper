#!/usr/bin/env python3
"""Two-panel GTDB census figure: ANI-SV spread and quality sensitivity."""

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


def save_figure(fig, stem, dpi, row_group_ids):
    require_matplotlib_panel_alignment(
        fig,
        json_out=f"{stem}.alignment.json",
        overlay_svg=f"{stem}.alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        row_groups=[{"id": "figure-row", "panels": row_group_ids}],
        exemptions=[{
            "panels": row_group_ids,
            "checks": ["panel-width"],
            "reason": "Panel a carries a colorbar and is intentionally wider than the quality-stratum bar panel.",
        }],
        strict=True,
    )
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--joint", required=True)
    ap.add_argument("--quality-sensitivity", required=True)
    ap.add_argument("--outstem", required=True)
    ap.add_argument("--min-pairs", type=int, default=100)
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

    joint = pd.read_csv(args.joint, sep="\t")
    complete = joint[
        (joint["ani_sv_layer"] == "complete")
        & (pd.to_numeric(joint["structural_coverage"], errors="coerce") == 1)
        & (pd.to_numeric(joint["ani_coverage"], errors="coerce") == 1)
        & (pd.to_numeric(joint["expected_pairs"], errors="coerce") >= args.min_pairs)
    ].copy()
    complete["mean_ani"] = pd.to_numeric(complete["mean_skani_ani"], errors="coerce")
    complete["sv_fraction"] = pd.to_numeric(complete["sv_ge2_fraction"], errors="coerce")
    complete = complete.dropna(subset=["mean_ani", "sv_fraction"])

    quality = pd.read_csv(args.quality_sensitivity, sep="\t").set_index("stratum")
    strata = [
        ("all", "All", "#4D4D4D"),
        ("HQ_comp90_contam5", "HQ\ncomp ≥ 90%\ncontam ≤ 5%", "#2166AC"),
        ("min_N50_ge10kb", "N50\n≥ 10 kb", "#B2182B"),
        ("HQ_N50_ge10kb", "HQ +\nN50\n≥ 10 kb", "#762A83"),
    ]

    fig = plt.figure(figsize=(180 / 25.4, 80 / 25.4))
    outer = fig.add_gridspec(
        1, 2, width_ratios=[1.28, 1.0], wspace=0.30,
        left=0.065, right=0.985, top=0.89, bottom=0.26,
    )

    ax_a = fig.add_subplot(outer[0])
    x = complete["mean_ani"].to_numpy()
    y = complete["sv_fraction"].to_numpy()
    hb = ax_a.hexbin(
        x, y, gridsize=(42, 18), bins="log", mincnt=1,
        cmap="viridis", linewidths=0, rasterized=True,
    )
    cb = fig.colorbar(hb, ax=ax_a, pad=0.012, extend="max")
    cb.set_label("Species pairs", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax_a.axhline(0.5, color="0.55", lw=0.6, ls="--")
    ax_a.set_xlabel("Mean skani ANI within species (%)")
    ax_a.set_ylabel("Fraction of pairs with ≥2 breakpoints")
    ax_a.text(
        0.02, 0.95,
        f"{len(complete):,} complete species\n≥ {args.min_pairs:,} expected pairs",
        transform=ax_a.transAxes, ha="left", va="top", fontsize=5.8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1},
    )

    ax_b = fig.add_subplot(outer[1])
    xpos = np.arange(len(strata))
    reported_fractions = []
    informative_fractions = []
    labels = []
    colors = []
    for i, (key, label, color) in enumerate(strata):
        row = quality.loc[key]
        reported_fractions.append(float(row["fraction_ge2_reported"]))
        informative_fractions.append(float(row["fraction_ge2_informative"]))
        labels.append(label)
        colors.append(color)
    width = 0.36
    reported_bars = ax_b.bar(
        xpos - width / 2, reported_fractions, width,
        color="0.72", edgecolor="white", lw=0.3,
        label="Reported pairs",
    )
    informative_bars = ax_b.bar(
        xpos + width / 2, informative_fractions, width,
        color=colors, edgecolor="white", lw=0.3,
        label="Breakpoint-informative pairs",
    )
    for bars in (reported_bars, informative_bars):
        for rect in bars:
            height = rect.get_height()
            ax_b.text(
                rect.get_x() + rect.get_width() / 2,
                height + 0.018,
                f"{height * 100:.1f}%",
                ha="center", va="bottom", fontsize=5.0, rotation=90,
            )
    ax_b.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
        handlelength=1.2, columnspacing=1.2,
    )
    ax_b.set_xticks(xpos, labels)
    ax_b.set_ylim(0, 1.08)
    ax_b.set_ylabel("Pairs with ≥2 breakpoints")
    ax_b.set_xlabel("GTDB quality stratum")

    fig.text(0.006, 0.96, "a", ha="left", va="top", fontsize=8, weight="bold")
    fig.text(0.625, 0.96, "b", ha="left", va="top", fontsize=8, weight="bold")
    fig.canvas.draw()
    save_figure(fig, args.outstem, args.dpi, ["a", "b"])

    qa = {
        "joint_table": args.joint,
        "complete_species": len(complete),
        "minimum_expected_pairs": args.min_pairs,
        "ani_range": [float(x.min()), float(x.max())] if len(x) else [],
        "sv_fraction_range": [float(y.min()), float(y.max())] if len(y) else [],
        "quality_strata": {
            key: {
                "pairs": int(quality.loc[key, "pairs"]),
                "pairs_ge2": int(quality.loc[key, "pairs_ge2"]),
                "fraction_ge2_reported": float(quality.loc[key, "fraction_ge2_reported"]),
                "fraction_ge2_informative": float(quality.loc[key, "fraction_ge2_informative"]),
            }
            for key, _, _ in strata
        },
    }
    Path(f"{args.outstem}.qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
