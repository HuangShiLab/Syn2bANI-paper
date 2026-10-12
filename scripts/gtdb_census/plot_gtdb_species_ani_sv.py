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


def save_figure(fig, stem, dpi):
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
        ("all", "All effective pairs", "#4D4D4D"),
        ("HQ_comp90_contam5", "Completeness ≥ 90%\ncontamination ≤ 5%", "#2166AC"),
        ("min_N50_ge10kb", "Both N50 ≥ 10 kb", "#B2182B"),
        ("HQ_N50_ge10kb", "HQ and N50 ≥ 10 kb", "#762A83"),
    ]

    fig = plt.figure(figsize=(180 / 25.4, 72 / 25.4))
    outer = fig.add_gridspec(
        1, 2, width_ratios=[1.28, 1.0], wspace=0.30,
        left=0.065, right=0.985, top=0.88, bottom=0.17,
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
    fractions = []
    labels = []
    colors = []
    for i, (key, label, color) in enumerate(strata):
        row = quality.loc[key]
        pairs = int(row["pairs"])
        ge2 = int(row["pairs_ge2"])
        fraction = ge2 / pairs if pairs else np.nan
        fractions.append(fraction)
        labels.append(label)
        colors.append(color)
        ax_b.text(
            i, fraction + 0.025, f"{fraction * 100:.1f}%\nn={pairs:,}",
            ha="center", va="bottom", fontsize=5.5,
        )
    bars = ax_b.bar(xpos, fractions, width=0.60, color=colors, edgecolor="white", lw=0.3)
    ax_b.set_xticks(xpos, labels)
    ax_b.set_ylim(0, 1.03)
    ax_b.set_ylabel("Pairs with ≥2 breakpoints")
    ax_b.set_xlabel("GTDB quality stratum")

    fig.text(0.006, 0.96, "a", ha="left", va="top", fontsize=8, weight="bold")
    fig.text(0.625, 0.96, "b", ha="left", va="top", fontsize=8, weight="bold")
    fig.canvas.draw()
    save_figure(fig, args.outstem, args.dpi)

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
                "fraction_ge2": float(quality.loc[key, "fraction_ge2"]),
            }
            for key, _, _ in strata
        },
    }
    Path(f"{args.outstem}.qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
