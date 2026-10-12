#!/usr/bin/env python3
"""Two-panel figure for OAPGC oral-airway SNP and SV effects."""

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


CATEGORY_COLORS = {
    "SNP-distance enriched only": "#2166AC",
    "SV/breakpoint enriched only": "#B2182B",
    "Both axes enriched": "#762A83",
    "Not significant / depleted": "#8C8C8C",
}


def classify(row):
    snp_high = float(row["SNP_q"]) < 0.05 and float(row["SNP_auc"]) > 0.5
    bp_high = float(row["BP_q"]) < 0.05 and float(row["BP_auc"]) > 0.5
    if snp_high and bp_high:
        return "Both axes enriched"
    if snp_high:
        return "SNP-distance enriched only"
    if bp_high:
        return "SV/breakpoint enriched only"
    return "Not significant / depleted"


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
    ap.add_argument("--species", required=True)
    ap.add_argument("--overall-all", required=True)
    ap.add_argument("--overall-n50-10kb", required=True)
    ap.add_argument("--overall-n50-100kb", required=True)
    ap.add_argument("--quality", required=True,
                    help="quality_stratified_summary.tsv with pair counts")
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

    species = pd.read_csv(args.species, sep="\t")
    species["category"] = species.apply(classify, axis=1)
    category_order = [
        "SNP-distance enriched only",
        "SV/breakpoint enriched only",
        "Both axes enriched",
        "Not significant / depleted",
    ]
    counts = species["category"].value_counts().reindex(category_order).fillna(0)

    strata = []
    quality = pd.read_csv(args.quality, sep="\t").set_index("stratum")
    quality_labels = {
        "ALL": "All informative",
        "N50_ge10kb": "N50 ≥ 10 kb",
        "N50_ge100kb": "N50 ≥ 100 kb",
    }
    for label, path, key in [
        (quality_labels["ALL"], args.overall_all, "ALL"),
        (quality_labels["N50_ge10kb"], args.overall_n50_10kb, "N50_ge10kb"),
        (quality_labels["N50_ge100kb"], args.overall_n50_100kb, "N50_ge100kb"),
    ]:
        metrics = pd.read_csv(path, sep="\t").set_index("metric")["value"]
        if key == "N50_ge10kb":
            pairs = int(quality.loc[["N50_10kb_100kb", "N50_ge100kb"], "pairs"].sum())
            ge2 = int(quality.loc[["N50_10kb_100kb", "N50_ge100kb"], "pairs_ge2"].sum())
            qrow = {"pairs": pairs, "pairs_ge2": ge2}
        else:
            qrow = quality.loc[key]
        strata.append({
            "label": label,
            "fraction": float(qrow["pairs_ge2"] / qrow["pairs"]),
            "pairs": int(qrow["pairs"]),
            "species": int(metrics["species_tested"]),
        })

    fig = plt.figure(figsize=(180 / 25.4, 72 / 25.4))
    outer = fig.add_gridspec(
        1, 2, width_ratios=[1.18, 1.0], wspace=0.30,
        left=0.07, right=0.985, top=0.88, bottom=0.16,
    )

    ax_a = fig.add_subplot(outer[0])
    for category in category_order:
        sub = species[species["category"] == category]
        x = pd.to_numeric(sub["SNP_auc"], errors="coerce")
        y = pd.to_numeric(sub["BP_auc"], errors="coerce")
        ax_a.scatter(
            x, y, s=13 if category.startswith("Not") else 17,
            c=CATEGORY_COLORS[category],
            alpha=0.68 if category.startswith("Not") else 0.88,
            linewidths=0, label=f"{category} (n={int(counts[category])})",
            rasterized=True,
        )
    lim = (0.28, 0.82)
    ax_a.plot(lim, lim, color="0.45", lw=0.7, ls="--", zorder=1)
    ax_a.axvline(0.5, color="0.70", lw=0.6, zorder=0)
    ax_a.axhline(0.5, color="0.70", lw=0.6, zorder=0)
    ax_a.set_xlim(lim)
    ax_a.set_ylim(lim)
    ax_a.set_xlabel("AUC: oral-airway > within-site SNP distance")
    ax_a.set_ylabel("AUC: oral-airway > within-site breakpoints")
    ax_a.legend(loc="upper left", handletextpad=0.15, borderaxespad=0.2)
    ax_a.text(
        0.98, 0.03,
        "200 eligible species; BH q < 0.05\nand AUC > 0.5 defines enrichment",
        transform=ax_a.transAxes, ha="right", va="bottom",
        fontsize=5.5, color="0.25",
    )

    ax_b = fig.add_subplot(outer[1])
    xpos = np.arange(len(strata))
    fractions = [x["fraction"] for x in strata]
    bars = ax_b.bar(
        xpos, fractions, width=0.58,
        color=["#4D4D4D", "#2166AC", "#B2182B"], edgecolor="white", lw=0.3,
    )
    for bar, item in zip(bars, strata):
        ax_b.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012,
            f"{item['fraction'] * 100:.1f}%\nn={item['pairs']:,}",
            ha="center", va="bottom", fontsize=5.5,
        )
    ax_b.set_xticks(xpos, [x["label"] for x in strata])
    ax_b.set_ylim(0, 0.86)
    ax_b.set_ylabel("Pairs with ≥2 breakpoints")
    ax_b.set_xlabel("Minimum pair-level contig N50")
    ax_b.text(
        0.03, 0.88, "Contiguity sensitivity", transform=ax_b.transAxes,
        ha="left", va="top", fontsize=6, weight="bold",
    )

    fig.text(0.008, 0.96, "a", ha="left", va="top", fontsize=8, weight="bold")
    fig.text(0.605, 0.96, "b", ha="left", va="top", fontsize=8, weight="bold")
    fig.canvas.draw()
    save_figure(fig, args.outstem, args.dpi)

    qa = {
        "species_rows": len(species),
        "category_counts": {key: int(counts[key]) for key in category_order},
        "strata": strata,
        "random_seed_for_future_jitter": 42,
    }
    Path(f"{args.outstem}.qa.json").write_text(json.dumps(qa, indent=2) + "\n")


if __name__ == "__main__":
    main()
