#!/usr/bin/env python3
"""Generate Supplementary Figure S3: 2,074-pair GTDB-R207 band-holdout CV.

Bar chart of per-band MAE for raw gamma, calibrated v5, skani, and FastANI.
FastANI is available only for a 363-pair subset of the 2,074 pairs.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

set_publication_style()

D = Path("results/panel_by_band")
OUT = Path("paper/figures/supplementary/fig_s3_anim_by_band")

BANDS = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]
BAND_LABELS = ["80–85", "85–90", "90–95", "95–99"]

METHODS = [
    ("Syn2bANI raw gamma", "s2b_ani", COLORS["sky_blue"]),
    ("Syn2bANI cal. v5", "ridge_pred", COLORS["orange"]),
    ("skani", "skani_ani", COLORS["bluish_green"]),
    ("FastANI (363-pair subset)", "FastANI_subset", COLORS["vermillion"]),
]


def load_data():
    ep = pd.read_csv(D / "eval_pairs.tsv", sep="\t")
    ep["pairid"] = ep["query_asm"] + "__" + ep["ref_asm"]

    # Band-holdout ridge predictions for the deployed v5 model.
    ridge = pd.read_csv(D / "ridge_cv_preds_v5.tsv", sep="\t")
    ridge["pairid"] = ridge["query"] + "__" + ridge["reference"]
    ep = ep.merge(ridge[["pairid", "ridge_pred"]], on="pairid", how="left")

    # Pre-computed FastANI subset metrics.
    fa = pd.read_csv(D / "anim_main_table.tsv", sep="\t")
    fa = fa[fa["method"] == "FastANI_subset"].set_index("band")
    return ep, fa


def band_mae(df, col, band):
    sub = df[(df["band"] == band)].dropna(subset=[col, "anim_ani"])
    if len(sub) == 0:
        return np.nan
    return (sub[col] - sub["anim_ani"]).abs().mean()


def main():
    ep, fa = load_data()
    print(f"Loaded {len(ep)} eval pairs; {ep['ridge_pred'].notna().sum()} with v5 CV prediction")

    fig, ax = plt.subplots(figsize=figure_size(12.0, aspect=0.65))

    x = np.arange(len(BANDS))
    width = 0.2

    for i, (label, col, color) in enumerate(METHODS):
        if col == "FastANI_subset":
            vals = [fa.loc[b, "MAE"] if b in fa.index else np.nan for b in BANDS]
        else:
            vals = [band_mae(ep, col, b) for b in BANDS]
        bars = ax.bar(x + (i - 1.5) * width, vals, width, label=label, color=color)
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.05,
                    f"{val:.2f}",
                    ha="center", va="bottom", fontsize=6, rotation=90,
                    color="#333333"
                )
        # overall MAE annotation in legend text handled below

    ax.set_xticks(x)
    ax.set_xticklabels(BAND_LABELS)
    ax.set_xlabel("ANIm band (%)")
    ax.set_ylabel("MAE (ANI points)")
    ax.set_ylim(0, 2.7)
    ax.legend(loc="upper right", fontsize=7, handletextpad=0.2)

    # Overall MAE text for non-FastANI methods.
    overall_texts = []
    for label, col, _ in METHODS[:3]:
        sub = ep.dropna(subset=[col, "anim_ani"])
        mae = (sub[col] - sub["anim_ani"]).abs().mean()
        overall_texts.append(f"{label}: MAE {mae:.3f} (n={len(sub)})")
    fa_all = fa.loc["all", "MAE"]
    overall_texts.append(f"FastANI (subset): MAE {fa_all:.3f} (n={int(fa.loc['all', 'n'])})")
    ax.text(
        0.02, 0.98, "\n".join(overall_texts),
        transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9)
    )

    fig.suptitle(
        "Band-holdout cross-validation on the 2,074-pair GTDB-R207 ANIm benchmark",
        fontsize=9, y=0.98
    )

    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
