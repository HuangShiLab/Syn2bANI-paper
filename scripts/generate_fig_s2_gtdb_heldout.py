#!/usr/bin/env python3
"""Generate Supplementary Figure S2: GTDB-R207 43,334 held-out pair benchmark.

Panels:
  (a-c) Calibrated Syn2bANI v5, skani, and FastANI vs ANIm truth.
  (d) Per-band MAE.
  (e) Signed-error distributions.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

set_publication_style()

D = Path("results/gtdb50k")
OUT = Path("paper/figures/supplementary/fig_s2_gtdb50k_heldout")

METHODS = [
    ("Syn2bANI cal. (v5)", "ani_cal", COLORS["orange"]),
    ("skani", "skani_ani", COLORS["bluish_green"]),
    ("FastANI", "fastani_ani", COLORS["vermillion"]),
]
BANDS = ["80-85", "85-90", "90-95", "95-100"]


def load_data():
    truth = pd.read_csv(D / "truth_50k.tsv", sep="\t")
    s2b = pd.read_csv(D / "s2b_50k.tsv", sep="\t")
    pairs = pd.read_csv(D / "pairs_50k.tsv", sep="\t")
    fastani = pd.read_csv(D / "fastani_50k.tsv", sep="\t")

    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]

    df = truth[["pairid", "anim_ani"]].merge(
        s2b[["pairid", "ani_cal"]], on="pairid", how="left"
    )
    df = df.merge(pairs[["pairid", "skani_ani", "band"]], on="pairid", how="left")
    df = df.merge(fastani[["pairid", "fastani_ani"]], on="pairid", how="left")
    return df


def mae(df, col):
    sub = df.dropna(subset=[col, "anim_ani"])
    return (sub[col] - sub["anim_ani"]).abs().mean()


def common_mae(df, cols):
    sub = df.dropna(subset=[*cols, "anim_ani"])
    return {c: (sub[c] - sub["anim_ani"]).abs().mean() for c in cols}


def panel_hexbin(ax, df, col, color, title, mae_val):
    sub = df.dropna(subset=[col, "anim_ani"])
    hb = ax.hexbin(
        sub["anim_ani"], sub[col], gridsize=60,
        cmap="viridis", bins="log", mincnt=1,
        rasterized=True, linewidths=0,
    )
    ax.plot([80, 100], [80, 100], "r--", lw=0.8, alpha=0.7)
    ax.set_xlim(80, 100)
    ax.set_ylim(80, 100)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ANIm truth (%)")
    ax.set_ylabel("Estimated ANI (%)")
    ax.set_title(f"{title}  MAE {mae_val:.3f}", fontsize=8)
    cbar = plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("pairs (log)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    return hb


def panel_band_mae(ax, df):
    common = df.dropna(subset=["ani_cal", "skani_ani", "fastani_ani", "anim_ani"])
    x = np.arange(len(BANDS))
    width = 0.25
    for i, (label, col, color) in enumerate(METHODS):
        vals = []
        for band in BANDS:
            sub = common[common["band"] == band]
            vals.append((sub[col] - sub["anim_ani"]).abs().mean())
        bars = ax.bar(x + (i - 1) * width, vals, width, label=label, color=color)
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.03,
                    f"{val:.2f}",
                    ha="center", va="bottom", fontsize=5.5, rotation=90,
                    color="#333333"
                )
    ax.set_xticks(x)
    ax.set_xticklabels(BANDS)
    ax.set_xlabel("ANI band (%)")
    ax.set_ylabel("MAE (ANI points)")
    ax.set_ylim(0, 2.3)
    ax.legend(loc="upper right", fontsize=7, handletextpad=0.2)


def panel_error_distribution(ax, df):
    common = df.dropna(subset=["ani_cal", "skani_ani", "fastani_ani", "anim_ani"])
    xs = np.linspace(-5, 5, 300)
    for label, col, color in METHODS:
        errors = common[col] - common["anim_ani"]
        bias = errors.mean()
        kde = gaussian_kde(errors)
        kde.set_bandwidth(0.25)
        ax.plot(xs, kde(xs), color=color, lw=1.5, label=f"{label} (bias {bias:+.2f})")
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("Error (estimate − ANIm, ANI points)")
    ax.set_ylabel("Density")
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.legend(loc="upper right", fontsize=7, handletextpad=0.2)


def main():
    df = load_data()
    common = df.dropna(subset=["ani_cal", "skani_ani", "fastani_ani", "anim_ani"])
    print(f"Loaded {len(df)} pairs; common subset n = {len(common)}")
    for label, col, _ in METHODS:
        print(f"  {label}: overall n = {df[[col, 'anim_ani']].dropna().shape[0]}, "
              f"common-subset MAE = {mae(common, col):.4f}")

    fig = plt.figure(figsize=figure_size(17.8, aspect=0.78))
    gs = GridSpec(2, 3, figure=fig, wspace=0.45, hspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1:])

    common_maes = common_mae(df, ["ani_cal", "skani_ani", "fastani_ani"])
    panel_hexbin(ax_a, df, "ani_cal", METHODS[0][2], "Syn2bANI cal. (v5)", common_maes["ani_cal"])
    panel_hexbin(ax_b, df, "skani_ani", METHODS[1][2], "skani", common_maes["skani_ani"])
    panel_hexbin(ax_c, df, "fastani_ani", METHODS[2][2], "FastANI", common_maes["fastani_ani"])

    panel_band_mae(ax_d, df)
    panel_error_distribution(ax_e, df)

    label_panel(ax_a, "a")
    label_panel(ax_b, "b")
    label_panel(ax_c, "c")
    label_panel(ax_d, "d")
    label_panel(ax_e, "e")

    fig.suptitle(
        "GTDB-R207 held-out benchmark: 43,334 same-genus pairs (scored subset n = 39,903)",
        fontsize=9, y=0.98
    )

    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
