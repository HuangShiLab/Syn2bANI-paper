#!/usr/bin/env python3
"""Generate Figure 3: Unified GTDB-R207 80-100% benchmark against ANIm truth.

Publication-quality 2x2 panel figure.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS
import matplotlib.pyplot as plt

set_publication_style()
D = Path("results/gtdb50k")
OUT = Path("paper/figures/main/fig3_gtdb_r207_benchmark")

# Consistent method styling across panels.
METHODS = [
    ("Syn2bANI raw", "syn2bani_raw", "ani_gated", COLORS['sky_blue']),
    ("Syn2bANI v6", "syn2bani_v6", "ani_cal", COLORS['orange']),
    ("Syn2bANI hybrid", "syn2bani_hybrid", "hybrid", COLORS['reddish_purple']),
    ("skani", "skani", "skani_ani", COLORS['bluish_green']),
    ("FastANI", "fastani", "fastani_ani", COLORS['vermillion']),
]

BANDS = ["80-85", "85-90", "90-95", "95-97", "97-100"]


def load_data():
    """Load and merge main 43,334 held-out pairs with 727 high-ANI test pairs."""
    s2b = pd.read_csv(D / "s2b_50k.tsv", sep="\t")
    pairs = pd.read_csv(D / "pairs_50k.tsv", sep="\t")
    fastani = pd.read_csv(D / "fastani_50k.tsv", sep="\t")
    truth = pd.read_csv(D / "truth_50k.tsv", sep="\t")
    ha = pd.read_csv(D / "high_ani_results.tsv", sep="\t")

    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]

    def make_hybrid(row):
        raw = row["ani"]
        cal = row["ani_cal"]
        if pd.isna(raw) or pd.isna(cal):
            return np.nan
        return raw if raw >= 98 else cal

    def unify_band(row):
        band = row["band"]
        anim = row["anim_ani"]
        if band == "95-100":
            return "95-97" if anim < 97 else "97-100"
        return band

    # Main benchmark pairs.
    df_main = truth[["pairid", "anim_ani"]].merge(
        s2b[["pairid", "ani", "ani_cal", "ani_gated"]], on="pairid", how="left"
    )
    df_main = df_main.merge(
        pairs[["pairid", "skani_ani", "band"]], on="pairid", how="left"
    )
    df_main = df_main.merge(
        fastani[["pairid", "fastani_ani"]], on="pairid", how="left"
    )
    df_main["hybrid"] = df_main.apply(make_hybrid, axis=1)
    df_main["band"] = df_main.apply(unify_band, axis=1)

    # High-ANI held-out test pairs only.
    test = ha[ha["split"] == "test"].copy()
    test["hybrid"] = test.apply(make_hybrid, axis=1)

    cols = ["pairid", "anim_ani", "ani", "ani_cal", "ani_gated",
            "skani_ani", "fastani_ani", "band", "hybrid"]
    df = pd.concat([df_main[cols], test[cols]], ignore_index=True)
    return df


def compute_mae(df, col):
    """Return overall MAE for a method column."""
    sub = df.dropna(subset=[col, "anim_ani"])
    return (sub[col] - sub["anim_ani"]).abs().mean()


def band_mae(df, col, band):
    """Return MAE for a method column within a band."""
    sub = df[(df["band"] == band)].dropna(subset=[col, "anim_ani"])
    if len(sub) == 0:
        return np.nan
    return (sub[col] - sub["anim_ani"]).abs().mean()


def panel_scatter_all(ax, df):
    """Panel (a): estimated vs ANIm truth across all pairs."""
    for label, key, col, color in METHODS:
        sub = df.dropna(subset=[col, "anim_ani"])
        mae = compute_mae(df, col)
        ax.scatter(
            sub["anim_ani"], sub[col],
            s=3, alpha=0.35, c=color, edgecolors="none",
            rasterized=True,
        )
        ax.scatter([], [], c=color, s=20, label=f"{label} (MAE {mae:.2f})")
    ax.plot([80, 100], [80, 100], "k--", lw=0.8)
    ax.set_xlim(80, 100)
    ax.set_ylim(80, 100)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ANIm truth (%)")
    ax.set_ylabel("Estimated ANI (%)")
    ax.set_title("Unified GTDB-R207 benchmark", fontsize=8)
    ax.legend(loc="lower right", fontsize=6.5, handletextpad=0.2)
    label_panel(ax, "a")


def panel_band_mae(ax, df):
    """Panel (b): per-band MAE bar chart with value labels."""
    x = np.arange(len(BANDS))
    width = 0.15
    for i, (label, key, col, color) in enumerate(METHODS):
        values = [band_mae(df, col, b) for b in BANDS]
        bars = ax.bar(x + (i - 2) * width, values, width, label=label, color=color)
        # Annotate each bar with its MAE (small font).
        for bar, val in zip(bars, values):
            if not np.isnan(val):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.03,
                    f"{val:.2f}",
                    ha="center", va="bottom", fontsize=4.5, rotation=90,
                    color="#333333"
                )

    ax.set_xticks(x)
    ax.set_xticklabels(BANDS)
    ax.set_xlabel("ANI band")
    ax.set_ylabel("MAE (ANI points)")
    ax.set_title("MAE by band", fontsize=8)
    ax.set_ylim(0, 2.2)
    ax.legend(loc="upper right", fontsize=6.5, handletextpad=0.2)
    label_panel(ax, "b")


def panel_error_distribution(ax, df):
    """Panel (c): signed-error distributions as histograms."""
    bins = np.linspace(-15, 15, 61)
    for label, key, col, color in METHODS:
        sub = df.dropna(subset=[col, "anim_ani"])
        errors = sub[col] - sub["anim_ani"]
        ax.hist(errors, bins=bins, alpha=0.5, color=color, label=label,
                histtype="stepfilled", edgecolor="none")
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("Error (estimate - truth, ANI points)")
    ax.set_ylabel("Count")
    ax.set_title("Signed-error distribution", fontsize=8)
    ax.set_xlim(-15, 15)
    ax.legend(loc="upper right", fontsize=7, handletextpad=0.2)
    label_panel(ax, "c")


def panel_high_ani_zoom(ax, df):
    """Panel (d): high-ANI zoom (95-100%)."""
    sub = df[df["band"].isin(["95-97", "97-100"])]
    for label, key, col, color in METHODS:
        s = sub.dropna(subset=[col, "anim_ani"])
        mae = (s[col] - s["anim_ani"]).abs().mean()
        ax.scatter(
            s["anim_ani"], s[col],
            s=4, alpha=0.5, c=color, edgecolors="none",
            rasterized=True,
        )
        ax.scatter([], [], c=color, s=20, label=f"{label} (MAE {mae:.2f})")
    ax.plot([95, 100], [95, 100], "k--", lw=0.8)
    ax.set_xlim(95, 100)
    ax.set_ylim(95, 100)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ANIm truth (%)")
    ax.set_ylabel("Estimated ANI (%)")
    ax.set_title("High-ANI regime (95-100%)", fontsize=8)
    ax.legend(loc="lower right", fontsize=6.5, handletextpad=0.2)
    label_panel(ax, "d")


def main():
    df = load_data()
    print(f"Loaded {len(df)} pairs")
    print("Band counts:")
    for band in BANDS:
        print(f"  {band}: {(df['band'] == band).sum()}")
    print("\nOverall MAE:")
    for label, key, col, color in METHODS:
        mae = compute_mae(df, col)
        n = df[[col, "anim_ani"]].dropna().shape[0]
        print(f"  {label}: n={n}, MAE={mae:.4f}")

    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.95))

    panel_scatter_all(axes[0, 0], df)
    panel_band_mae(axes[0, 1], df)
    panel_error_distribution(axes[1, 0], df)
    panel_high_ani_zoom(axes[1, 1], df)

    plt.tight_layout()
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
