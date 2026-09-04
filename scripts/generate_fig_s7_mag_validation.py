#!/usr/bin/env python3
"""Generate Supplementary Figure S7: Accuracy on 695 CAMI2 MAGs.

Caption:
"695 CAMI2 bins vs dnadiff ANIm truth (raw gated syn2bani; skani; FastANI).
(a) Syn2bANI raw gated estimate vs dnadiff ANIm truth, colored by contamination
class. (b) Absolute-error distributions by tool. (c) Syn2bANI error by CheckM2
quality tier."
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_style


def load_data():
    """Load per-bin MAG validation data.

    Returns a DataFrame with one row per CAMI2 bin and columns for
    Syn2bANI (raw gated), skani, FastANI, dnadiff ANIm truth, CheckM2 tier,
    and contamination class.
    """
    repo = Path(__file__).parent.parent
    collect = repo / "results" / "mag_validation" / "collect"

    # Per-bin quality metadata (CheckM2 tier and contamination class)
    bins = pd.read_csv(collect / "bins.tsv", sep="\t")

    # dnadiff ANIm truth values
    truth = pd.read_csv(collect / "truth_dnadiff.tsv", sep="\t")

    # Tool estimates: keep only the anchor role to obtain one row per bin
    fast = pd.read_csv(collect / "ani_fast_tools.tsv", sep="\t")
    fast = fast[fast["role"] == "anchor"].copy()

    # Merge everything on bin id
    df = fast.merge(truth[["bin", "anim_ani"]], on="bin", how="left")
    df = df.merge(
        bins[["bin", "tier", "class", "checkm2_comp", "checkm2_cont"]],
        on="bin",
        how="left",
    )

    # Rename for clarity
    df = df.rename(
        columns={
            "s2b_ani_gated": "syn2bani",
            "skani_ani": "skani",
            "fastani_ani": "fastani",
            "anim_ani": "anim_truth",
        }
    )

    # Ensure ordering of categorical variables matches conventional labels
    df["tier"] = pd.Categorical(df["tier"], categories=["HQ", "MQ", "LQ"], ordered=True)
    df["class"] = pd.Categorical(
        df["class"],
        categories=["clean", "strain-mixed", "cross-species"],
        ordered=True,
    )

    # Absolute errors in percentage points (truth is already in percent)
    df["err_syn2bani"] = (df["syn2bani"] - df["anim_truth"]).abs()
    df["err_skani"] = (df["skani"] - df["anim_truth"]).abs()
    df["err_fastani"] = (df["fastani"] - df["anim_truth"]).abs()

    return df


def plot_fig_s7(df, outpath):
    plot_style.set_publication_style()

    # 1x3 layout, double-column width, moderate height
    fig, axes = plt.subplots(1, 3, figsize=plot_style.figure_size(17.5, aspect=0.42))
    fig.subplots_adjust(wspace=0.32)

    # ---- Panel (a): syn2bani vs ANIm truth, colored by contamination class ----
    ax = axes[0]
    class_colors = {
        "clean": plot_style.COLORS["blue"],
        "strain-mixed": plot_style.COLORS["orange"],
        "cross-species": plot_style.COLORS["vermillion"],
    }
    class_labels = {
        "clean": "clean",
        "strain-mixed": "strain-mixed",
        "cross-species": "cross-species",
    }

    # Plot in reverse order so rarer classes appear on top
    for cls in ["clean", "strain-mixed", "cross-species"]:
        sub = df[df["class"] == cls]
        ax.scatter(
            sub["anim_truth"],
            sub["syn2bani"],
            c=class_colors[cls],
            s=8,
            alpha=0.55,
            edgecolors="none",
            label=class_labels[cls],
            rasterized=True,
        )

    # Identity reference line
    lims = [94.5, 100.5]
    ax.plot(lims, lims, color=plot_style.COLORS["black"], linestyle="--", linewidth=0.8, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ANIm truth (dnadiff, %)")
    ax.set_ylabel("Syn2bANI raw gated estimate (%)")
    ax.legend(title="contamination class", loc="lower right", handletextpad=0.2)
    plot_style.label_panel(ax, "a")

    # ---- Panel (b): absolute error distributions by tool ----
    ax = axes[1]
    tool_data = [
        df["err_syn2bani"].dropna().values,
        df["err_skani"].dropna().values,
        df["err_fastani"].dropna().values,
    ]
    tool_labels = ["syn2bani", "skani", "FastANI"]
    tool_colors = [
        plot_style.COLORS["blue"],
        plot_style.COLORS["orange"],
        plot_style.COLORS["vermillion"],
    ]

    bp = ax.boxplot(
        tool_data,
        tick_labels=tool_labels,
        patch_artist=True,
        widths=0.55,
        medianprops={"color": plot_style.COLORS["black"], "linewidth": 1.0},
        whiskerprops={"color": plot_style.COLORS["black"], "linewidth": 0.8},
        capprops={"color": plot_style.COLORS["black"], "linewidth": 0.8},
        flierprops={
            "marker": "o",
            "markerfacecolor": plot_style.COLORS["grey"],
            "markeredgecolor": "none",
            "markersize": 2,
            "alpha": 0.4,
        },
    )
    for patch, color in zip(bp["boxes"], tool_colors):
        patch.set_facecolor(color)
        patch.set_edgecolor(plot_style.COLORS["black"])
        patch.set_linewidth(0.8)
        patch.set_alpha(0.75)

    ax.set_ylabel("|error| vs ANIm (pp)")
    ax.set_yscale("log")
    ax.set_ylim(8e-5, 2.0)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.minorticks_on()
    plot_style.label_panel(ax, "b")

    # ---- Panel (c): syn2bani error by CheckM2 quality tier ----
    ax = axes[2]
    tier_data = [
        df[df["tier"] == "HQ"]["err_syn2bani"].dropna().values,
        df[df["tier"] == "MQ"]["err_syn2bani"].dropna().values,
        df[df["tier"] == "LQ"]["err_syn2bani"].dropna().values,
    ]
    tier_labels = ["HQ", "MQ", "LQ"]
    tier_color = plot_style.COLORS["sky_blue"]

    bp = ax.boxplot(
        tier_data,
        tick_labels=tier_labels,
        patch_artist=True,
        widths=0.55,
        medianprops={"color": plot_style.COLORS["black"], "linewidth": 1.0},
        whiskerprops={"color": plot_style.COLORS["black"], "linewidth": 0.8},
        capprops={"color": plot_style.COLORS["black"], "linewidth": 0.8},
        flierprops={
            "marker": "o",
            "markerfacecolor": plot_style.COLORS["grey"],
            "markeredgecolor": "none",
            "markersize": 2,
            "alpha": 0.4,
        },
    )
    for patch in bp["boxes"]:
        patch.set_facecolor(tier_color)
        patch.set_edgecolor(plot_style.COLORS["black"])
        patch.set_linewidth(0.8)
        patch.set_alpha(0.75)

    ax.set_ylabel("|syn2bani error| (pp)")
    ax.set_yscale("log")
    ax.set_ylim(8e-5, 2.0)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.minorticks_on()
    plot_style.label_panel(ax, "c")

    # Save
    plot_style.save_figure(fig, outpath, formats=("png", "pdf"))
    plt.close(fig)


def main():
    df = load_data()
    n_total = len(df)
    n_fastani = df["fastani"].notna().sum()
    print(f"Loaded {n_total} CAMI2 MAGs ({n_fastani} with FastANI estimates)")

    repo = Path(__file__).parent.parent
    outdir = repo / "paper" / "figures" / "supplementary"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "fig_s7_mag_validation"

    plot_fig_s7(df, outpath)


if __name__ == "__main__":
    main()
