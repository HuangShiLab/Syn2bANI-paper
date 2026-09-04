#!/usr/bin/env python3
"""Generate Supplementary Figure S1: Synteny benchmark inversion ladder.

E. coli MG1655 evolved to ANI 95.00/98.00 (counted substitutions) with 0-32
non-overlapping inversions (100-400 kb; 2 true breakpoints each).
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import plot_style


def load_data():
    results = Path(__file__).parent.parent / "results" / "synteny_bench" / "synteny_ladder_results.tsv"
    return pd.read_csv(results, sep="\t")


def plot_s1(df, outpath):
    plot_style.set_publication_style()

    # Use a wide 1x3 layout (double-column width); aspect tuned for labels + legend
    fig, axes = plt.subplots(1, 3, figsize=plot_style.figure_size(17.5, aspect=0.42))

    colors = {
        0.95: plot_style.COLORS["blue"],
        0.98: plot_style.COLORS["vermillion"],
    }
    labels = {
        0.95: "ANI 95.00",
        0.98: "ANI 98.00",
    }

    inv_counts = sorted(df["n_inv"].unique())

    # Panel (a): breakpoint_count vs true breakpoints (identity line)
    ax = axes[0]
    for ani, sub in df.groupby("ani"):
        ax.plot(sub["true_breakpoints"], sub["breakpoint_count"],
                marker="o", linestyle="-", color=colors[ani],
                label=labels[ani], clip_on=False)
    lims = [min(df["true_breakpoints"].min(), df["breakpoint_count"].min()) - 2,
            max(df["true_breakpoints"].max(), df["breakpoint_count"].max()) + 2]
    ax.plot(lims, lims, color=plot_style.COLORS["grey"], linestyle="--", linewidth=1.0, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("True breakpoints")
    ax.set_ylabel("syn2bani breakpoint count")
    ax.set_xticks(np.arange(0, df["true_breakpoints"].max() + 1, 16))
    ax.set_yticks(np.arange(0, df["breakpoint_count"].max() + 1, 16))
    ax.legend(loc="upper left", frameon=False)
    plot_style.label_panel(ax, "a")

    # Panel (b): anchor_adjacency vs inversion count
    ax = axes[1]
    for ani, sub in df.groupby("ani"):
        ax.plot(sub["n_inv"], sub["anchor_adjacency"],
                marker="o", linestyle="-", color=colors[ani],
                label=labels[ani], clip_on=False)
    ax.set_xlabel("Inversions")
    ax.set_ylabel("Anchor adjacency")
    ax.set_xscale("symlog", linthresh=1, base=2)
    ax.set_xticks(inv_counts)
    ax.set_xticklabels([str(i) for i in inv_counts])
    ax.set_xlim(-0.5, 40)
    ax.set_ylim(0.975, 1.005)
    ax.set_yticks([0.98, 0.99, 1.00])
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    plot_style.label_panel(ax, "b")

    # Panel (c): ANI estimate vs inversion count
    ax = axes[2]
    for ani, sub in df.groupby("ani"):
        ax.plot(sub["n_inv"], sub["ani_est"],
                marker="o", linestyle="-", color=colors[ani],
                label=labels[ani], clip_on=False)
        ax.axhline(y=ani * 100, color=colors[ani], linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlabel("Inversions")
    ax.set_ylabel("ANI estimate (%)")
    ax.set_xscale("symlog", linthresh=1, base=2)
    ax.set_xticks(inv_counts)
    ax.set_xticklabels([str(i) for i in inv_counts])
    ax.set_xlim(-0.5, 40)
    ax.set_ylim(94.85, 98.15)
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    plot_style.label_panel(ax, "c")

    fig.tight_layout()
    plot_style.save_figure(fig, str(outpath.with_suffix("")), formats=("png", "pdf"))
    plt.close(fig)


def main():
    df = load_data()
    outdir = Path(__file__).parent.parent / "paper" / "figures" / "supplementary"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "fig_s1_synteny_ladder"
    plot_s1(df, outpath)


if __name__ == "__main__":
    main()
