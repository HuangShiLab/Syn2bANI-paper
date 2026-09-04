#!/usr/bin/env python3
"""Generate Figure 4: Computational efficiency on GTDB-R207 subsets.

Publication-quality 2x2 panel figure using the measured efficiency
benchmarks in results/efficiency_v8/.

Panels:
  (a) Wall time vs. number of pairs (log-log).
  (b) Peak memory vs. number of pairs (log-log).
  (c) Throughput (comparisons per second) vs. number of pairs (log-log).
  (d) Time breakdown at the largest subset (n=22): sketch vs. compare.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

set_publication_style()
D = Path("results/efficiency_v8")
OUT = Path("paper/figures/main/fig4_efficiency")


def load_runtime():
    """Load runtime_scaling.tsv and aggregate replicates by median."""
    df = pd.read_csv(D / "runtime_scaling.tsv", sep="\t")
    # Summarise across the 3 replicates per (tool, mode, n_genomes).
    summary = (
        df.groupby(["tool", "mode", "n_genomes", "n_pairs"])
        .agg(wall_s=("wall_s", "median"), peak_rss_mb=("peak_rss_mb", "median"))
        .reset_index()
    )
    summary["throughput"] = summary["n_pairs"] / summary["wall_s"]
    return summary


def load_sketch():
    """Load sketch_benchmark.tsv and aggregate replicates by median."""
    df = pd.read_csv(D / "sketch_benchmark.tsv", sep="\t")
    # Some replicates share the same artifact path (deleted .tmp files);
    # all reps for a (tool, n_genomes) used identical inputs, so median is valid.
    summary = (
        df.groupby(["tool", "n_genomes"])
        .agg(sketch_s=("wall_s", "median"), total_size_kb=("total_size_kb", "median"))
        .reset_index()
    )
    return summary


def panel_a(ax, rt):
    """Wall time versus number of pairs."""
    series = [
        ("syn2bani", "ani_fasta", "syn2bANI (FASTA)", COLORS["blue"], "o"),
        ("syn2bani", "ani_sketches", "syn2bANI (sketch reuse)", COLORS["sky_blue"], "^"),
        ("skani", "dist", "skani dist", COLORS["bluish_green"], "s"),
        ("fastani", "all_vs_all", "FastANI", COLORS["reddish_purple"], "D"),
    ]
    for tool, mode, label, color, marker in series:
        sub = rt[(rt["tool"] == tool) & (rt["mode"] == mode)].sort_values("n_pairs")
        ax.plot(
            sub["n_pairs"], sub["wall_s"],
            marker=marker, color=color, label=label,
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of pairs (all-vs-all, n\u00b2)")
    ax.set_ylabel("Wall time (s)")
    ax.set_xlim(3, 700)
    label_panel(ax, "a")


def panel_b(ax, rt):
    """Peak memory versus number of pairs."""
    series = [
        ("syn2bani", "ani_fasta", "syn2bANI (FASTA)", COLORS["blue"], "o"),
        ("syn2bani", "ani_sketches", "syn2bANI (sketch reuse)", COLORS["sky_blue"], "^"),
        ("skani", "dist", "skani dist", COLORS["bluish_green"], "s"),
        ("fastani", "all_vs_all", "FastANI", COLORS["reddish_purple"], "D"),
    ]
    for tool, mode, label, color, marker in series:
        sub = rt[(rt["tool"] == tool) & (rt["mode"] == mode)].sort_values("n_pairs")
        ax.plot(
            sub["n_pairs"], sub["peak_rss_mb"],
            marker=marker, color=color, label=label,
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of pairs (all-vs-all, n\u00b2)")
    ax.set_ylabel("Peak RSS (MiB)")
    ax.set_xlim(3, 700)
    label_panel(ax, "b")


def panel_c(ax, rt):
    """Throughput versus number of pairs."""
    series = [
        ("syn2bani", "ani_fasta", "syn2bANI (FASTA)", COLORS["blue"], "o"),
        ("syn2bani", "ani_sketches", "syn2bANI (sketch reuse)", COLORS["sky_blue"], "^"),
        ("skani", "dist", "skani dist", COLORS["bluish_green"], "s"),
        ("fastani", "all_vs_all", "FastANI", COLORS["reddish_purple"], "D"),
    ]
    for tool, mode, label, color, marker in series:
        sub = rt[(rt["tool"] == tool) & (rt["mode"] == mode)].sort_values("n_pairs")
        ax.plot(
            sub["n_pairs"], sub["throughput"],
            marker=marker, color=color, label=label,
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of pairs (all-vs-all, n\u00b2)")
    ax.set_ylabel("Throughput (comparisons / s)")
    ax.set_xlim(3, 700)
    label_panel(ax, "c")


def panel_d(ax, rt, sk):
    """Time breakdown at the largest subset (n=22).

    For sketch-reuse workflows the total cost is the one-off sketch plus the
    all-vs-all comparison.  FastANI and syn2bani FASTA have no separately
    timed sketch step, so they are shown as single compute bars.
    """
    n = 22
    # Comparison/dist times at n=22.
    cmp = rt[rt["n_genomes"] == n].set_index(["tool", "mode"])
    # Sketch times at n=22.
    skt = sk[sk["n_genomes"] == n].set_index("tool")

    bars = []
    # Workflow labels and measured components.
    items = [
        ("syn2bANI\nFASTA", [("compute", cmp.loc[("syn2bani", "ani_fasta"), "wall_s"])], COLORS["blue"]),
        ("syn2bANI\nsketch +\ncompare", [
            ("sketch", skt.loc["syn2bani", "sketch_s"]),
            ("compare", cmp.loc[("syn2bani", "ani_sketches"), "wall_s"]),
        ], COLORS["sky_blue"]),
        ("skani\nsketch +\ndist", [
            ("sketch", skt.loc["skani", "sketch_s"]),
            ("compare", cmp.loc[("skani", "dist"), "wall_s"]),
        ], COLORS["bluish_green"]),
        ("FastANI", [("compute", cmp.loc[("fastani", "all_vs_all"), "wall_s"])], COLORS["reddish_purple"]),
    ]

    x = np.arange(len(items))
    width = 0.55
    for i, (label, segments, base_color) in enumerate(items):
        bottom = 0.0
        total = sum(v for _, v in segments)
        for j, (seg_name, value) in enumerate(segments):
            # Alternate slightly within the same workflow for sketch/compare.
            color = base_color if len(segments) == 1 else (
                base_color if j == 0 else _lighten(base_color, 0.35)
            )
            ax.bar(i, value, width, bottom=bottom, color=color, edgecolor="white", linewidth=0.4)
            # Label segment if it occupies a reasonable fraction of the bar.
            if value / total > 0.12:
                ax.text(i, bottom + value / 2, seg_name, ha="center", va="center",
                        fontsize=6, color="white", fontweight="bold")
            bottom += value
        bars.append(bottom)
        # Annotate very short bars with their total value above the bar.
        if total < max(bars) * 0.08:
            ax.text(i, total + max(bars) * 0.02, f"{total:.2f}s", ha="center", va="bottom",
                    fontsize=6, color=COLORS["black"])

    ax.set_xticks(x)
    ax.set_xticklabels([label for label, _, _ in items], fontsize=7, linespacing=0.85)
    ax.set_ylabel("Wall time (s)")
    ax.set_ylim(0, max(bars) * 1.25)
    label_panel(ax, "d")

    # Manual legend for panel d.
    from matplotlib.patches import Patch
    legend_elems = [
        Patch(facecolor=COLORS["blue"], label="compute"),
        Patch(facecolor=COLORS["sky_blue"], label="sketch"),
        Patch(facecolor=_lighten(COLORS["sky_blue"], 0.35), label="compare"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", fontsize=7, frameon=False)


def _lighten(hexcolor, fraction):
    """Blend hexcolor with white by fraction (0 = original, 1 = white)."""
    hexcolor = hexcolor.lstrip("#")
    rgb = tuple(int(hexcolor[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    light = tuple(c + (1.0 - c) * fraction for c in rgb)
    return "#{:02x}{:02x}{:02x}".format(*(int(v * 255) for v in light))


def main():
    rt = load_runtime()
    sk = load_sketch()

    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.85))

    panel_a(axes[0, 0], rt)
    panel_b(axes[0, 1], rt)
    panel_c(axes[1, 0], rt)
    panel_d(axes[1, 1], rt, sk)

    # Shared legend for panels a-c, placed below the figure.
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.02),
               ncol=4, fontsize=7, frameon=False)

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
