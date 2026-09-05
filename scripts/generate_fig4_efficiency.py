#!/usr/bin/env python3
"""Generate Figure 4: Computational efficiency on GTDB-R207 subsets.

Publication-quality 2x2 panel figure combining ANI-only benchmarks
(all-vs-all and one-to-all) with SV-inclusive workflow benchmarks.
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
    summary = (
        df.groupby(["tool", "mode", "n_genomes", "n_pairs"])
        .agg(wall_s=("wall_s", "median"), peak_rss_mb=("peak_rss_mb", "median"))
        .reset_index()
    )
    summary["throughput"] = summary["n_pairs"] / summary["wall_s"]
    return summary


def load_one_to_all():
    """Load one_to_all_scaling.tsv and aggregate replicates by median."""
    df = pd.read_csv(D / "one_to_all_scaling.tsv", sep="\t")
    summary = (
        df.groupby(["tool", "mode", "n_genomes"])
        .agg(wall_s=("wall_s", "median"))
        .reset_index()
    )
    return summary


def load_sketch():
    """Load sketch_benchmark.tsv and aggregate replicates by median."""
    df = pd.read_csv(D / "sketch_benchmark.tsv", sep="\t")
    summary = (
        df.groupby(["tool", "n_genomes"])
        .agg(sketch_s=("wall_s", "median"), total_size_kb=("total_size_kb", "median"))
        .reset_index()
    )
    return summary


def load_sv():
    """Load SV benchmark and aggregate replicates by median."""
    df = pd.read_csv(D / "sv_benchmark.tsv", sep="\t")
    summary = (
        df.groupby(["mode", "n_genomes", "n_pairs"])
        .agg(skani_wall_s=("skani_wall_s", "median"),
             dnadiff_wall_s=("dnadiff_wall_s", "median"),
             n_ok=("n_ok", "median"))
        .reset_index()
    )
    summary["total_wall_s"] = summary["skani_wall_s"] + summary["dnadiff_wall_s"]
    return summary


def load_struct():
    """Load syn2bani struct benchmark and aggregate replicates by median."""
    df = pd.read_csv(D / "syn2b_struct_benchmark.tsv", sep="\t")
    summary = (
        df.groupby(["mode", "n_genomes", "n_pairs"])
        .agg(struct_wall_s=("struct_wall_s", "median"), n_ok=("n_ok", "median"))
        .reset_index()
    )
    return summary


def panel_a(ax, rt):
    """All-vs-all ANI wall time versus number of pairs."""
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


def panel_b(ax, oa):
    """One-to-all ANI wall time versus number of references."""
    series = [
        ("syn2bani", "one_to_all_fasta", "syn2bANI (FASTA)", COLORS["blue"], "o"),
        ("syn2bani", "one_to_all_sketches", "syn2bANI (sketch reuse)", COLORS["sky_blue"], "^"),
        ("skani", "one_to_all_dist", "skani dist", COLORS["bluish_green"], "s"),
        ("fastani", "one_to_all", "FastANI", COLORS["reddish_purple"], "D"),
    ]
    for tool, mode, label, color, marker in series:
        sub = oa[(oa["tool"] == tool) & (oa["mode"] == mode)].sort_values("n_genomes")
        ax.plot(
            sub["n_genomes"], sub["wall_s"],
            marker=marker, color=color, label=label,
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of reference genomes")
    ax.set_ylabel("Wall time (s)")
    ax.set_xlim(1.5, 30)
    label_panel(ax, "b")


def panel_c(ax, rt):
    """Peak memory versus number of pairs (all-vs-all)."""
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
    label_panel(ax, "c")


def panel_d(ax, sv, st, rt, sk):
    """SV-inclusive workflow wall time versus number of pairs.

    Workflows shown:
      - dnadiff alone (alignment-based SV truth)
      - skani + dnadiff (fast ANI pre-filter + alignment-based SV)
      - syn2bANI struct alone (SV from FASTA)
      - syn2bANI FASTA (ANI + SV from FASTA)
      - syn2bANI sketch (one-off sketch + ANI + SV)
    """
    dnadiff = sv[sv["mode"] == "dnadiff"].sort_values("n_pairs")
    skani_dna = sv[sv["mode"] == "skani_dnadiff"].sort_values("n_pairs")
    struct = st.sort_values("n_pairs").copy()
    struct["total_wall_s"] = struct["struct_wall_s"]

    # syn2bANI FASTA workflow = ani_fasta + struct
    ani_fasta = rt[(rt["tool"] == "syn2bani") & (rt["mode"] == "ani_fasta")][["n_genomes", "n_pairs", "wall_s"]]
    s2b_fasta = struct.merge(ani_fasta, on=["n_genomes", "n_pairs"], how="inner").copy()
    s2b_fasta["total_wall_s"] = s2b_fasta["struct_wall_s"] + s2b_fasta["wall_s"]

    # syn2bANI sketch workflow = sketch + ani_sketches + struct
    ani_sk = rt[(rt["tool"] == "syn2bani") & (rt["mode"] == "ani_sketches")][["n_genomes", "n_pairs", "wall_s"]]
    sk_syn = sk[sk["tool"] == "syn2bani"][["n_genomes", "sketch_s"]]
    s2b_sk = struct.merge(ani_sk, on=["n_genomes", "n_pairs"], how="inner").merge(sk_syn, on="n_genomes", how="inner").copy()
    s2b_sk["total_wall_s"] = s2b_sk["struct_wall_s"] + s2b_sk["wall_s"] + s2b_sk["sketch_s"]

    series = [
        (dnadiff, "dnadiff", COLORS["vermillion"], "o"),
        (skani_dna, "skani + dnadiff", COLORS["reddish_purple"], "s"),
        (struct, "syn2bANI struct", COLORS["bluish_green"], "^"),
        (s2b_fasta, "syn2bANI FASTA (ANI+SV)", COLORS["blue"], "D"),
        (s2b_sk, "syn2bANI sketch (ANI+SV)", COLORS["sky_blue"], "v"),
    ]
    for sub, label, color, marker in series:
        ax.plot(
            sub["n_pairs"], sub["total_wall_s"],
            marker=marker, color=color, label=label,
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of pairs (all-vs-all, n\u00b2)")
    ax.set_ylabel("Wall time (s)")
    ax.set_xlim(3, 700)
    label_panel(ax, "d")


def main():
    rt = load_runtime()
    oa = load_one_to_all()
    sk = load_sketch()
    sv = load_sv()
    st = load_struct()

    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.85))

    panel_a(axes[0, 0], rt)
    panel_b(axes[0, 1], oa)
    panel_c(axes[1, 0], rt)
    panel_d(axes[1, 1], sv, st, rt, sk)

    # Shared legend for panels a-c, placed below the figure.
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.27, 0.02),
               ncol=2, fontsize=7, frameon=False)

    # Separate legend for panel d (SV workflows), below right.
    handles_d, labels_d = axes[1, 1].get_legend_handles_labels()
    fig.legend(handles_d, labels_d, loc="upper center", bbox_to_anchor=(0.75, 0.02),
               ncol=2, fontsize=7, frameon=False)

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
