#!/usr/bin/env python3
"""Generate Figure 2: Accuracy and robustness under exact truth.

Reads raw simulation outputs and produces a publication-quality composite
figure with five panels:
  (a) estimated vs true ANI for Syn2bANI, skani, and FastANI;
  (b) signed error vs true ANI for the same three tools;
  (c) indel sweep at fixed 95% ANI (0–4 deletions per 100 kb);
  (d) simulated fragmentation (20–200 contigs) at 95% ANI;
  (e) accessory-content sweep at 95% core ANI, with ANI error and aligned
      fraction on twin axes.

Data provenance
---------------
The exact-truth simulations were generated in the companion Syn2bANI
repository under prototype/ and then copied/rerun into this paper repo.
This script looks for the needed TSVs in the following order:
  1. results/simulations/           (preferred, self-contained)
  2. results/sv_validation/         (copies already in the paper repo)
  3. ../Syn2bANI/prototype/         (original raw outputs)

The panel (a,b) cross-tool comparison requires simindel_cross_tool_4e.tsv,
which is currently only available in the Syn2bANI prototype directory; the
other four families are present in results/sv_validation/.

Outputs
-------
paper/figures/main/fig2_simulations.png
paper/figures/main/fig2_simulations.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]  # repository root
OUTDIR = ROOT / "paper" / "figures" / "main"

# Search roots in priority order for the raw simulation TSVs.
SEARCH_ROOTS = [
    ROOT / "results" / "simulations",
    ROOT / "results" / "sv_validation",
    Path("/Users/macstudio/Downloads/Syn2bANI/prototype"),
]


def find_file(name):
    """Return the first existing path for *name* across SEARCH_ROOTS."""
    for root in SEARCH_ROOTS:
        candidate = root / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Could not find {name} in any of: "
        + ", ".join(str(r) for r in SEARCH_ROOTS)
    )


# Data files
ANI_LADDER_TSV = "simindel_results_4e.tsv"
ANI_LADDER_MANIFEST_TSV = "simindel/manifest.tsv"
CROSS_TOOL_TSV = "simindel_cross_tool_4e.tsv"
INDEL_SWEEP_TSV = "simindel_sweep_results_4e.tsv"
INDEL_SWEEP_MANIFEST_TSV = "simindel_sweep/manifest.tsv"
FRAGMENTATION_TSV = "simfrag_fixed.tsv"  # includes fwd + shuffled replicates
ACCESSORY_TSV = "simacc_fixed.tsv"

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------
C_SYN2B = COLORS["vermillion"]       # Syn2bANI primary estimate
C_SKANI = COLORS["bluish_green"]     # skani
C_FASTANI = COLORS["reddish_purple"] # FastANI
C_TRUTH = COLORS["black"]            # identity / truth lines
C_GREY = COLORS["grey"]


def parse_ani_ladder():
    """Return DataFrame with true ANI and estimates for Syn2bANI/skani/FastANI."""
    res = pd.read_csv(find_file(ANI_LADDER_TSV), sep="\t")
    man = pd.read_csv(find_file(ANI_LADDER_MANIFEST_TSV), sep="\t")
    df = res.merge(man[["name", "true_ani"]], left_on="query", right_on="name", how="left")
    df["true_ani_pct"] = df["true_ani"] * 100.0
    df = df.sort_values("true_ani_pct").reset_index(drop=True)

    cross = pd.read_csv(find_file(CROSS_TOOL_TSV), sep="\t")
    cross = cross.merge(man[["name", "true_ani"]], on="name", how="left")
    cross["true_ani_pct"] = cross["true_ani"] * 100.0
    cross = cross.sort_values("true_ani_pct").reset_index(drop=True)

    df = df.merge(
        cross[["name", "skani_ani", "fastani_ani"]],
        on="name",
        how="left",
    )
    return df


def parse_indel_sweep():
    """Return indel-sweep DataFrame with indel rate and true ANI."""
    res = pd.read_csv(find_file(INDEL_SWEEP_TSV), sep="\t")
    man = pd.read_csv(find_file(INDEL_SWEEP_MANIFEST_TSV), sep="\t")
    df = res.merge(man[["name", "indel_rate", "true_ani"]],
                   left_on="query", right_on="name", how="left")
    df["true_ani_pct"] = df["true_ani"] * 100.0
    return df.sort_values("indel_rate").reset_index(drop=True)


def parse_fragmentation():
    """Return fragmentation DataFrame with contig counts."""
    df = pd.read_csv(find_file(FRAGMENTATION_TSV), sep="\t")
    df["n_contigs"] = df["query"].str.extract(r"q95_c(\d+)_")[0].astype(int)
    # Keep one replicate per contig count for plotting (the fixed TSV contains
    # both shuffled and forward-only assemblies; use the mean per count).
    grouped = df.groupby("n_contigs").agg({"ani": "mean", "af_query": "mean"}).reset_index()
    grouped = grouped.sort_values("n_contigs").reset_index(drop=True)
    return grouped


def parse_accessory():
    """Return accessory-sweep DataFrame with shuffled fraction."""
    df = pd.read_csv(find_file(ACCESSORY_TSV), sep="\t")
    df["accessory_frac"] = df["query"].str.extract(r"acc(\d+\.\d+)")[0].astype(float) * 100.0
    return df.sort_values("accessory_frac").reset_index(drop=True)


def mae(pred, truth):
    return np.mean(np.abs(np.asarray(pred) - np.asarray(truth)))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_panel_a(ax, df):
    """Estimated vs true ANI."""
    lo, hi = 84.0, 101.0
    ax.plot([lo, hi], [lo, hi], color=C_TRUTH, lw=0.8, ls="--", zorder=1)

    ax.scatter(df["true_ani_pct"], df["ani"], color=C_SYN2B, marker="o",
               s=20, label="Syn2bANI", zorder=3)
    ax.scatter(df["true_ani_pct"], df["skani_ani"], color=C_SKANI, marker="s",
               s=18, label="skani", zorder=3)
    ax.scatter(df["true_ani_pct"], df["fastani_ani"], color=C_FASTANI,
               marker="D", s=16, label="FastANI", zorder=3)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("True ANI (%)")
    ax.set_ylabel("Estimated ANI (%)")

    err_s2b = df["ani"] - df["true_ani_pct"]
    err_skani = df["skani_ani"] - df["true_ani_pct"]
    err_fastani = df["fastani_ani"] - df["true_ani_pct"]
    text = (
        f"MAE Syn2bANI {mae(df['ani'], df['true_ani_pct']):.3f}\n"
        f"MAE skani {mae(df['skani_ani'], df['true_ani_pct']):.3f}\n"
        f"MAE FastANI {mae(df['fastani_ani'], df['true_ani_pct']):.3f}"
    )
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=7, linespacing=1.3)

    ax.legend(loc="lower right", frameon=False)
    label_panel(ax, "a")


def plot_panel_b(ax, df):
    """Signed error vs true ANI."""
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)

    ax.scatter(df["true_ani_pct"], df["ani"] - df["true_ani_pct"],
               color=C_SYN2B, marker="o", s=20, label="Syn2bANI", zorder=3)
    ax.scatter(df["true_ani_pct"], df["skani_ani"] - df["true_ani_pct"],
               color=C_SKANI, marker="s", s=18, label="skani", zorder=3)
    ax.scatter(df["true_ani_pct"], df["fastani_ani"] - df["true_ani_pct"],
               color=C_FASTANI, marker="D", s=16, label="FastANI", zorder=3)

    ax.set_xlim(84.0, 101.0)
    ax.set_xlabel("True ANI (%)")
    ax.set_ylabel("Error (est − truth, ANI points)")

    ax.legend(loc="lower right", frameon=False)
    label_panel(ax, "b")


def plot_panel_c(ax, df):
    """Indel sweep at 95% ANI."""
    truth = df["true_ani_pct"].iloc[0]
    err = df["ani"] - truth

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    ax.plot(df["indel_rate"], err, color=C_SYN2B, marker="o", ms=5,
            label="Syn2bANI")

    ax.set_xlabel("Deletions per 100 kb")
    ax.set_ylabel("Error (est − 95.000, ANI points)")
    ax.set_xticks([0.0, 0.5, 1.0, 2.0, 4.0])
    ax.set_xticklabels(["0", "0.5", "1", "2", "4"])

    mae_val = mae(df["ani"], truth)
    ax.text(0.03, 0.97, f"MAE {mae_val:.3f}", transform=ax.transAxes,
            va="top", ha="left", fontsize=8)
    label_panel(ax, "c")


def plot_panel_d(ax, df):
    """Fragmentation at 95% ANI."""
    truth = 95.0
    err = df["ani"] - truth

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    ax.plot(df["n_contigs"], err, color=C_SYN2B, marker="o", ms=5)

    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200])
    ax.set_xticklabels([20, 50, 100, 200])
    ax.minorticks_off()
    ax.set_xlabel("Number of contigs")
    ax.set_ylabel("Error (est − 95.000, ANI points)")

    mae_val = mae(df["ani"], truth)
    ax.text(0.03, 0.97, f"MAE {mae_val:.3f}", transform=ax.transAxes,
            va="top", ha="left", fontsize=8)
    label_panel(ax, "d")


def plot_panel_e(ax, df):
    """Accessory-content sweep at 95% core ANI."""
    truth = 95.0
    frac = df["accessory_frac"]
    err = df["ani"] - truth

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    ax.plot(frac, err, color=C_SYN2B, marker="o", ms=5,
            label="ANI error")
    ax.set_xlabel("Accessory fraction (%)")
    ax.set_ylabel("Error (ANI points)", color=C_SYN2B)
    ax.tick_params(axis="y", labelcolor=C_SYN2B)
    ax.set_ylim(-0.3, 0.35)

    ax2 = ax.twinx()
    ax2.plot(frac, df["af_query"], color=C_GREY, marker="s", ms=4,
             ls=":", label="af_query")
    ax2.plot(frac, 1.0 - frac / 100.0, color=C_TRUTH, lw=0.8, ls="--",
             label="1 − F (truth)")
    ax2.set_ylabel("Aligned fraction (af_query)", color=C_GREY)
    ax2.tick_params(axis="y", labelcolor=C_GREY)
    ax2.set_ylim(0.4, 1.05)

    mae_val = mae(df["ani"], truth)
    ax.text(0.03, 0.97, f"MAE {mae_val:.3f}", transform=ax.transAxes,
            va="top", ha="left", fontsize=8)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=False, fontsize=6.5)
    label_panel(ax, "e")


def main():
    set_publication_style()

    # Load data
    df_ladder = parse_ani_ladder()
    df_indel = parse_indel_sweep()
    df_frag = parse_fragmentation()
    df_acc = parse_accessory()

    # Figure: 2 rows x 3 columns, leave bottom-right empty
    fig = plt.figure(figsize=figure_size(17.8, aspect=0.62))
    gs = fig.add_gridspec(2, 3, wspace=0.38, hspace=0.42)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    # gs[1, 2] intentionally left blank

    plot_panel_a(ax_a, df_ladder)
    plot_panel_b(ax_b, df_ladder)
    plot_panel_c(ax_c, df_indel)
    plot_panel_d(ax_d, df_frag)
    plot_panel_e(ax_e, df_acc)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    save_figure(fig, OUTDIR / "fig2_simulations")
    plt.close(fig)
    print("Figure 2 generated successfully.")


if __name__ == "__main__":
    main()
