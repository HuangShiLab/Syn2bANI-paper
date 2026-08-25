#!/usr/bin/env python3
"""Generate supplementary figures for Task 0b simulation families.

Reads rerun *_4e.tsv files from /Users/macstudio/Downloads/Syn2bANI/prototype and
writes publication-quality PNGs (and PDFs) to figures/report.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAPER = Path("/Users/macstudio/Downloads/Syn2bANI-paper")
S2B = Path("/Users/macstudio/Downloads/Syn2bANI")
FIGDIR = PAPER / "figures" / "report"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Data files from the rerun
SIMINDEL = S2B / "prototype" / "simindel_results_4e.tsv"
SIMINDEL_MANIFEST = S2B / "prototype" / "simindel" / "manifest.tsv"
SIMINDEL_SWEEP = S2B / "prototype" / "simindel_sweep_results_4e.tsv"
SIMINDEL_SWEEP_MANIFEST = S2B / "prototype" / "simindel_sweep" / "manifest.tsv"
SIMFRAG = S2B / "prototype" / "simfrag_results_4e.tsv"
SIMACC = S2B / "prototype" / "simacc_results_4e.tsv"
SIMMOSAIC = S2B / "prototype" / "simmosaic_results_4e.tsv"
SIMMOSAIC_MANIFEST = S2B / "prototype" / "simmosaic" / "manifest.tsv"

# Historical GC-sweep values from ALGORITHM_MLE.md section 4.8. These could not
# be regenerated because the source genomes for F. nucleatum, S. mutans,
# B. longum and S. coelicolor are not present in the repository.
GC_SWEEP = {
    "F. nucleatum": (27.2, 0.162, 0.125),
    "S. mutans": (36.8, 0.135, 0.063),
    "E. coli K-12": (50.8, 0.074, 0.066),
    "B. longum": (60.1, 0.356, 0.312),
    "S. coelicolor": (72.1, 0.166, 0.200),
}

C_UNIFORM = "#0072B2"
C_GAMMA = "#D55E00"
C_TRUTH = "#000000"
C_GREY = "#7F7F7F"

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.labelsize": 8,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})


def save(fig, name):
    png = FIGDIR / f"{name}.png"
    pdf = FIGDIR / f"{name}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png.name} / {pdf.name}")


def panel_label(ax, letter):
    ax.text(-0.18, 1.10, f"({letter})", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="left")


# ---------------------------------------------------------------------------
# fig_s5: indel family — ANI ladder + indel sweep
# ---------------------------------------------------------------------------

def fig_s5_indel():
    print("fig_s5_simulation_indel")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.1))

    # (a) ANI ladder with indels
    res = pd.read_csv(SIMINDEL, sep="\t")
    man = pd.read_csv(SIMINDEL_MANIFEST, sep="\t")
    df = res.merge(man[["name", "true_ani"]], left_on="query", right_on="name", how="left")
    df["truth"] = df["true_ani"] * 100.0
    df = df.sort_values("truth")

    lo, hi = 84, 101
    ax1.plot([lo, hi], [lo, hi], color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    ax1.scatter(df["truth"], df["ani"], color=C_GAMMA, marker="o",
                label="syn2bani (gamma)", zorder=3)
    ax1.scatter(df["truth"], df["ani_uniform"], color=C_UNIFORM, marker="^",
                facecolors="none", edgecolors=C_UNIFORM,
                label="syn2bani (uniform)", zorder=3)
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    ax1.set_xlabel("True ANI (%)")
    ax1.set_ylabel("Estimated ANI (%)")
    ax1.grid(True)
    err_g = (df["ani"] - df["truth"]).abs().mean()
    err_u = (df["ani_uniform"] - df["truth"]).abs().mean()
    ax1.text(0.03, 0.97, f"MAE gamma = {err_g:.3f}\nMAE uniform = {err_u:.3f}",
             transform=ax1.transAxes, va="top", ha="left", fontsize=6.5)
    ax1.legend(loc="lower right")
    panel_label(ax1, "a")

    # (b) indel sweep at 95% ANI
    res = pd.read_csv(SIMINDEL_SWEEP, sep="\t")
    man = pd.read_csv(SIMINDEL_SWEEP_MANIFEST, sep="\t")
    df = res.merge(man[["name", "indel_rate", "true_ani"]],
                   left_on="query", right_on="name").sort_values("indel_rate")
    truth = df["true_ani"].iloc[0] * 100.0

    ax2.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax2.plot(df["indel_rate"], df["ani"] - truth, color=C_GAMMA, marker="o",
             label="syn2bani (gamma)")
    ax2.plot(df["indel_rate"], df["ani_uniform"] - truth, color=C_UNIFORM,
             marker="^", label="syn2bani (uniform)")
    ax2.set_xlabel("Indel rate (deletions per 100 kb)")
    ax2.set_ylabel("Error (est − truth, ANI points)")
    ax2.grid(True)
    ax2.legend(loc="best")
    panel_label(ax2, "b")

    fig.tight_layout()
    save(fig, "fig_s5_simulation_indel")


# ---------------------------------------------------------------------------
# fig_s6: GC family — historical sweep (source genomes unavailable)
# ---------------------------------------------------------------------------

def fig_s6_gc():
    print("fig_s6_simulation_gc")
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.1))

    gc = np.array([v[0] for v in GC_SWEEP.values()])
    mae4 = np.array([v[1] for v in GC_SWEEP.values()])
    mae5 = np.array([v[2] for v in GC_SWEEP.values()])

    ax.plot(gc, mae4, color=C_GAMMA, marker="o", label="current 4-enzyme panel")
    ax.plot(gc, mae5, color=C_UNIFORM, marker="^", label="balanced 5-enzyme panel")
    for x, y, name in zip(gc, mae4, GC_SWEEP):
        ax.annotate(name, (x, y), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=5.5, color=C_GREY)
    ax.set_xlabel("Genome GC content (%)")
    ax.set_ylabel("MAE vs known truth (ANI points)")
    ax.set_ylim(0, None)
    ax.grid(True)
    ax.legend(loc="upper left")
    panel_label(ax, "a")

    fig.tight_layout()
    save(fig, "fig_s6_simulation_gc")


# ---------------------------------------------------------------------------
# fig_s7: fragmentation family
# ---------------------------------------------------------------------------

def fig_s7_fragment():
    print("fig_s7_simulation_fragment")
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.1))

    frag = pd.read_csv(SIMFRAG, sep="\t")
    frag["n_contigs"] = frag["query"].str.extract(r"q95_c(\d+)_")[0].astype(int)
    frag = frag.sort_values("n_contigs")
    truth = 95.000

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.plot(frag["n_contigs"], frag["ani"] - truth, color=C_GAMMA, marker="o",
            label="syn2bani (gamma)")
    ax.plot(frag["n_contigs"], frag["ani_uniform"] - truth, color=C_UNIFORM,
            marker="^", label="syn2bani (uniform)")
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200])
    ax.set_xticklabels([20, 50, 100, 200])
    ax.minorticks_off()
    ax.set_xlabel("Number of contigs")
    ax.set_ylabel("Error (est − 95.000, ANI points)")
    ax.grid(True)
    ax.legend(loc="best")
    panel_label(ax, "a")

    fig.tight_layout()
    save(fig, "fig_s7_simulation_fragment")


# ---------------------------------------------------------------------------
# fig_s8: accessory family
# ---------------------------------------------------------------------------

def fig_s8_accessory():
    print("fig_s8_simulation_accessory")
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.1))

    acc = pd.read_csv(SIMACC, sep="\t")
    acc["frac"] = acc["query"].str.extract(r"acc(\d+\.\d+)")[0].astype(float) * 100.0
    acc = acc.sort_values("frac")
    truth = 95.000

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.plot(acc["frac"], acc["ani"] - truth, color=C_GAMMA, marker="o",
            label="ANI error (gamma)")
    ax.plot(acc["frac"], acc["ani_uniform"] - truth, color=C_UNIFORM, marker="^",
            label="ANI error (uniform)")
    ax.set_xlabel("Accessory fraction (%)")
    ax.set_ylabel("Error (est − 95.000, ANI points)", color=C_GAMMA)
    ax.tick_params(axis="y", labelcolor=C_GAMMA)
    ax.set_ylim(-0.3, 0.35)
    ax.grid(True)

    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.plot(acc["frac"], acc["af_query"], color=C_GREY, marker="s", ls=":",
             label="af_query")
    ax2.plot(acc["frac"], 1 - acc["frac"] / 100.0, color=C_TRUTH, lw=0.8,
             ls="--", label="1 − F (truth)")
    ax2.set_ylabel("Aligned fraction (af_query)", color=C_GREY)
    ax2.tick_params(axis="y", labelcolor=C_GREY)
    ax2.set_ylim(0.4, 1.05)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=6)
    panel_label(ax, "a")

    fig.tight_layout()
    save(fig, "fig_s8_simulation_accessory")


# ---------------------------------------------------------------------------
# fig_s9: mosaic family — gamma/bimodal rate heterogeneity
# ---------------------------------------------------------------------------

def fig_s9_mosaic():
    print("fig_s9_simulation_mosaic")
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.1))

    res = pd.read_csv(SIMMOSAIC, sep="\t")
    man = pd.read_csv(SIMMOSAIC_MANIFEST, sep="\t")
    df = res.merge(man[["name", "true_ani", "regime", "param"]],
                   left_on="query", right_on="name", how="left")
    df["truth"] = df["true_ani"] * 100.0
    df = df.sort_values("truth")

    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    for regime, marker, color in (("gamma", "o", C_GAMMA), ("bimodal", "s", C_UNIFORM)):
        sub = df[df["regime"] == regime]
        ax.scatter(sub["truth"], sub["ani"] - sub["truth"], color=color,
                   marker=marker, label=f"{regime} (gamma est.)", zorder=3)
    ax.set_xlabel("True ANI (%)")
    ax.set_ylabel("Error (est − truth, ANI points)")
    ax.grid(True)
    ax.legend(loc="best")
    panel_label(ax, "a")

    fig.tight_layout()
    save(fig, "fig_s9_simulation_mosaic")


def main():
    fig_s5_indel()
    fig_s6_gc()
    fig_s7_fragment()
    fig_s8_accessory()
    fig_s9_mosaic()
    print("done.")


if __name__ == "__main__":
    main()
