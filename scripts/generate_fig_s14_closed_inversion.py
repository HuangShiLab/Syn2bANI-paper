#!/usr/bin/env python3
"""Generate Supplementary Fig. S14: closed-genome all-vs-all inversion
fraction, raw vs orientation-corrected.

Context
-------
701 near-complete GTDB-R207 genomes (4 species) were compared all-vs-all.
The reference-oriented raw inverted fraction has median 0.36 and q90 0.93,
which is biologically implausible: in undirected all-vs-all mode the two
assemblies are randomly oriented, so whole chromosomes get classified as
inverted. The corrected column ``syn2b_inverted_fraction = min(raw, 1-raw)``
removes this global-orientation artifact (median 0.184).

Panels
------
(a) Distribution of the raw inverted fraction; median and q90 marked,
    with the raw > 0.5 mirror region shaded (38.9% of pairs).
(b) Distribution of the corrected inverted fraction (bounded at 0.5);
    median marked.
(c) Corrected vs raw scatter: points lie on the y = x identity arm
    (61.1% of pairs) or the y = 1 - x mirror arm (38.9% of pairs).

Outputs
-------
paper/figures/supplementary/fig_s14_closed_genome_inversion.png
paper/figures/supplementary/fig_s14_closed_genome_inversion.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Ensure the repository plot-style module is importable when running from root.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import plot_style  # noqa: E402

REPO_ROOT = SCRIPT_DIR.parent
DATA_PATH = REPO_ROOT / "results" / "gtdb50k" / "syn2b_inverted_fraction_closed.tsv"
GENOMES_PATH = REPO_ROOT / "results" / "gtdb50k" / "closed_inversions_genomes.tsv"
OUT_DIR = REPO_ROOT / "paper" / "figures" / "supplementary"
OUT_STEM = OUT_DIR / "fig_s14_closed_genome_inversion"


def load_pairs() -> pd.DataFrame:
    """Load the closed-genome all-vs-all pairs, keeping status == 'ok' rows."""
    df = pd.read_csv(DATA_PATH, sep="\t", low_memory=False)
    df = df[df["status"] == "ok"].copy()
    df["syn2b_raw_inverted_fraction"] = df["syn2b_raw_inverted_fraction"].astype(float)
    df["syn2b_inverted_fraction"] = df["syn2b_inverted_fraction"].astype(float)

    # Sanity check: the cohort spans 4 species; confirm every kept pair is
    # same-species (the diagnostic report quotes 61,537 same-species pairs).
    genomes = pd.read_csv(GENOMES_PATH, sep="\t")
    species = dict(zip(genomes["acc"], genomes["species"]))
    q_sp = df["q_acc"].map(species)
    r_sp = df["r_acc"].map(species)
    if not (q_sp == r_sp).all():
        raise ValueError("Unexpected: cohort contains cross-species pairs.")
    df["species"] = q_sp
    return df


def _label_above_axes(ax, x: float, text: str, ha: str) -> None:
    """Place a small label just above the axes top edge at data coordinate x."""
    dx = 3 if ha == "left" else -3
    ax.annotate(text, xy=(x, 1.0), xycoords=("data", "axes fraction"),
                xytext=(dx, 4), textcoords="offset points",
                ha=ha, va="bottom", fontsize=7)


def plot_panel_raw(ax, raw: pd.Series) -> None:
    """(a) Raw inverted-fraction distribution with median / q90 lines."""
    median = raw.median()
    q90 = raw.quantile(0.9)
    n = len(raw)
    n_mirror = int((raw > 0.5).sum())
    frac_mirror = 100.0 * n_mirror / n

    bins = np.linspace(0, 1, 51)
    ax.hist(raw, bins=bins, color=plot_style.COLORS["vermillion"],
            alpha=0.8, edgecolor="white", linewidth=0.3)

    # Mirror region: raw > 0.5 means the two genomes are globally
    # reverse-complemented, so the whole chromosome is scored as inverted.
    ax.axvspan(0.5, 1.0, color=plot_style.COLORS["light_grey"], alpha=0.45, zorder=0)

    ax.axvline(median, color="black", linewidth=1.2, zorder=3)
    ax.axvline(q90, color="black", linewidth=1.0, linestyle="--", zorder=3)
    _label_above_axes(ax, median, f"median = {median:.4f}", ha="left")
    ymax = ax.get_ylim()[1]
    ax.text(0.03, ymax * 0.97, f"n = {n:,}", ha="left", va="top", fontsize=7)
    ax.text(0.91, ymax * 0.78, f"q90 = {q90:.4f}", ha="right", va="top", fontsize=7)
    ax.text(0.75, ymax * 0.97, f"{frac_mirror:.1f}% of pairs\nraw > 0.5 (mirror)",
            ha="center", va="top", fontsize=7)

    ax.set_xlim(0, 1)
    ax.set_xlabel("Raw inverted fraction (reference-oriented)")
    ax.set_ylabel("Pairs")
    ax.set_title("Raw metric", fontsize=9, pad=13)


def plot_panel_corrected(ax, corrected: pd.Series) -> None:
    """(b) Corrected inverted-fraction distribution with median line."""
    median = corrected.median()
    n = len(corrected)

    bins = np.linspace(0, 0.5, 51)
    ax.hist(corrected, bins=bins, color=plot_style.COLORS["bluish_green"],
            alpha=0.8, edgecolor="white", linewidth=0.3)

    ax.axvline(median, color="black", linewidth=1.2, zorder=3)
    _label_above_axes(ax, median, f"median = {median:.4f}", ha="left")
    ymax = ax.get_ylim()[1]
    ax.text(0.03, ymax * 0.97, f"n = {n:,}", ha="left", va="top", fontsize=7)
    ax.text(0.485, ymax * 0.97, "min(raw, 1 − raw)", ha="right", va="top", fontsize=7)

    ax.set_xlim(0, 0.5)
    ax.set_xlabel("Corrected inverted fraction")
    ax.set_ylabel("Pairs")
    ax.set_title("Corrected metric", fontsize=9, pad=13)


def plot_panel_scatter(ax, raw: pd.Series, corrected: pd.Series) -> None:
    """(c) Corrected vs raw: identity arm vs mirror arm."""
    n = len(raw)
    n_mirror = int((raw > 0.5).sum())
    frac_mirror = 100.0 * n_mirror / n
    frac_identity = 100.0 - frac_mirror

    ax.scatter(raw, corrected, s=1, alpha=0.12, linewidths=0,
               color=plot_style.COLORS["blue"], rasterized=True, zorder=2)

    xs = np.linspace(0, 1, 100)
    line_id, = ax.plot(xs, xs, color=plot_style.COLORS["grey"], linestyle="--",
                       linewidth=0.9, zorder=3)
    line_mirror, = ax.plot(xs, 1.0 - xs, color=plot_style.COLORS["orange"],
                           linestyle="--", linewidth=0.9, zorder=3)

    ax.legend(
        [line_id, line_mirror],
        [f"identity arm: y = x ({frac_identity:.1f}%)",
         f"mirror arm: y = 1 − x ({frac_mirror:.1f}%)"],
        loc="upper left", fontsize=7, handlelength=1.8,
        bbox_to_anchor=(0.01, 0.99),
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 0.62)
    ax.set_xlabel("Raw inverted fraction")
    ax.set_ylabel("Corrected inverted fraction")
    ax.set_title("Correction folds the mirror arm", fontsize=9, pad=13)


def main() -> int:
    plot_style.set_publication_style()

    df = load_pairs()
    raw = df["syn2b_raw_inverted_fraction"]
    corrected = df["syn2b_inverted_fraction"]

    # Recomputed statistics (printed for the record; the diagnostic report
    # quotes the same numbers).
    n = len(df)
    n_mirror = int((raw > 0.5).sum())
    print(f"n pairs (status == ok): {n}")
    print(f"raw:      median {raw.median():.4f}  q90 {raw.quantile(0.9):.4f}  "
          f"q95 {raw.quantile(0.95):.4f}  mean {raw.mean():.4f}  max {raw.max():.4f}")
    print(f"corrected: median {corrected.median():.4f}  q90 {corrected.quantile(0.9):.4f}  "
          f"mean {corrected.mean():.4f}  max {corrected.max():.4f}")
    print(f"mirror fraction (raw > 0.5): {n_mirror} / {n} = {100.0 * n_mirror / n:.2f}%")
    max_dev = float((corrected - np.minimum(raw, 1 - raw)).abs().max())
    print(f"max |corrected - min(raw, 1-raw)|: {max_dev:.2e}")

    fig = plt.figure(figsize=plot_style.figure_size(17.8, aspect=0.42))
    gs = fig.add_gridspec(1, 3, left=0.065, right=0.985, top=0.88, bottom=0.20,
                          wspace=0.38)

    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_raw(ax_a, raw)
    plot_style.label_panel(ax_a, "a")

    ax_b = fig.add_subplot(gs[0, 1])
    plot_panel_corrected(ax_b, corrected)
    plot_style.label_panel(ax_b, "b")

    ax_c = fig.add_subplot(gs[0, 2])
    plot_panel_scatter(ax_c, raw, corrected)
    plot_style.label_panel(ax_c, "c")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_style.save_figure(fig, str(OUT_STEM), formats=("png", "pdf"))
    plt.close(fig)

    print("Figure S14 generated successfully.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
