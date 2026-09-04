#!/usr/bin/env python3
"""Generate Supplementary Figure S5: GC-content substitution-ladder sweep.

The GC sweep source genomes (*F. nucleatum*, *S. mutans*, *B. longum*,
*S. coelicolor*) are not present in the repository, so this script reproduces
the historical 4-enzyme sweep reported in `ALGORITHM_MLE.md` §4.8.  The 5-enzyme
balanced-panel line is included for reference.

Outputs
-------
paper/figures/supplementary/fig_s5_simulation_gc.{png,pdf}
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

# ---------------------------------------------------------------------------
# Historical GC-sweep values
# ---------------------------------------------------------------------------
# (genome GC%, MAE current 4-enzyme panel, MAE balanced 5-enzyme panel)
GC_SWEEP = {
    "F. nucleatum": (27.2, 0.162, 0.125),
    "S. mutans": (36.8, 0.135, 0.063),
    "E. coli K-12": (50.8, 0.074, 0.066),
    "B. longum": (60.1, 0.356, 0.312),
    "S. coelicolor": (72.1, 0.166, 0.200),
}

# Order left-to-right by increasing GC
ORDER = ["F. nucleatum", "S. mutans", "E. coli K-12", "B. longum", "S. coelicolor"]


def main():
    set_publication_style()

    gc = np.array([GC_SWEEP[name][0] for name in ORDER])
    mae4 = np.array([GC_SWEEP[name][1] for name in ORDER])
    mae5 = np.array([GC_SWEEP[name][2] for name in ORDER])

    fig, ax = plt.subplots(1, 1, figsize=figure_size(8.5, aspect=0.80))

    ax.plot(gc, mae4, color=COLORS["vermillion"], marker="o", ms=6,
            label="current 4-enzyme panel", zorder=3)
    ax.plot(gc, mae5, color=COLORS["sky_blue"], marker="^", ms=6,
            label="balanced 5-enzyme panel", zorder=3)

    # Genome labels next to the 4-enzyme points
    for x, y, name in zip(gc, mae4, ORDER):
        ax.annotate(
            name,
            (x, y),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            va="bottom",
            fontsize=7,
            color=COLORS["grey"],
        )

    ax.set_xlabel("Genome GC content (%)")
    ax.set_ylabel("MAE vs known truth (ANI points)")
    ax.set_ylim(0.0, max(mae4.max(), mae5.max()) * 1.18)

    # Summary annotation tied to the 4-enzyme panel (the current estimator)
    range_text = (
        f"Syn2bANI error: {mae4.min():.3f}–{mae4.max():.3f} ANI points\n"
        f"worst: B. longum (GC {GC_SWEEP['B. longum'][0]:.1f}%)"
    )
    ax.text(
        0.97, 0.03, range_text,
        transform=ax.transAxes,
        va="bottom", ha="right",
        fontsize=7,
        linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor=COLORS["light_grey"], alpha=0.95),
    )

    ax.legend(loc="upper left", frameon=False, fontsize=7)
    label_panel(ax, "a")

    outdir = Path(__file__).resolve().parents[1] / "paper" / "figures" / "supplementary"
    outdir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, outdir / "fig_s5_simulation_gc")
    plt.close(fig)
    print("Figure S5 (GC ladder) generated successfully.")


if __name__ == "__main__":
    main()
