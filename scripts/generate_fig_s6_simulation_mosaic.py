#!/usr/bin/env python3
"""Generate Supplementary Figure S6: mosaic / rate-heterogeneity simulations.

Per-block divergence rates sampled from Gamma(α, α) with α = 0.5/1/2
(mean ANI 90–98) and deliberately misspecified bimodal 50/70%-core cases.
The figure compares the gamma and uniform Syn2bANI estimators against exact
ground truth.

Data provenance
---------------
Reads rerun outputs from the Syn2bANI prototype directory:
  prototype/simmosaic_results_4e.tsv
  prototype/simmosaic/manifest.tsv

Outputs
-------
paper/figures/supplementary/fig_s6_simulation_mosaic.{png,pdf}
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
PROTOTYPE = Path("/Users/macstudio/Downloads/Syn2bANI/prototype")
RESULTS_TSV = PROTOTYPE / "simmosaic_results_4e.tsv"
MANIFEST_TSV = PROTOTYPE / "simmosaic" / "manifest.tsv"


def load_data():
    """Load mosaic results and merge with manifest to obtain truth/regime."""
    res = pd.read_csv(RESULTS_TSV, sep="\t")
    man = pd.read_csv(MANIFEST_TSV, sep="\t")
    df = res.merge(man, left_on="query", right_on="name", how="left")
    df["true_ani_pct"] = df["true_ani"] * 100.0
    df["err_gamma"] = df["ani"] - df["true_ani_pct"]
    df["err_uniform"] = df["ani_uniform"] - df["true_ani_pct"]
    return df.sort_values("true_ani_pct").reset_index(drop=True)


def regime_mae_bias(df, regime, estimator):
    """Return MAE and mean bias for a regime and estimator column."""
    sub = df[df["regime"] == regime]
    errs = sub[estimator].values
    return float(np.mean(np.abs(errs))), float(np.mean(errs))


def main():
    set_publication_style()

    if not RESULTS_TSV.exists():
        raise FileNotFoundError(
            f"Mosaic results not found at {RESULTS_TSV}. "
            "Regenerate them in the Syn2bANI prototype directory."
        )
    if not MANIFEST_TSV.exists():
        raise FileNotFoundError(f"Mosaic manifest not found at {MANIFEST_TSV}.")

    df = load_data()

    # -----------------------------------------------------------------------
    # Figure: two panels, gamma vs uniform estimator
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=figure_size(17.8, aspect=0.42))
    fig.subplots_adjust(wspace=0.30)

    regime_style = {
        "gamma": {
            "color": COLORS["vermillion"],
            "marker": "o",
            "label": "gamma regime",
        },
        "bimodal": {
            "color": COLORS["sky_blue"],
            "marker": "s",
            "label": "bimodal regime",
        },
    }

    for ax, estimator, title in [
        (axes[0], "err_gamma", "Gamma estimator"),
        (axes[1], "err_uniform", "Uniform estimator"),
    ]:
        ax.axhline(0, color=COLORS["black"], lw=0.8, ls="--", zorder=1)

        for regime in ("gamma", "bimodal"):
            sub = df[df["regime"] == regime]
            style = regime_style[regime]
            ax.scatter(
                sub["true_ani_pct"],
                sub[estimator],
                color=style["color"],
                marker=style["marker"],
                s=45,
                label=style["label"],
                zorder=3,
            )

        ax.set_xlabel("True ANI (%)")
        ax.set_ylabel("Error (est − truth, ANI points)")
        ax.set_title(title, fontsize=9)

        # Per-regime MAE/bias annotation
        lines = []
        for regime in ("gamma", "bimodal"):
            mae, bias = regime_mae_bias(df, regime, estimator)
            lines.append(
                f"{regime:<7} MAE {mae:5.2f}, bias {bias:+5.2f}"
            )
        ax.text(
            0.97, 0.97, "\n".join(lines),
            transform=ax.transAxes,
            va="top", ha="right",
            fontsize=7,
            linespacing=1.3,
            zorder=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=COLORS["light_grey"], alpha=0.95),
        )

        ax.legend(loc="lower left", frameon=False, fontsize=7, ncol=1)

    label_panel(axes[0], "a")
    label_panel(axes[1], "b")

    outdir = ROOT / "paper" / "figures" / "supplementary"
    outdir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, outdir / "fig_s6_simulation_mosaic")
    plt.close(fig)
    print("Figure S6 (mosaic / rate heterogeneity) generated successfully.")


if __name__ == "__main__":
    main()
