#!/usr/bin/env python3
"""Diagnostic scatter: af_query vs anchor_adjacency, coloured by dnadiff synteny.

Shows that the two Syn2bANI structural metrics measure different things and
identifies four interpretable quadrants/failure modes on GTDB-R207 held-out pairs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "gtdb50k"
OUT = ROOT / "figures" / "report"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    s2b = pd.read_csv(DATA / "s2b_50k.tsv", sep="\t")
    dna = pd.read_csv(DATA / "sv_truth_50k.tsv", sep="\t")
    df = s2b[["pairid", "ani_gated", "af_query", "anchor_adjacency",
              "synteny_blocks", "breakpoint_count", "n_anchors"]].merge(dna, on="pairid")

    # Failure-mode masks
    high_adj_low_syn = (df.anchor_adjacency > 0.99) & (df.dnadiff_anchor_adjacency < 0.5)
    low_adj_high_syn = (df.anchor_adjacency < 0.95) & (df.dnadiff_anchor_adjacency > 0.8)
    low_af_high_adj = (df.af_query < 0.3) & (df.anchor_adjacency > 0.99)
    high_af_low_syn = (df.af_query > 0.8) & (df.dnadiff_anchor_adjacency < 0.5)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # --- Panel (a): density scatter coloured by dnadiff synteny ---
    ax = axes[0]
    # subsample for visibility if needed
    n = len(df)
    if n > 20000:
        rng = np.random.default_rng(7)
        idx = rng.choice(n, size=20000, replace=False)
        plot_df = df.iloc[idx]
    else:
        plot_df = df

    sc = ax.scatter(plot_df["af_query"], plot_df["anchor_adjacency"],
                    c=plot_df["dnadiff_anchor_adjacency"],
                    cmap="viridis_r", s=8, alpha=0.6, edgecolors="none",
                    vmin=0.0, vmax=1.0)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("dnadiff synteny (alignment-based)", fontsize=10)
    ax.set_xlabel("Syn2bANI af_query (base-pair chain coverage)", fontsize=11)
    ax.set_ylabel("Syn2bANI anchor_adjacency", fontsize=11)
    ax.set_title("(a) Two Syn2bANI metrics capture different structures", fontsize=12)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    # --- Panel (b): same axes, highlight failure modes ---
    ax = axes[1]
    ax.scatter(df["af_query"], df["anchor_adjacency"],
               c="#dddddd", s=4, alpha=0.5, edgecolors="none", label=f"all pairs (n={n:,})")

    colors = {
        "A: high adj, low dnadiff synteny": ("#d62728", high_adj_low_syn),
        "B: high af, low dnadiff synteny": ("#ff7f0e", high_af_low_syn),
        "C: low af, high adj": ("#2ca02c", low_af_high_adj),
        "D: low adj, high dnadiff synteny": ("#9467bd", low_adj_high_syn),
    }
    for label, (color, mask) in colors.items():
        cnt = mask.sum()
        ax.scatter(df.loc[mask, "af_query"], df.loc[mask, "anchor_adjacency"],
                   c=color, s=18, alpha=0.7, edgecolors="none",
                   label=f"{label} (n={cnt:,})")

    ax.set_xlabel("Syn2bANI af_query (base-pair chain coverage)", fontsize=11)
    ax.set_ylabel("Syn2bANI anchor_adjacency", fontsize=11)
    ax.set_title("(b) Interpretable deviation quadrants", fontsize=12)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="lower right", fontsize=8, markerscale=2)

    plt.tight_layout()
    out_png = OUT / "fig_diagnostic_af_query_vs_anchor_adjacency.png"
    out_pdf = OUT / "fig_diagnostic_af_query_vs_anchor_adjacency.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png} and {out_pdf}")

    # Text summary for manuscript
    summary = (
        "# af_query vs anchor_adjacency diagnostic summary\n\n"
        f"Total pairs: {n:,}\n\n"
        "Failure modes:\n"
        f"- A. high anchor_adjacency (>0.99) but low dnadiff synteny (<0.5): {high_adj_low_syn.sum():,}\n"
        f"- B. high af_query (>0.8) but low dnadiff synteny (<0.5): {high_af_low_syn.sum():,}\n"
        f"- C. low af_query (<0.3) but high anchor_adjacency (>0.99): {low_af_high_adj.sum():,}\n"
        f"- D. low anchor_adjacency (<0.95) but high dnadiff synteny (>0.8): {low_adj_high_syn.sum():,}\n\n"
        "Correlations with dnadiff synteny:\n"
        f"- af_query: Pearson {np.corrcoef(df.af_query, df.dnadiff_anchor_adjacency)[0,1]:.3f}\n"
        f"- anchor_adjacency: Pearson {np.corrcoef(df.anchor_adjacency, df.dnadiff_anchor_adjacency)[0,1]:.3f}\n"
    )
    out_txt = OUT / "fig_diagnostic_af_query_vs_anchor_adjacency_summary.txt"
    out_txt.write_text(summary)
    print(f"wrote {out_txt}")


if __name__ == "__main__":
    main()
