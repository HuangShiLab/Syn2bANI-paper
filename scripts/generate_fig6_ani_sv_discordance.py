#!/usr/bin/env python3
"""Generate Figure 6: ANI and SV are decoupled in GTDB-R207.

Publication-quality 2-panel figure.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS

set_publication_style()
D = Path("results/gtdb50k")
OUT = Path("paper/figures/main/fig6_ani_sv_discordance")


def main():
    # 336 high-AF same-species pairs
    pairs = pd.read_csv(D / "all_high_ani_high_af_pairs.tsv", sep="\t")
    pairs = pairs.rename(columns={"synteny_score": "anchor_adjacency"})

    # Top-20 discordant pairs for highlight (lowest anchor adjacency among high-ANI pairs)
    sv = pd.read_csv(D / "struct_top_discordant_summary.tsv", sep="\t")
    discordant_ids = set(sv["case"].head(20))
    pairs["discordant"] = pairs["pairid"].isin(discordant_ids)

    sv_top = sv.head(10).copy()
    # Keep full taxonomic labels on a single line; use small font to avoid overlap.
    sv_top["pair_label"] = sv_top["query_species"]

    fig, axes = plt.subplots(1, 2, figsize=figure_size(17.8, aspect=0.55))

    # (a) ANIm vs anchor adjacency
    ax = axes[0]
    ax.scatter(
        pairs["anim_ani"], pairs["anchor_adjacency"],
        c=COLORS["grey"], s=12, alpha=0.6, edgecolors="none", label="High-AF same-species pairs"
    )
    disc = pairs[pairs["discordant"]]
    ax.scatter(
        disc["anim_ani"], disc["anchor_adjacency"],
        c=COLORS["vermillion"], s=20, alpha=0.9, edgecolors="none", label="Top 20 discordant pairs"
    )
    ax.axhline(0.95, color=COLORS["black"], ls="--", lw=0.8)
    ax.axvline(99, color=COLORS["black"], ls="--", lw=0.8)
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("Anchor adjacency")
    ax.set_xlim(98.8, 100.05)
    ax.set_ylim(0.78, 1.005)
    label_panel(ax, "a")
    ax.set_title(f"GTDB-R207 high-AF same-species pairs (n = {len(pairs)})", fontsize=8)
    ax.legend(loc="lower left", fontsize=7, frameon=False)

    # (b) SV composition of top-10 discordant pairs (horizontal bars)
    ax = axes[1]
    y = np.arange(len(sv_top))
    height = 0.6
    left = np.zeros(len(sv_top))
    colors = [COLORS["blue"], COLORS["orange"], COLORS["bluish_green"]]
    labels = ["Inversions", "Translocations", "Indels"]
    for col, color, label in zip(["n_inversions", "n_translocations", "n_indels"], colors, labels):
        vals = sv_top[col].values
        ax.barh(y, vals, height, left=left, color=color, label=label, edgecolor="white", linewidth=0.3)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(sv_top["pair_label"], fontsize=5)
    ax.set_xlabel("SV count")
    ax.set_ylabel("Top discordant pairs")
    ax.set_xlim(0, left.max() * 1.15)
    ax.invert_yaxis()
    label_panel(ax, "b")
    ax.set_title("SV composition of 10 most discordant pairs", fontsize=8)
    ax.legend(loc="lower right", fontsize=7, frameon=False)

    plt.tight_layout()
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
