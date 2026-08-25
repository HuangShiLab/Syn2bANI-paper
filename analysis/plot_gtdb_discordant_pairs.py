#!/usr/bin/env python3
"""Publication-quality figure for GTDB-R207 high-ANI discordant pairs.

Shows that database-scale ANI-based searches can return top hits with
near-clonal ANI yet substantially degraded synteny / extensive rearrangements.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "gtdb50k"
OUT_FIG = ROOT / "figures" / "gtdb50k"
OUT_FIG.mkdir(parents=True, exist_ok=True)


def abbrev_name(name, max_len=35):
    """Shorten species name for plot labels."""
    if len(name) <= max_len:
        return name
    return name[:max_len - 3] + "..."


def main():
    # All high-ANI high-AF pairs (background)
    all_pairs = pd.read_csv(DATA / "all_high_ani_high_af_pairs.tsv", sep="\t")

    # Top discordant pairs with struct-level SV counts
    struct = pd.read_csv(DATA / "struct_top_discordant_summary.tsv", sep="\t")
    struct["discordance"] = struct["ani"] - 100 * struct["synteny_score"]
    struct = struct.sort_values("discordance", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # --- Panel (a): ANI vs synteny_score scatter ---
    ax = axes[0]
    ax.scatter(all_pairs["ani"], all_pairs["synteny_score"],
               c="#cccccc", s=25, alpha=0.6, edgecolors="none", label="GTDB high-AF pairs (n=336)")
    ax.scatter(struct["ani"], struct["synteny_score"],
               c="#d62728", s=60, alpha=0.85, edgecolors="black", linewidths=0.5,
               label=f"Top discordant pairs (n={len(struct)})")
    ax.set_xlabel("Syn2bANI ANI (%)", fontsize=11)
    ax.set_ylabel("Synteny score", fontsize=11)
    ax.set_title("(a) Near-clonal ANI can coexist with low synteny", fontsize=12)
    ax.set_xlim(98.95, 100.005)
    ax.set_ylim(0.815, 1.005)
    ax.axhline(0.95, color="gray", ls="--", lw=0.8, alpha=0.7)
    ax.text(99.0, 0.948, "synteny = 0.95", fontsize=8, color="gray", va="top")
    ax.legend(loc="lower left", fontsize=8)

    # annotate top 2 most discordant with non-overlapping offsets
    top2 = struct.head(2)
    offsets = [(30, 22), (-30, -22)]
    for i, (_, r) in enumerate(top2.iterrows()):
        ox, oy = offsets[i]
        label = f"{r['query_species'][:22]}\nvs {r['reference_species'][:22]}"
        ax.annotate(label, xy=(r["ani"], r["synteny_score"]),
                    xytext=(ox, oy), textcoords="offset points",
                    fontsize=6, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="black", alpha=0.85))

    # --- Panel (b): SV breakdown for top 10 discordant pairs ---
    ax = axes[1]
    top10 = struct.head(10).copy()
    # Build readable labels (species only, truncated)
    def short_label(qs, rs, max_len=28):
        q = qs.split()[-1] if len(qs.split()) > 1 else qs
        r = rs.split()[-1] if len(rs.split()) > 1 else rs
        base = f"{q} vs {r}"
        return base if len(base) <= max_len else base[:max_len - 3] + "..."

    top10["label"] = [
        short_label(r["query_species"], r["reference_species"])
        for _, r in top10.iterrows()
    ]
    y_pos = np.arange(len(top10))

    inv = top10["n_inversions"].values
    tra = top10["n_translocations"].values
    indel = top10["n_indels"].values

    ax.barh(y_pos, inv, color="#1f77b4", label="Inversions")
    ax.barh(y_pos, tra, left=inv, color="#ff7f0e", label="Translocations")
    ax.barh(y_pos, indel, left=inv + tra, color="#2ca02c", label="Indels")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top10["label"].values, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Number of SV calls", fontsize=11)
    ax.set_title("(b) SV composition of the 10 most discordant pairs", fontsize=12)
    ax.legend(loc="lower right", fontsize=8)

    # annotate total breakpoints
    for i, (_, r) in enumerate(top10.iterrows()):
        total = inv[i] + tra[i] + indel[i]
        ax.text(total + max(top10["breakpoint_count"].max() * 0.02, 10), i,
                f"{int(r['breakpoint_count'])} breaks", va="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(OUT_FIG / "gtdb_discordant_high_ani.png", dpi=300)
    fig.savefig(OUT_FIG / "gtdb_discordant_high_ani.pdf")
    plt.close(fig)

    print(f"Wrote {OUT_FIG / 'gtdb_discordant_high_ani.png'}")
    print("Top 5 discordant GTDB pairs:")
    print(top10.head(5)[["case", "query_species", "reference_species", "ani",
                         "synteny_score", "breakpoint_count", "n_inversions",
                         "n_translocations", "n_indels"]].to_string(index=False))


if __name__ == "__main__":
    main()
