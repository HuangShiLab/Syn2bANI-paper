#!/usr/bin/env python3
"""Generate Supplementary Figure S10: FDA-ARGOS S. aureus negative control.

Uses the self-contained Syn2bANI v0.1.1 triangle rerun. The collection has
67 deduplicated genomes and 2,211 non-self pairs. Panel (a) colors pairs whose
two assemblies share an isolation country and greys all cross-country or
missing-country pairs. Panel (b) shows the breakpoint-count distribution.
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "case_studies" / "fda_argos_s_aureus" / "results"
TRIANGLE_PATH = RESULTS_DIR / "triangle_rerun.tsv"
METADATA_PATH = RESULTS_DIR / "assembly_metadata.tsv"
OUT = ROOT / "paper" / "figures" / "supplementary" / "fig_s10_saureus_breakpoints"


def load_data():
    tri = pd.read_csv(TRIANGLE_PATH, sep="\t", low_memory=False)
    meta = pd.read_csv(METADATA_PATH, sep="\t", low_memory=False)
    country = meta.set_index("assembly_acc")["country"].to_dict()

    def shared(row):
        a = country.get(str(row["q_acc"]), "missing")
        b = country.get(str(row["r_acc"]), "missing")
        if a in ("", "missing") or b in ("", "missing") or a != b:
            return "different/missing country"
        return a

    tri["country_label"] = tri.apply(shared, axis=1)
    return tri


def main():
    set_publication_style()
    tri = load_data()
    if tri.duplicated(["query", "reference"]).any():
        raise ValueError("Duplicate unordered pairs in triangle_rerun.tsv")
    if len(tri) != 2211:
        raise ValueError(f"Expected 2,211 rerun pairs; found {len(tri)}")

    fig, axes = plt.subplots(1, 2, figsize=figure_size(22.0, aspect=0.60))

    ax = axes[0]
    for label, sub in tri.groupby("country_label", sort=True):
        grey = label == "different/missing country"
        ax.scatter(
            sub["breakpoint_count"], sub["ani"],
            c=[{"color": "lightgrey"} if grey else {} for _ in [0]][0]["color"] if grey else None,
            s=12, alpha=0.65 if grey else 0.8, edgecolors="none", label=label,
        )
    ax.set_xlabel("Breakpoint count")
    ax.set_ylabel("Syn2bANI ANI (%)")
    ax.set_ylim(92.5, 100.5)
    ax.set_xlim(-0.4, 8.4)
    ax.set_xticks(range(0, 9, 1))
    ax.legend(loc="upper left", fontsize=6, frameon=False, title="Country",
              title_fontsize=6.5, ncol=2)
    label_panel(ax, "a")
    ax.set_title("ANI versus breakpoint count", fontsize=9)

    ax = axes[1]
    ax.hist(tri["breakpoint_count"], bins=np.arange(-0.5, 9.5, 1),
            color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Breakpoint count")
    ax.set_ylabel("Pairs")
    ax.set_xticks(range(0, 9, 1))
    label_panel(ax, "b")
    ax.set_title(
        f"{(tri['breakpoint_count'] == 0).mean():.1%} of pairs have zero breakpoints",
        fontsize=9,
    )

    outdir = OUT.parent
    outdir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, str(OUT))
    print(
        f"n={len(tri)}, zero={(tri['breakpoint_count'] == 0).sum()}, "
        f"max={tri['breakpoint_count'].max()}, "
        f"ANI={tri['ani'].min():.1f}-{tri['ani'].max():.1f}"
    )


if __name__ == "__main__":
    main()
