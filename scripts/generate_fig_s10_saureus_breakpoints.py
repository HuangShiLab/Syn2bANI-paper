#!/usr/bin/env python3
"""Generate Supplementary Figure S10: High-ANI FDA-ARGOS *S. aureus* pairs
show wide breakpoint variation.

One hundred and twenty-two genomes, 7,381 refined pairs.
(a) ANI vs breakpoint count colored by country.
(b) ANI vs breakpoint count colored by isolation source.
Pairs at ~100% ANI still carry >150 breakpoints.
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure

set_publication_style()

# Inputs
RESULTS_DIR = Path("case_studies/fda_argos_s_aureus/results")
GENOMES_DIR = Path("/Volumes/MoneyCat/Data/fda_argos_staphylococcus_aureus/genomes")
TRIANGLE_PATH = RESULTS_DIR / "triangle.tsv"
METADATA_PATH = RESULTS_DIR / "assembly_metadata.tsv"

# Output
OUT = Path("paper/figures/supplementary/fig_s10_saureus_breakpoints")

MIN_CATEGORY_COUNT = 5


def parse_fasta_headers(fasta_path: Path):
    """Return all sequence IDs from a FASTA file."""
    seq_ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                seq_ids.append(seq_id)
    return seq_ids


def build_seq_to_assembly_map(genomes_dir: Path):
    """Map every sequence ID to its assembly accession (FASTA stem)."""
    mapping = {}
    for fna in sorted(genomes_dir.glob("*.fna")):
        acc = fna.stem
        for seq_id in parse_fasta_headers(fna):
            mapping[seq_id] = acc
    return mapping


def group_categories(series: pd.Series, min_count: int = MIN_CATEGORY_COUNT):
    """Collapse rare categories into an 'other' bin."""
    counts = series.value_counts()
    keep = set(counts[counts >= min_count].index)
    return series.apply(lambda x: x if x in keep else "other")


def load_and_annotate(triangle_path: Path, metadata_path: Path, genomes_dir: Path):
    """Load triangle output and attach query-side metadata."""
    seq_map = build_seq_to_assembly_map(genomes_dir)
    meta = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    meta["assembly_acc"] = meta["assembly_acc"].astype(str)

    tri = pd.read_csv(triangle_path, sep="\t", low_memory=False)

    # Handle duplicated column name ('flag' appears twice in some outputs).
    if isinstance(tri.columns, pd.Index):
        cols = list(tri.columns)
        seen = set()
        renamed = []
        for c in cols:
            if c in seen:
                renamed.append(f"{c}_dup")
            else:
                renamed.append(c)
                seen.add(c)
        tri.columns = renamed

    tri["query_acc"] = tri["query"].map(seq_map)
    tri["ref_acc"] = tri["reference"].map(seq_map)

    missing_query = tri["query_acc"].isna().sum()
    missing_ref = tri["ref_acc"].isna().sum()
    if missing_query or missing_ref:
        print(f"Warning: {missing_query} query / {missing_ref} reference mapping failures")
        tri = tri.dropna(subset=["query_acc", "ref_acc"]).copy()

    # Attach query metadata as the representative label for each pair.
    tri = tri.merge(
        meta,
        left_on="query_acc",
        right_on="assembly_acc",
        how="left",
        suffixes=("", "_query"),
    )

    tri["country_grouped"] = group_categories(tri["country"].fillna("unknown"))
    tri["source_grouped"] = group_categories(tri["isolation_source"].fillna("unknown"))
    return tri


def scatter_by_category(
    ax,
    df,
    category_col,
    xlabel,
    legend_title,
    legend_loc="upper left",
    legend_anchor=(1.02, 1.0),
    ncol=1,
    legend_fontsize=6,
):
    """Scatter ANI vs breakpoint_count colored by a categorical column."""
    categories = sorted(df[category_col].unique())
    # Use a categorical colormap with enough distinct colors.
    cmap = plt.cm.tab20
    colors = cmap(np.linspace(0, 1, len(categories)))
    color_map = dict(zip(categories, colors))

    for cat in categories:
        sub = df[df[category_col] == cat]
        ax.scatter(
            sub["breakpoint_count"],
            sub["ani"],
            c=[color_map[cat]],
            s=12,
            alpha=0.7,
            edgecolors="none",
            label=cat,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("ANI (%)")
    ax.set_ylim(92.5, 100.5)
    ax.set_xlim(-10, 430)
    leg = ax.legend(
        loc=legend_loc,
        bbox_to_anchor=legend_anchor,
        fontsize=legend_fontsize,
        handletextpad=0.15,
        borderaxespad=0.15,
        labelspacing=0.25,
        markerscale=0.9,
        title=legend_title,
        title_fontsize=legend_fontsize + 0.5,
        ncol=ncol,
    )
    leg.set_clip_on(False)


def main():
    tri = load_and_annotate(TRIANGLE_PATH, METADATA_PATH, GENOMES_DIR)
    print(f"Loaded {len(tri)} refined pairs")
    print(f"  ANI range: {tri['ani'].min():.4f} - {tri['ani'].max():.4f}")
    print(f"  Breakpoint range: {tri['breakpoint_count'].min()} - {tri['breakpoint_count'].max()}")

    n_high = ((tri["ani"] >= 99.9) & (tri["breakpoint_count"] > 150)).sum()
    print(f"  Pairs with ANI >= 99.9 and breakpoints > 150: {n_high}")

    # Supplementary figure: wide enough for two-column legends below each panel.
    fig, axes = plt.subplots(1, 2, figsize=figure_size(22.0, aspect=0.60))

    # Panel (a): colored by country
    scatter_by_category(
        axes[0],
        tri,
        "country_grouped",
        "Breakpoint count",
        legend_title="Country",
        legend_loc="upper center",
        legend_anchor=(0.5, -0.18),
        ncol=2,
        legend_fontsize=5.5,
    )
    label_panel(axes[0], "a")

    # Panel (b): colored by isolation source
    scatter_by_category(
        axes[1],
        tri,
        "source_grouped",
        "Breakpoint count",
        legend_title="Isolation source",
        legend_loc="upper center",
        legend_anchor=(0.5, -0.18),
        ncol=2,
        legend_fontsize=5.5,
    )
    label_panel(axes[1], "b")

    fig.subplots_adjust(left=0.06, right=0.95, bottom=0.22, top=0.92, wspace=0.25)
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
