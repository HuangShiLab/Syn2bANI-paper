#!/usr/bin/env python3
"""Generate Supplementary Figure S9: High-ANI E. coli O157:H7 breakpoint diversity.

Caption
-------
"Seventy-four genomes from Fitzgerald et al. (2021), 2,701 non-self pairs.
(a) ANI vs breakpoint count colored by lineage (I/II, II, Ia, Ic).
(b) Same data colored by host category (bovine, human, other/unknown).
All pairwise ANIs exceed 99.886% yet breakpoints range from 171 to >1,100."

Inputs
------
case_studies/ecoli_o157_fitzgerald_2021/results/triangle.tsv
    Pairwise syn2bani results (query, reference, ani, breakpoint_count, ...).
case_studies/ecoli_o157_fitzgerald_2021/results/metadata_with_lineage.tsv
    Per-genome metadata including assigned_lineage and host_category.

Outputs
-------
paper/figures/supplementary/fig_s9_ecoli_o157_breakpoints.png
paper/figures/supplementary/fig_s9_ecoli_o157_breakpoints.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import COLORS, figure_size, label_panel, save_figure, set_publication_style

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
CASE_STUDY = Path("case_studies/ecoli_o157_fitzgerald_2021")
TRIANGLE = CASE_STUDY / "results" / "triangle.tsv"
METADATA = CASE_STUDY / "results" / "metadata_with_lineage.tsv"
OUT_STEM = Path("paper/figures/supplementary/fig_s9_ecoli_o157_breakpoints")

# Colorblind-safe lineage palette.
LINEAGE_ORDER = ["I/II", "II", "Ia", "Ic"]
LINEAGE_COLORS = {
    "I/II": COLORS["orange"],
    "II": COLORS["sky_blue"],
    "Ia": COLORS["bluish_green"],
    "Ic": COLORS["vermillion"],
}

# Host category palette.
HOST_ORDER = ["bovine", "human", "other/unknown"]
HOST_COLORS = {
    "bovine": COLORS["blue"],
    "human": COLORS["vermillion"],
    "other/unknown": COLORS["grey"],
}

# Marker and plot settings.
MARKER_SIZE = 10
MARKER_ALPHA = 0.65
EDGE_COLOR = "white"
EDGE_WIDTH = 0.3


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def normalize_accession(acc: str | float) -> str | None:
    """Strip version suffix and NZ_ prefix to match metadata accessions."""
    if pd.isna(acc):
        return None
    acc = str(acc).split(".")[0]
    if acc.startswith("NZ_"):
        acc = acc[3:]
    return acc


def load_data() -> pd.DataFrame:
    """Load triangle results and merge query metadata for lineage/host."""
    tri = pd.read_csv(TRIANGLE, sep="\t")
    meta = pd.read_csv(METADATA, sep="\t")

    # Normalize accessions.
    tri["query"] = tri["query"].apply(normalize_accession)
    tri["reference"] = tri["reference"].apply(normalize_accession)
    meta = meta.rename(columns={"nucleotide_acc": "query"})
    meta["query"] = meta["query"].apply(normalize_accession)

    required_tri = {"query", "reference", "ani", "breakpoint_count"}
    missing_tri = required_tri - set(tri.columns)
    if missing_tri:
        raise ValueError(f"triangle.tsv missing columns: {missing_tri}")

    required_meta = {"query", "assigned_lineage", "host_category"}
    missing_meta = required_meta - set(meta.columns)
    if missing_meta:
        raise ValueError(f"metadata_with_lineage.tsv missing columns: {missing_meta}")

    df = tri.merge(
        meta[["query", "assigned_lineage", "host_category"]],
        on="query",
        how="left",
    )

    # Convert fractional ANI (0-1) if necessary; file appears to already be percent.
    if df["ani"].max() <= 1.0:
        df["ani_pct"] = df["ani"] * 100.0
    else:
        df["ani_pct"] = df["ani"].astype(float)

    # Validate caption claims.
    n_pairs = len(df)
    min_ani = df["ani_pct"].min()
    max_bp = df["breakpoint_count"].max()
    min_bp = df["breakpoint_count"].min()
    print(
        f"Loaded {n_pairs} pairs from {TRIANGLE}; "
        f"ANI range {min_ani:.4f}%–{df['ani_pct'].max():.4f}%; "
        f"breakpoint range {min_bp}–{max_bp}"
    )

    unmatched = df["assigned_lineage"].isna().sum()
    if unmatched:
        print(f"Warning: {unmatched} pairs lack query lineage metadata.")

    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _scatter_by_category(
    ax: plt.Axes,
    df: pd.DataFrame,
    category_col: str,
    category_order: list[str],
    color_map: dict[str, str],
    title: str,
) -> None:
    """Scatter ANI vs breakpoint_count colored by a categorical column."""
    x_col = "log10_breakpoint_count"
    plot_df = df.copy()
    plot_df[x_col] = np.log10(plot_df["breakpoint_count"])

    for cat in category_order:
        sub = plot_df[plot_df[category_col] == cat]
        if sub.empty:
            continue
        ax.scatter(
            sub[x_col],
            sub["ani_pct"],
            c=color_map[cat],
            s=MARKER_SIZE,
            alpha=MARKER_ALPHA,
            edgecolors=EDGE_COLOR,
            linewidths=EDGE_WIDTH,
            label=cat,
            zorder=2,
        )

    ax.set_xlabel("Breakpoint count")
    ax.set_ylabel("ANI (%)")
    ax.set_title(title, fontsize=8)

    # Log-scaled x ticks covering the observed range (~170–1100).
    x_ticks = [200, 300, 400, 500, 700, 1000]
    ax.set_xticks(np.log10(x_ticks))
    ax.set_xticklabels([str(t) for t in x_ticks])
    ax.set_xlim(np.log10(150), np.log10(1300))

    # Tight y-axis around the high-ANI cloud.
    ax.set_ylim(99.84, 100.02)
    ax.set_yticks([99.90, 99.95, 100.00])
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))

    ax.legend(
        title=category_col.replace("_", " ").title(),
        loc="lower right",
        frameon=False,
        handletextpad=0.2,
    )


def plot_fig_s9(df: pd.DataFrame) -> None:
    """Create and save the 1x2 supplementary figure."""
    set_publication_style()

    fig, axes = plt.subplots(
        1, 2, figsize=figure_size(17.5, aspect=0.42), sharey=True
    )
    fig.subplots_adjust(wspace=0.22)

    _scatter_by_category(
        axes[0],
        df,
        "assigned_lineage",
        LINEAGE_ORDER,
        LINEAGE_COLORS,
        "By lineage",
    )
    label_panel(axes[0], "a")

    _scatter_by_category(
        axes[1],
        df,
        "host_category",
        HOST_ORDER,
        HOST_COLORS,
        "By host category",
    )
    label_panel(axes[1], "b")

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, OUT_STEM, formats=("png", "pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    df = load_data()
    plot_fig_s9(df)


if __name__ == "__main__":
    main()
