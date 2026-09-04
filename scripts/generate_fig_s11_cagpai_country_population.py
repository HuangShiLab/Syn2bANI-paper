#!/usr/bin/env python3
"""Generate Supplementary Figure S11: cagPAI extended state distribution.

Panels
------
(a) Stacked-bar distribution of extended cagPAI states by country of isolation.
(b) Stacked-bar distribution of extended cagPAI states by phylogenetic population.

Outputs
-------
paper/figures/supplementary/fig_s11_cagpai_country_population.png
paper/figures/supplementary/fig_s11_cagpai_country_population.pdf
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Make the repository plot-style module importable when running from root.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import plot_style  # noqa: E402

REPO_ROOT = SCRIPT_DIR.parent
CASE_DIR = REPO_ROOT / "case_studies" / "h_pylori_cagpai" / "results"
OUT_DIR = REPO_ROOT / "paper" / "figures" / "supplementary"
OUT_STEM = OUT_DIR / "fig_s11_cagpai_country_population"

ASSOC_FILE = CASE_DIR / "cagpai_association_filtered.tsv"

STATE_ORDER = ["empty", "partial", "complete_collinear", "complete_rearranged"]
STATE_LABELS = {
    "empty": "empty",
    "partial": "partial",
    "complete_collinear": "complete collinear",
    "complete_rearranged": "complete rearranged",
}

STATE_COLORS = {
    "empty": plot_style.COLORS["vermillion"],
    "partial": plot_style.COLORS["orange"],
    "complete_collinear": plot_style.COLORS["bluish_green"],
    "complete_rearranged": plot_style.COLORS["blue"],
}


def parse_count_table(path: Path, section_name: str) -> dict[str, dict[str, int]]:
    """Parse a contingency count table from `cagpai_association_filtered.tsv`.

    Returns
    -------
    dict[str, dict[str, int]]
        ``{group: {state: count}}`` for the requested section.
    """
    text = path.read_text()
    pattern = rf"^## {re.escape(section_name)}\n(.*?)\n(?:## |chi2=|\\Z)"
    match = re.search(pattern, text, re.MULTILINE | re.DOTALL)
    if not match:
        raise ValueError(f"Section '{section_name}' not found in {path}")

    lines = match.group(1).strip().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Section '{section_name}' has no data rows")

    header = lines[0].split("\t")
    if header[0] != "group":
        raise ValueError(f"Unexpected header in {section_name}: {header}")

    table: dict[str, dict[str, int]] = {}
    for line in lines[1:]:
        cols = line.split("\t")
        group = cols[0]
        table[group] = {}
        for state, cell in zip(STATE_ORDER, cols[1 : 1 + len(STATE_ORDER)]):
            m = re.match(r"^(\d+)", cell.strip())
            if not m:
                raise ValueError(f"Cannot parse count from '{cell}' in {section_name}")
            table[group][state] = int(m.group(1))
    return table


def plot_stacked_bar(ax, table: dict[str, dict[str, int]], order: list[str],
                     title: str, ylabel: str) -> None:
    """Plot a horizontal stacked bar of state proportions."""
    counts = np.array(
        [[table[g].get(s, 0) for s in STATE_ORDER] for g in order], dtype=float
    )
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    fracs = counts / totals

    left = np.zeros(len(order))
    for state in STATE_ORDER:
        ax.barh(
            np.arange(len(order)),
            fracs[:, STATE_ORDER.index(state)],
            left=left,
            color=STATE_COLORS[state],
            label=STATE_LABELS[state],
            height=0.65,
            edgecolor="white",
            linewidth=0.5,
        )
        left += fracs[:, STATE_ORDER.index(state)]

    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of isolates")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.invert_yaxis()

    # Sample-size labels to the right of each bar.
    for j, g in enumerate(order):
        ax.text(
            1.02,
            j,
            f"n={int(totals[j].item())}",
            va="center",
            ha="left",
            fontsize=7,
            color="#333333",
        )
    ax.set_xlim(0, 1.18)


def main() -> int:
    plot_style.set_publication_style()

    country_table = parse_count_table(ASSOC_FILE, "country")
    pop_table = parse_count_table(ASSOC_FILE, "phylogenetic_population")

    # Order categories by descending sample size so the panels are readable.
    country_order = sorted(
        country_table.keys(), key=lambda g: sum(country_table[g].values()), reverse=True
    )
    pop_order = sorted(
        pop_table.keys(), key=lambda g: sum(pop_table[g].values()), reverse=True
    )

    # Wide, moderately tall layout for many categorical rows.
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=plot_style.figure_size(17.8, aspect=0.85),
        gridspec_kw={"wspace": 0.45},
    )

    plot_stacked_bar(
        ax_a,
        country_table,
        country_order,
        title="Country of isolation",
        ylabel="Country",
    )
    plot_style.label_panel(ax_a, "a")

    plot_stacked_bar(
        ax_b,
        pop_table,
        pop_order,
        title="Phylogenetic population",
        ylabel="Population",
    )
    plot_style.label_panel(ax_b, "b")

    # Shared legend below the panels.
    legend_handles = [
        Patch(facecolor=STATE_COLORS[s], edgecolor="none", label=STATE_LABELS[s])
        for s in STATE_ORDER
    ]
    fig.legend(
        legend_handles,
        [h.get_label() for h in legend_handles],
        loc="lower center",
        ncol=4,
        frameon=False,
        title="Extended cagPAI state",
        title_fontsize=8,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_style.save_figure(fig, str(OUT_STEM), formats=("png", "pdf"))
    plt.close(fig)

    print("Supplementary Figure S11 generated successfully.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
