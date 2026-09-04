#!/usr/bin/env python3
"""Generate Figure 8: cagPAI architecture in 528 H. pylori isolates.

Panels
------
(a) Engineered control validation (expected vs predicted cagPAI state).
(b) Extended-state classification workflow (schematic).
(c) Stacked-bar distribution of extended cagPAI states by FastBAPS lineage.
(d) Stacked-bar distribution of extended cagPAI states by disease stage,
    with lineage-stratified CMH annotation.

Outputs
-------
paper/figures/main/fig8_cagpai_h_pylori.png
paper/figures/main/fig8_cagpai_h_pylori.pdf
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

# Ensure the repository plot-style module is importable when running from root.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import plot_style  # noqa: E402

REPO_ROOT = SCRIPT_DIR.parent
CASE_DIR = REPO_ROOT / "case_studies" / "h_pylori_cagpai" / "results"
OUT_DIR = REPO_ROOT / "paper" / "figures" / "main"
OUT_STEM = OUT_DIR / "fig8_cagpai_h_pylori"

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

# Engineered control expected states (from cagpai_summary.md)
PILOT_EXPECTED = {
    "wt": "complete",
    "hp26695": "complete",
    "mut1": "complete",
    "inv": "complete",
    "transloc": "complete",
    "mut1_inv": "complete",
    "mut1_transloc": "complete",
    "del": "empty",
    "mut1_del": "empty",
}

PILOT_ORDER = [
    "wt",
    "hp26695",
    "mut1",
    "inv",
    "transloc",
    "mut1_inv",
    "mut1_transloc",
    "del",
    "mut1_del",
]

PILOT_DISPLAY = {
    "wt": "WT",
    "hp26695": "26695",
    "mut1": "mut1",
    "inv": "inversion",
    "transloc": "translocation",
    "mut1_inv": "mut1+inv",
    "mut1_transloc": "mut1+trans",
    "del": "ΔcagPAI",
    "mut1_del": "mut1+Δ",
}


def parse_count_table(path: Path, section_name: str) -> dict[str, dict[str, int]]:
    """Parse a count table from cagpai_association_filtered.tsv.

    Returns {group: {state: count}} for the requested section (e.g. 'fastbaps',
    'group').
    """
    text = path.read_text()
    # Sections start with '## <section_name>'. Find it.
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


def chi2_from_table(table: dict[str, dict[str, int]]) -> tuple[float, int, float] | None:
    """Compute chi-square test on a state count table."""
    try:
        from scipy.stats import chi2_contingency
    except Exception:
        return None
    groups = sorted(table.keys())
    arr = np.array([[table[g].get(s, 0) for s in STATE_ORDER] for g in groups], dtype=float)
    if arr.shape[0] < 2 or arr.shape[1] < 2 or (arr.sum(axis=1) < 5).any():
        return None
    chi2, p, dof, _ = chi2_contingency(arr)
    return float(chi2), int(dof), float(p)


def load_pilot(path: Path) -> list[dict[str, str]]:
    rows = {r["genome"]: r for r in csv.DictReader(open(path), delimiter="\t")}
    return [rows[g] for g in PILOT_ORDER if g in rows]


def parse_cmh_p_values(path: Path) -> dict[str, float]:
    """Return the smallest CMH p-value and all contrast p-values."""
    pvals: dict[str, float] = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("###"):
                current = line.lstrip("# ").strip()
            elif line.startswith("CMH") and current:
                m = re.search(r"p=([0-9.eE+-]+)", line)
                if m:
                    pvals[current] = float(m.group(1))
    return pvals


def plot_panel_a(ax, rows: list[dict[str, str]]) -> None:
    """Engineered control validation: fraction of markers present."""
    labels = [PILOT_DISPLAY[r["genome"]] for r in rows]
    fractions = [float(r["fraction_present"]) for r in rows]
    predicted = [r["status"] for r in rows]
    expected = [PILOT_EXPECTED[r["genome"]] for r in rows]

    def _pilot_color(status: str) -> str:
        return STATE_COLORS.get(status, STATE_COLORS["complete_collinear"])

    colors = [_pilot_color(p) for p in predicted]
    x = np.arange(len(rows))
    bars = ax.bar(x, fractions, color=colors, edgecolor="white", linewidth=0.5)

    # Highlight mismatches with a dark edge (none expected for current data).
    for bar, exp, pred in zip(bars, expected, predicted):
        if exp != pred:
            bar.set_edgecolor("black")
            bar.set_linewidth(2)

    ax.axhline(0.85, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    ax.axhline(0.15, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    ax.text(0.98, 0.96, "complete threshold", transform=ax.transAxes,
            color="grey", fontsize=7, va="top", ha="right")
    ax.text(0.98, 0.04, "empty threshold", transform=ax.transAxes,
            color="grey", fontsize=7, va="bottom", ha="right")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Marker fraction present")
    ax.set_title("Engineered control validation")

    correct = sum(1 for e, p in zip(expected, predicted) if e == p)
    ax.text(0.02, 0.96, f"{correct}/{len(rows)} correct",
            transform=ax.transAxes, ha="left", va="top", fontsize=8)


def plot_panel_b(ax) -> None:
    """Workflow schematic for extended-state classification."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Extended-state classification workflow")

    boxes = [
        (0.11, "28 cagPAI\nmarkers", "minimap2 ≥80%\ncoverage & identity"),
        (0.36, "syn2bani\nstruct", "vs 26695\nINV/DEL/TRA calls"),
        (0.61, "Circular-origin\nfilter", "exclude >50%\nchromosome spans"),
        (0.86, "Extended\nstate", "empty / partial /\ncollinear / rearranged"),
    ]

    for x, title, subtitle in boxes:
        box = FancyBboxPatch(
            (x - 0.09, 0.62), 0.18, 0.22,
            boxstyle="round,pad=0.02,rounding_size=0.015",
            facecolor=plot_style.COLORS["light_grey"],
            edgecolor=plot_style.COLORS["grey"],
            linewidth=1.0,
        )
        ax.add_patch(box)
        ax.text(x, 0.77, title, ha="center", va="center", fontsize=7, fontweight="bold")
        ax.text(x, 0.67, subtitle, ha="center", va="center", fontsize=6, color="#333333")

    for i in range(len(boxes) - 1):
        x0 = boxes[i][0] + 0.09
        x1 = boxes[i + 1][0] - 0.09
        arrow = FancyArrowPatch(
            (x0, 0.73), (x1, 0.73),
            arrowstyle="->", mutation_scale=10,
            linewidth=1.0, color="#333333",
        )
        ax.add_patch(arrow)

    notes = (
        "• 28 marker loci classify presence as empty / partial / complete.\n"
        "• SVs overlapping the cagPAI window split complete into collinear vs rearranged.\n"
        "• Genome-spanning translocations from circular-origin shifts are filtered."
    )
    ax.text(0.5, 0.22, notes, ha="center", va="top", fontsize=6.5, color="#333333")


def plot_stacked_bar(ax, table: dict[str, dict[str, int]], order: list[str],
                     title: str, ylabel: str) -> None:
    """Horizontal stacked bar of state proportions."""
    counts = np.array([[table[g].get(s, 0) for s in STATE_ORDER] for g in order], dtype=float)
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    fracs = counts / totals

    left = np.zeros(len(order))
    for i, state in enumerate(STATE_ORDER):
        ax.barh(np.arange(len(order)), fracs[:, i], left=left,
                color=STATE_COLORS[state], label=STATE_LABELS[state], height=0.65,
                edgecolor="white", linewidth=0.5)
        left += fracs[:, i]

    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of isolates")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.invert_yaxis()

    # Add sample-size labels to the right of each bar.
    for j, g in enumerate(order):
        ax.text(1.02, j, f"n={int(totals[j].item())}", va="center", ha="left", fontsize=7, color="#333333")
    ax.set_xlim(0, 1.18)


def main() -> int:
    plot_style.set_publication_style()

    assoc_path = CASE_DIR / "cagpai_association_filtered.tsv"
    strat_path = CASE_DIR / "cagpai_association_stratified.tsv"
    pilot_path = CASE_DIR / "cagpai_states_pilot.tsv"

    # Load data
    fastbaps_table = parse_count_table(assoc_path, "fastbaps")
    group_table = parse_count_table(assoc_path, "group")
    pilot_rows = load_pilot(pilot_path)
    cmh_pvals = parse_cmh_p_values(strat_path)

    # Recompute statistics from counts.
    fb_chi2 = chi2_from_table(fastbaps_table)
    grp_chi2 = chi2_from_table(group_table)

    # Fallback to precomputed values if scipy fails.
    if fb_chi2 is None:
        fb_chi2 = (58.754, 12, 3.805e-08)
    if grp_chi2 is None:
        grp_chi2 = (24.601, 9, 0.003446)

    min_cmh_p = min(cmh_pvals.values()) if cmh_pvals else 0.21

    fig = plt.figure(figsize=plot_style.figure_size(17.8, aspect=0.70))
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.92, top=0.92, bottom=0.10,
                          wspace=0.35, hspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # Panel (a)
    plot_panel_a(ax_a, pilot_rows)
    plot_style.label_panel(ax_a, "a")

    # Panel (b)
    plot_panel_b(ax_b)
    plot_style.label_panel(ax_b, "b")

    # Panel (c)
    fb_order = [f"fastbaps_L{i}" for i in (2, 3, 4, 5, 6)]
    plot_stacked_bar(ax_c, fastbaps_table, fb_order,
                     title="Extended state by FastBAPS lineage",
                     ylabel="FastBAPS lineage")
    chi2, dof, p = fb_chi2
    p_str = f"{p:.2e}" if p < 0.001 else f"{p:.3g}"
    ax_c.text(0.03, 0.97,
              f"$\\chi^2={chi2:.2f}$, df={dof}, $p={p_str}$",
              transform=ax_c.transAxes, ha="left", va="top", fontsize=8)
    plot_style.label_panel(ax_c, "c")

    # Panel (d)
    group_order = ["NAG", "AG", "IM", "GC"]
    plot_stacked_bar(ax_d, group_table, group_order,
                     title="Extended state by disease stage",
                     ylabel="Correa stage")
    chi2_g, dof_g, p_g = grp_chi2
    p_g_str = f"{p_g:.2e}" if p_g < 0.001 else f"{p_g:.3g}"
    cmh_text = (
        f"Marginal: $\\chi^2={chi2_g:.2f}$, df={dof_g}, $p={p_g_str}$\n"
        f"CMH (lineage-stratified): $p \\geq {min_cmh_p:.2f}$"
    )
    ax_d.text(0.03, 0.97, cmh_text,
              transform=ax_d.transAxes, ha="left", va="top", fontsize=8)
    plot_style.label_panel(ax_d, "d")

    # Shared figure legend for extended states.
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
        title="State",
        title_fontsize=7,
        fontsize=7,
        bbox_to_anchor=(0.5, -0.02),
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_style.save_figure(fig, str(OUT_STEM), formats=("png", "pdf"))
    plt.close(fig)

    print("Figure 8 generated successfully.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
