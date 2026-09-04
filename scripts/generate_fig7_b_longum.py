#!/usr/bin/env python3
"""Generate Figure 7: B. longum abfA locus chain coverage.

Panel (a): Syn2bANI chain coverage across the abfA locus on reference
FSHHK16M1_ctg_contig10, stratified by curated abfA cluster status.
Panel (b): Per-strain coverage of the cluster periphery (17.2-37.1 kb),
with Mann-Whitney U p-value.

Inputs
------
results/b_longum_abfA/metadata.tsv
    accession, abfA_status (complete/deleted), hypba_status, phenotype
results/b_longum_abfA/paf_vs_ref/<accession>.paf
    PAF alignments of each genome vs the abfA+ reference FSHHK16M1_ctg.
results/b_longum_abfA/abfA_region_coverage.tsv
    Per-strain core / periphery / downstream coverage fractions.

Outputs
-------
paper/figures/main/fig7_b_longum_abfa.png
paper/figures/main/fig7_b_longum_abfa.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from plot_style import COLORS, figure_size, label_panel, save_figure, set_publication_style

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
WORK = Path("results/b_longum_abfA")
OUT_STEM = Path("paper/figures/main/fig7_b_longum_abfa")

REFERENCE_CONTIG = "FSHHK16M1_ctg_contig10"
CONTIG_LEN = 113_170

# Locus coordinates on the reference contig (bp).
LOCUS_S, LOCUS_E = 8_546, 37_075
CORE_S, CORE_E = 8_840, 17_163
PERI_S, PERI_E = 17_163, 37_075

BIN_SIZE = 1_000  # bp

# Group labels and display properties.
GROUP_ORDER = ["complete", "deleted"]
GROUP_COLORS = {
    "complete": COLORS["sky_blue"],
    "deleted": COLORS["vermillion"],
}


def load_metadata(path: Path) -> pd.DataFrame:
    """Load metadata and keep only complete/deleted abfA statuses."""
    meta = pd.read_csv(path, sep="\t")
    required = {"accession", "abfA_status"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"metadata.tsv missing columns: {missing}")
    return meta[meta["abfA_status"].isin(GROUP_ORDER)].copy()


def load_paf_coverage(paf_dir: Path, meta: pd.DataFrame) -> pd.DataFrame:
    """Compute per-genome 1-kb binned coverage on the reference contig from PAF.

    Returns a DataFrame with one row per genome and one column per bin.
    """
    n_bins = (CONTIG_LEN + BIN_SIZE - 1) // BIN_SIZE
    records = []
    for paf in sorted(paf_dir.glob("*.paf")):
        acc = paf.stem
        if acc not in set(meta["accession"]):
            continue
        cov = np.zeros(CONTIG_LEN, dtype=bool)
        with open(paf) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12 or f[5] != REFERENCE_CONTIG:
                    continue
                ts, te = int(f[7]), int(f[8])
                if te > ts:
                    cov[ts:te] = True
        binned = np.array(
            [cov[i * BIN_SIZE : (i + 1) * BIN_SIZE].mean() for i in range(n_bins)],
            dtype=float,
        )
        records.append((acc, *binned.tolist()))

    cols = ["accession"] + [f"bin_{i}" for i in range(n_bins)]
    return pd.DataFrame(records, columns=cols)


def panel_a(ax: plt.Axes, cov_df: pd.DataFrame) -> None:
    """Plot group-wise chain coverage across the abfA locus."""
    bin_cols = [c for c in cov_df.columns if c.startswith("bin_")]
    bin_edges = np.arange(len(bin_cols) + 1) * BIN_SIZE / 1_000  # kb

    # Background shading for functional regions.
    ax.axvspan(
        CORE_S / 1_000, CORE_E / 1_000, color=COLORS["grey"], alpha=0.15, zorder=0
    )
    ax.axvspan(
        PERI_S / 1_000, PERI_E / 1_000, color=COLORS["orange"], alpha=0.12, zorder=0
    )
    # Locus boundary lines.
    ax.axvline(LOCUS_S / 1_000, color=COLORS["grey"], linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(CORE_E / 1_000, color=COLORS["grey"], linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(LOCUS_E / 1_000, color=COLORS["grey"], linestyle="--", linewidth=0.8, zorder=1)

    for status in GROUP_ORDER:
        sub = cov_df[cov_df["abfA_status"] == status]
        frac = sub[bin_cols].mean(axis=0).to_numpy() * 100
        n = len(sub)
        label = f"abfA cluster {status} ($n$={n})"
        ax.step(
            bin_edges[:-1],
            frac,
            where="post",
            color=GROUP_COLORS[status],
            linewidth=1.4,
            label=label,
            zorder=2,
        )

    # Region annotation (above the plot area, in axes coordinates).
    core_ax_x = ((CORE_S + CORE_E) / 2 / 1_000) / (CONTIG_LEN / 1_000)
    peri_ax_x = ((PERI_S + PERI_E) / 2 / 1_000) / (CONTIG_LEN / 1_000)
    ax.text(
        core_ax_x,
        1.07,
        "core genes",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6,
        color="#333333",
    )
    ax.text(
        peri_ax_x,
        1.01,
        "periphery",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6,
        color="#333333",
    )

    ax.set_xlim(0, CONTIG_LEN / 1_000)
    ax.set_ylim(0, 125)
    ax.set_ylabel("Genomes with chained anchors (%)")
    ax.set_xticklabels([])  # shared x-label added by main()
    ax.legend(loc="lower left", frameon=True, edgecolor="none", facecolor="white")
    label_panel(ax, "a", x=-0.09, y=1.04)


def panel_b(ax: plt.Axes, region_df: pd.DataFrame) -> None:
    """Boxplot of periphery coverage by abfA status with MWU p-value."""
    data = [region_df[region_df["abfA_status"] == g]["periphery"].to_numpy() for g in GROUP_ORDER]
    positions = [1, 2]
    n_complete = len(data[0])
    n_deleted = len(data[1])

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#333333", "linewidth": 1.2},
        whiskerprops={"color": "#333333", "linewidth": 0.8},
        capprops={"color": "#333333", "linewidth": 0.8},
    )
    for patch, status in zip(bp["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[status])
        patch.set_edgecolor("#333333")
        patch.set_alpha(0.85)
        patch.set_linewidth(0.8)

    # Jittered individual points.
    rng = np.random.default_rng(42)
    for pos, vals, status in zip(positions, data, GROUP_ORDER):
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(
            np.full(len(vals), pos) + jitter,
            vals,
            color=GROUP_COLORS[status],
            edgecolor="none",
            alpha=0.45,
            s=8,
            zorder=3,
        )

    # One-sided MWU: complete > deleted (directional hypothesis in the text).
    _, pval = stats.mannwhitneyu(
        data[0], data[1], alternative="greater", method="auto"
    )
    exponent = int(np.floor(np.log10(pval)))
    mantissa = pval / (10 ** exponent)
    p_text = f"MWU $p$ = {mantissa:.1f} × 10$^{{{exponent}}}$"

    y_max = max(np.max(d) for d in data)
    ax.annotate(
        p_text,
        xy=(1.5, y_max),
        xytext=(1.5, y_max + 0.08),
        ha="center",
        va="bottom",
        fontsize=8,
    )

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"complete\n($n$={n_complete})", f"deleted\n($n$={n_deleted})"]
    )
    ax.set_ylabel("Periphery coverage")
    ax.set_ylim(-0.05, 1.25)
    ax.set_xlim(0.4, 2.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    label_panel(ax, "b", x=-0.09, y=1.04)


def main() -> None:
    set_publication_style()

    meta = load_metadata(WORK / "metadata.tsv")
    cov_df = load_paf_coverage(WORK / "paf_vs_ref", meta)
    cov_df = cov_df.merge(meta, on="accession", how="left")

    region_df = pd.read_csv(WORK / "abfA_region_coverage.tsv", sep="\t")
    region_df = region_df[region_df["abfA_status"].isin(GROUP_ORDER)].copy()

    # Use the same sample set for both panels.
    used_accs = set(cov_df["accession"])
    region_df = region_df[region_df["accession"].isin(used_accs)].copy()

    fig, axes = plt.subplots(
        2, 1, figsize=figure_size(8.5, aspect=1.15), gridspec_kw={"hspace": 0.28}
    )

    panel_a(axes[0], cov_df)
    panel_b(axes[1], region_df)

    fig.text(
        0.55,
        0.01,
        "Position on FSHHK16M1 contig10 (kb)",
        ha="center",
        fontsize=9,
    )

    fig.subplots_adjust(bottom=0.16)
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, OUT_STEM, formats=("png", "pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
