#!/usr/bin/env python3
"""Lineage-aware analysis of the Syn2bANI triangle output for Fitzgerald et al. 2021."""
import csv
import re
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd

WORKDIR = Path(__file__).resolve().parent

# Consistent lineage colors
LINEAGE_COLORS = {
    "Ia": "#1f77b4",
    "Ic": "#ff7f0e",
    "I/II": "#2ca02c",
    "II": "#d62728",
    "unassigned": "#7f7f7f",
}


def normalize_acc(s):
    """Map versioned RefSeq/GenBank accessions to the base nucleotide accession."""
    s = re.sub(r"^NZ_", "", s)
    return s.split(".")[0]


def load_triangle():
    rows = []
    with (WORKDIR / "triangle.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({
                "query": normalize_acc(row["query"]),
                "reference": normalize_acc(row["reference"]),
                "ani": float(row["ani"]),
                "breakpoint_count": int(row["breakpoint_count"]),
                "anchor_adjacency": float(row["anchor_adjacency"]),
                "af_query": float(row["af_query"]),
                "af_reference": float(row["af_reference"]),
            })
    return pd.DataFrame(rows)


def load_metadata():
    with (WORKDIR / "metadata_with_lineage.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [row for row in reader]
    df = pd.DataFrame(rows)
    # keep key columns
    return df[["nucleotide_acc", "strain", "assigned_lineage", "host_category", "isolation_source", "country"]]


def main():
    tri = load_triangle()
    meta = load_metadata()

    # Attach metadata for query and reference
    tri = tri.merge(
        meta.rename(columns={
            "nucleotide_acc": "query",
            "strain": "q_strain",
            "assigned_lineage": "q_lineage",
            "host_category": "q_host",
            "isolation_source": "q_source",
            "country": "q_country",
        }),
        on="query", how="left"
    )
    tri = tri.merge(
        meta.rename(columns={
            "nucleotide_acc": "reference",
            "strain": "r_strain",
            "assigned_lineage": "r_lineage",
            "host_category": "r_host",
        }),
        on="reference", how="left"
    )

    # Thresholds
    ani_threshold = 95.0
    bp_threshold = 30
    high_ani_high_bp = tri[(tri["ani"] >= ani_threshold) & (tri["breakpoint_count"] >= bp_threshold)]
    n_high = len(high_ani_high_bp)
    total_pairs = len(tri)

    # Per-lineage breakpoint stats (pairwise observations)
    lineage_stats = (
        tri.groupby("q_lineage")["breakpoint_count"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
        .sort_values("q_lineage")
    )

    # Host category breakpoint stats
    host_stats = (
        tri.groupby("q_host")["breakpoint_count"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
        .sort_values("mean", ascending=False)
    )

    # Top 20 high-ANI/high-breakpoint pairs
    top20 = high_ani_high_bp.sort_values(["breakpoint_count", "ani"], ascending=[False, False]).head(20)

    # --- Summary report ---
    summary_lines = [
        "# E. coli O157:H7 (Fitzgerald et al. 2021) – lineage-aware Syn2bANI summary",
        "",
        f"- Genomes analyzed: {meta['nucleotide_acc'].nunique()}",
        f"- Total non-self pairs: {total_pairs}",
        f"- ANI range: {tri['ani'].min():.4f} – {tri['ani'].max():.4f}",
        f"- Breakpoint count range: {tri['breakpoint_count'].min()} – {tri['breakpoint_count'].max()}",
        f"- Pairs with ANI ≥ {ani_threshold}% AND breakpoints ≥ {bp_threshold}: {n_high} ({100*n_high/total_pairs:.2f}%)",
        "",
        "## Lineage distribution",
        "",
    ]
    lineage_counts = Counter(meta["assigned_lineage"].fillna("unassigned"))
    for lineage, cnt in sorted(lineage_counts.items()):
        summary_lines.append(f"- {lineage}: {cnt} genomes")
    summary_lines.append("")

    summary_lines.append("## Per-lineage breakpoint statistics (pairwise observations)")
    summary_lines.append("")
    summary_lines.append("| lineage | pairs | mean | median | std | min | max |")
    summary_lines.append("|---|---|---|---|---|---|---|")
    for _, r in lineage_stats.iterrows():
        summary_lines.append(
            f"| {r['q_lineage']} | {int(r['count'])} | {r['mean']:.1f} | {r['median']:.1f} | {r['std']:.1f} | {int(r['min'])} | {int(r['max'])} |"
        )
    summary_lines.append("")

    summary_lines.append("## Breakpoint burden by host category")
    summary_lines.append("")
    summary_lines.append("| host_category | pairs | mean | median | std | min | max |")
    summary_lines.append("|---|---|---|---|---|---|---|")
    for _, r in host_stats.iterrows():
        summary_lines.append(
            f"| {r['q_host']} | {int(r['count'])} | {r['mean']:.1f} | {r['median']:.1f} | {r['std']:.1f} | {int(r['min'])} | {int(r['max'])} |"
        )
    summary_lines.append("")

    summary_lines.append("## Top 20 high-ANI/high-breakpoint pairs")
    summary_lines.append("")
    summary_lines.append("| query | q_lineage | reference | r_lineage | ANI | breakpoints | anchor_adjacency | q_host | r_host |")
    summary_lines.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in top20.iterrows():
        summary_lines.append(
            f"| {r['query']} | {r['q_lineage']} | {r['reference']} | {r['r_lineage']} | {r['ani']:.4f} | {r['breakpoint_count']} | {r['anchor_adjacency']:.4f} | {r['q_host']} | {r['r_host']} |"
        )
    summary_lines.append("")

    (WORKDIR / "summary_report_lineage.md").write_text("\n".join(summary_lines))
    print(f"Wrote {WORKDIR / 'summary_report_lineage.md'}")

    # --- Plot 1: ANI vs breakpoint_count, colored by assigned_lineage ---
    fig, ax = plt.subplots(figsize=(9, 6))
    for lineage in sorted(tri["q_lineage"].dropna().unique()):
        sub = tri[tri["q_lineage"] == lineage]
        ax.scatter(
            sub["breakpoint_count"], sub["ani"], s=25, alpha=0.7,
            label=lineage, color=LINEAGE_COLORS.get(lineage, "#7f7f7f")
        )
    ax.axhline(ani_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(bp_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("breakpoint_count")
    ax.set_ylabel("ANI (%)")
    ax.set_title("ANI vs breakpoint_count colored by query lineage")
    ax.legend(title="assigned_lineage", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    fig.tight_layout()
    fig.savefig(WORKDIR / "figures" / "ani_vs_breakpoints_lineage.png", dpi=300)
    plt.close(fig)
    print("Wrote figures/ani_vs_breakpoints_lineage.png")

    # --- Plot 2: breakpoint_count distribution by lineage (violin + swarm-like overlay) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    lineages = sorted(tri["q_lineage"].dropna().unique())
    data = [tri[tri["q_lineage"] == lin]["breakpoint_count"].values for lin in lineages]
    parts = ax.violinplot(data, positions=range(len(lineages)), showmeans=True, showmedians=True)
    for i, lin in enumerate(lineages):
        color = LINEAGE_COLORS.get(lin, "#7f7f7f")
        if i < len(parts["bodies"]):
            parts["bodies"][i].set_facecolor(color)
            parts["bodies"][i].set_alpha(0.5)
        # jittered scatter overlay
        y = tri[tri["q_lineage"] == lin]["breakpoint_count"]
        x = i + pd.Series(range(len(y))).apply(lambda v: (v % 7 - 3) * 0.04)
        ax.scatter(x, y, s=10, alpha=0.5, color=color)

    ax.set_xticks(range(len(lineages)))
    ax.set_xticklabels(lineages)
    ax.set_xlabel("assigned_lineage")
    ax.set_ylabel("breakpoint_count")
    ax.set_title("Distribution of breakpoint counts by lineage")
    fig.tight_layout()
    fig.savefig(WORKDIR / "figures" / "breakpoints_by_lineage.png", dpi=300)
    plt.close(fig)
    print("Wrote figures/breakpoints_by_lineage.png")

    # --- Plot 3: breakpoint_count distribution by host category ---
    fig, ax = plt.subplots(figsize=(8, 5))
    hosts = sorted(tri["q_host"].dropna().unique())
    data = [tri[tri["q_host"] == h]["breakpoint_count"].values for h in hosts]
    parts = ax.violinplot(data, positions=range(len(hosts)), showmeans=True, showmedians=True)
    host_colors = {"human": "#9467bd", "bovine": "#8c564b", "other/unknown": "#7f7f7f"}
    for i, host in enumerate(hosts):
        color = host_colors.get(host, "#7f7f7f")
        if i < len(parts["bodies"]):
            parts["bodies"][i].set_facecolor(color)
            parts["bodies"][i].set_alpha(0.5)
        y = tri[tri["q_host"] == host]["breakpoint_count"]
        x = i + pd.Series(range(len(y))).apply(lambda v: (v % 7 - 3) * 0.04)
        ax.scatter(x, y, s=10, alpha=0.5, color=color)

    ax.set_xticks(range(len(hosts)))
    ax.set_xticklabels(hosts)
    ax.set_xlabel("host_category")
    ax.set_ylabel("breakpoint_count")
    ax.set_title("Distribution of breakpoint counts by host category")
    fig.tight_layout()
    fig.savefig(WORKDIR / "figures" / "breakpoints_by_host.png", dpi=300)
    plt.close(fig)
    print("Wrote figures/breakpoints_by_host.png")


if __name__ == "__main__":
    main()
