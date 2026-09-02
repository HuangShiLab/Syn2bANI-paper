#!/usr/bin/env python3
"""Analyze the Syn2bANI triangle output for the Fitzgerald et al. 2021 E. coli O157:H7 dataset."""
import csv
import re
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

WORKDIR = Path(__file__).resolve().parent


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
    with (WORKDIR / "metadata.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [row for row in reader]
    df = pd.DataFrame(rows)
    df["lineage"] = df["serotype"].where(df["serotype"].astype(str).str.strip() != "", df["status"])
    return df


def main():
    tri = load_triangle()
    meta = load_metadata()

    # Attach metadata for query and reference; color by query lineage (proxy)
    tri = tri.merge(
        meta[["nucleotide_acc", "lineage", "status", "strain", "country", "assembly_acc"]].rename(
            columns={"nucleotide_acc": "query", "lineage": "q_lineage", "status": "q_status",
                     "strain": "q_strain", "country": "q_country", "assembly_acc": "q_assembly"}
        ),
        on="query", how="left"
    )
    tri = tri.merge(
        meta[["nucleotide_acc", "assembly_acc", "strain"]].rename(
            columns={"nucleotide_acc": "reference", "assembly_acc": "r_assembly", "strain": "r_strain"}
        ),
        on="reference", how="left"
    )

    # Summary statistics
    total_pairs = len(tri)
    ani_threshold = 95.0
    bp_threshold = 30
    high_ani_high_bp = tri[(tri["ani"] >= ani_threshold) & (tri["breakpoint_count"] >= bp_threshold)]
    n_high = len(high_ani_high_bp)
    top20 = high_ani_high_bp.sort_values(["breakpoint_count", "ani"], ascending=[False, False]).head(20)

    summary_lines = [
        "# E. coli O157:H7 (Fitzgerald et al. 2021) – Syn2bANI summary",
        "",
        f"- Genomes analyzed: {meta['assembly_acc'].nunique()}",
        f"- Total non-self pairs: {total_pairs}",
        f"- ANI range: {tri['ani'].min():.4f} – {tri['ani'].max():.4f}",
        f"- Breakpoint count range: {tri['breakpoint_count'].min()} – {tri['breakpoint_count'].max()}",
        f"- Pairs with ANI ≥ {ani_threshold}% AND breakpoints ≥ {bp_threshold}: {n_high} ({100*n_high/total_pairs:.2f}%)",
        "",
        "## Top 20 high-ANI/high-breakpoint pairs",
        "",
        "| query | reference | ANI | breakpoints | anchor_adjacency | q_lineage |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in top20.iterrows():
        summary_lines.append(
            f"| {r['query']} | {r['reference']} | {r['ani']:.4f} | {r['breakpoint_count']} | {r['anchor_adjacency']:.4f} | {r['q_lineage']} |"
        )
    summary_lines.append("")

    # Lineage/status composition
    summary_lines.append("## Color/group counts")
    summary_lines.append("")
    for group, cnt in Counter(tri["q_lineage"].fillna("unknown")).most_common():
        summary_lines.append(f"- {group}: {cnt} query observations")
    summary_lines.append("")

    (WORKDIR / "summary_report.md").write_text("\n".join(summary_lines))
    print(f"Wrote {WORKDIR / 'summary_report.md'}")

    # Plot 1: ANI vs breakpoint_count, colored by lineage/status
    fig, ax = plt.subplots(figsize=(9, 6))
    groups = tri["q_lineage"].fillna("unknown").unique()
    cmap = plt.cm.tab10
    for i, grp in enumerate(sorted(groups)):
        sub = tri[tri["q_lineage"].fillna("unknown") == grp]
        ax.scatter(sub["breakpoint_count"], sub["ani"], s=25, alpha=0.7,
                   label=grp, color=cmap(i % 10))
    ax.axhline(ani_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(bp_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("breakpoint_count")
    ax.set_ylabel("ANI (%)")
    ax.set_title("ANI vs breakpoint_count (colored by lineage/status proxy)")
    ax.legend(title="lineage/status", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    fig.tight_layout()
    fig.savefig(WORKDIR / "figures" / "ani_vs_breakpoints.png", dpi=300)
    plt.close(fig)
    print("Wrote figures/ani_vs_breakpoints.png")

    # Plot 2: breakpoint distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(tri["breakpoint_count"], bins=50, edgecolor="black")
    ax.axvline(bp_threshold, color="red", linestyle="--", linewidth=0.8, label=f"breakpoints ≥ {bp_threshold}")
    ax.set_xlabel("breakpoint_count")
    ax.set_ylabel("number of pairs")
    ax.set_title("Distribution of breakpoint counts")
    ax.legend()
    fig.tight_layout()
    fig.savefig(WORKDIR / "figures" / "breakpoint_distribution.png", dpi=300)
    plt.close(fig)
    print("Wrote figures/breakpoint_distribution.png")


if __name__ == "__main__":
    main()
