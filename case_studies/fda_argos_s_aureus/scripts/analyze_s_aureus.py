#!/usr/bin/env python3
"""Correlate FDA-ARGOS S. aureus metadata with syn2bani triangle output."""

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

WORK_DIR = Path("/Volumes/MoneyCat/Data/fda_argos_staphylococcus_aureus")
GENOMES_DIR = WORK_DIR / "genomes"
FIGURES_DIR = WORK_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

MIN_CATEGORY_COUNT = 5
HIGH_ANI_THRESHOLD = 99.9


def parse_fasta_headers(fasta_path):
    seq_ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                # header format: >seqid description
                seq_id = line[1:].split()[0]
                seq_ids.append(seq_id)
    return seq_ids


def build_seq_to_assembly_map(genomes_dir):
    mapping = {}
    for fna in sorted(genomes_dir.glob("*.fna")):
        acc = fna.stem
        for seq_id in parse_fasta_headers(fna):
            mapping[seq_id] = acc
    return mapping


def group_categories(series, min_count=MIN_CATEGORY_COUNT):
    counts = series.value_counts()
    keep = set(counts[counts >= min_count].index)
    return series.apply(lambda x: x if x in keep else "other")


def load_triangle(triangle_path):
    df = pd.read_csv(triangle_path, sep="\t", low_memory=False)
    # Drop duplicated 'flag' column if present
    if isinstance(df.columns, pd.Index):
        cols = list(df.columns)
        seen = set()
        keep_cols = []
        for c in cols:
            if c in seen:
                keep_cols.append(f"{c}_dup")
            else:
                keep_cols.append(c)
                seen.add(c)
        df.columns = keep_cols
    return df


def main():
    triangle_path = WORK_DIR / "triangle.tsv"
    metadata_path = WORK_DIR / "assembly_metadata.tsv"

    if not triangle_path.exists():
        raise FileNotFoundError(triangle_path)
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)

    seq_map = build_seq_to_assembly_map(GENOMES_DIR)
    meta = pd.read_csv(metadata_path, sep="\t")
    meta["assembly_acc"] = meta["assembly_acc"].astype(str)

    tri = load_triangle(triangle_path)
    tri["query_acc"] = tri["query"].map(seq_map)
    tri["ref_acc"] = tri["reference"].map(seq_map)

    missing_query = tri["query_acc"].isna().sum()
    missing_ref = tri["ref_acc"].isna().sum()

    tri = tri.dropna(subset=["query_acc", "ref_acc"]).copy()
    tri = tri.merge(meta, left_on="query_acc", right_on="assembly_acc", how="left", suffixes=("", "_query"))
    tri = tri.merge(meta, left_on="ref_acc", right_on="assembly_acc", how="left", suffixes=("", "_ref"))

    # Use query metadata as the representative label for the pair
    tri["country_grouped"] = group_categories(tri["country"].fillna("unknown"))
    tri["source_grouped"] = group_categories(tri["isolation_source"].fillna("unknown"))

    # Basic stats
    n_pairs = len(tri)
    mean_ani = tri["ani"].mean()
    mean_breakpoints = tri["breakpoint_count"].mean()
    median_breakpoints = tri["breakpoint_count"].median()

    # Breakpoint burden by metadata
    burden_country = (
        tri.groupby("country_grouped")["breakpoint_count"]
        .agg(["count", "mean", "median", "std"])
        .sort_values("mean", ascending=False)
    )
    burden_source = (
        tri.groupby("source_grouped")["breakpoint_count"]
        .agg(["count", "mean", "median", "std"])
        .sort_values("mean", ascending=False)
    )

    # Top 20 high-ANI / high-breakpoint pairs
    high_ani = tri[tri["ani"] >= HIGH_ANI_THRESHOLD]
    top_pairs = high_ani.nlargest(20, "breakpoint_count")[
        [
            "query",
            "reference",
            "ani",
            "breakpoint_count",
            "af_query",
            "af_reference",
            "strain",
            "country",
            "isolation_source",
            "collection_date",
            "strain_ref",
            "country_ref",
            "isolation_source_ref",
            "collection_date_ref",
        ]
    ]

    # Plots
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=tri,
        x="breakpoint_count",
        y="ani",
        hue="country_grouped",
        alpha=0.7,
        s=40,
    )
    plt.title("ANI vs breakpoint count colored by country")
    plt.xlabel("Breakpoint count")
    plt.ylabel("ANI (%)")
    plt.legend(title="Country", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "ani_vs_breakpoints_country.png", dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=tri,
        x="breakpoint_count",
        y="ani",
        hue="source_grouped",
        alpha=0.7,
        s=40,
    )
    plt.title("ANI vs breakpoint count colored by isolation source")
    plt.xlabel("Breakpoint count")
    plt.ylabel("ANI (%)")
    plt.legend(title="Isolation source", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "ani_vs_breakpoints_source.png", dpi=300)
    plt.close()

    # Distribution plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sns.boxplot(data=tri, x="country_grouped", y="breakpoint_count", ax=axes[0])
    axes[0].set_title("Breakpoint count by country")
    axes[0].tick_params(axis="x", rotation=45)
    sns.boxplot(data=tri, x="source_grouped", y="breakpoint_count", ax=axes[1])
    axes[1].set_title("Breakpoint count by isolation source")
    axes[1].tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "breakpoint_burden_boxplots.png", dpi=300)
    plt.close()

    # Report
    def fmt_df(df, float_cols=None):
        df = df.copy()
        if float_cols is None:
            float_cols = [c for c in df.columns if pd.api.types.is_float_dtype(df[c])]
        for c in float_cols:
            df[c] = df[c].round(2)
        lines = ["| " + " | ".join(map(str, df.columns)) + " |"]
        lines.append("| " + " | ".join(["---"] * len(df.columns)) + " |")
        for _, row in df.iterrows():
            lines.append("| " + " | ".join(map(str, row)) + " |")
        return "\n".join(lines)

    burden_country = burden_country.round(2)
    burden_source = burden_source.round(2)
    top_pairs = top_pairs.round({c: 4 for c in top_pairs.columns if pd.api.types.is_float_dtype(top_pairs[c])})

    report_lines = [
        "# FDA-ARGOS *Staphylococcus aureus* analysis summary",
        "",
        f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d')}",
        f"**Genomes:** {len(meta)}",
        f"**Refined pairs:** {n_pairs}",
        f"**Sequence-ID mapping failures:** query={missing_query}, reference={missing_ref}",
        "",
        "## Pairwise statistics",
        "",
        f"- Mean ANI: {mean_ani:.4f}",
        f"- Mean breakpoint count: {mean_breakpoints:.2f}",
        f"- Median breakpoint count: {median_breakpoints:.2f}",
        "",
        "## Breakpoint burden by country",
        "",
        fmt_df(burden_country.reset_index()),
        "",
        "## Breakpoint burden by isolation source",
        "",
        fmt_df(burden_source.reset_index()),
        "",
        f"## Top 20 high-ANI (≥{HIGH_ANI_THRESHOLD}) / high-breakpoint pairs",
        "",
        fmt_df(top_pairs),
        "",
        "## Figures",
        "",
        "- `figures/ani_vs_breakpoints_country.png`",
        "- `figures/ani_vs_breakpoints_source.png`",
        "- `figures/breakpoint_burden_boxplots.png`",
        "",
    ]

    report_path = WORK_DIR / "summary_report.md"
    with open(report_path, "w") as fh:
        fh.write("\n".join(report_lines))

    print(f"Report saved to {report_path}")
    print(f"Figures saved to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
