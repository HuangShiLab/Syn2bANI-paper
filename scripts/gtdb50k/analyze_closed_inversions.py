#!/usr/bin/env python3
"""Summarize the closed-genome all-vs-all Syn2b structural run.

Inputs:
    --pairs      closed_inversion_pairs.tsv (the 371 near-closed pairs used as seed)
    --syn2b      syn2b_inverted_fraction_closed.tsv (all-vs-all on 680 genomes)
    --out-tsv    output ranking TSV
    --out-md     output markdown report
"""
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", required=True)
    p.add_argument("--syn2b", required=True)
    p.add_argument("--genomes", help="genome metadata TSV with acc/species")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    syn = pd.read_csv(args.syn2b, sep="\t", low_memory=False)
    syn = syn[syn["status"] == "ok"].copy()

    # Basic distribution
    n_total = len(syn)
    inv = syn["syn2b_raw_inverted_fraction"].astype(float)
    bp = syn["syn2b_breakpoints"].astype(int)

    lines = [
        "# Closed-genome all-vs-all structural summary\n",
        f"**Date:** generated from `{Path(args.syn2b).name}`\n",
        "\n## 1. Run overview\n",
        f"- Total pairs: {n_total:,}\n",
        f"- Reverse-direction columns included: {'syn2b_rev_raw_inverted_fraction' in syn.columns}\n",
        f"- Pairs with non-zero raw_inverted_fraction: {(inv > 0).sum():,} ({(inv > 0).mean()*100:.1f}%)\n",
        f"- Pairs with non-zero breakpoints: {(bp > 0).sum():,} ({(bp > 0).mean()*100:.1f}%)\n",
        "\n## 2. raw_inverted_fraction distribution\n",
        "| quantile | value |\n|---|---:|\n",
    ]
    for q in [0.00, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00]:
        lines.append(f"| {q*100:.0f}% | {inv.quantile(q):.4f} |\n")

    lines += [
        "\n## 3. Breakpoint count distribution\n",
        "| quantile | value |\n|---|---:|\n",
    ]
    for q in [0.00, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00]:
        lines.append(f"| {q*100:.0f}% | {int(bp.quantile(q))} |\n")

    # Add taxonomy if a genome metadata file is provided
    if args.genomes and Path(args.genomes).exists():
        meta = pd.read_csv(args.genomes, sep="\t")
        meta = meta.rename(columns={"acc": "q_acc", "species": "q_species"})
        syn = syn.merge(meta[["q_acc", "q_species"]], on="q_acc", how="left")
        meta = meta.rename(columns={"q_acc": "r_acc", "q_species": "r_species"})
        syn = syn.merge(meta[["r_acc", "r_species"]], on="r_acc", how="left")
        syn["same_species"] = syn["q_species"] == syn["r_species"]

    # Overlap with the original near-closed seed pairs
    seed_ids = set(pairs["pairid"])
    syn_ids = set(syn["pairid"])
    overlap = seed_ids & syn_ids
    lines += [
        "\n## 4. Overlap with near-closed seed pairs\n",
        f"- Seed pairs: {len(seed_ids):,}\n",
        f"- Seed pairs present in all-vs-all output: {len(overlap):,} ({len(overlap)/len(seed_ids)*100:.1f}%)\n",
    ]

    if "same_species" in syn.columns:
        same = syn[syn["same_species"] == True]
        diff = syn[syn["same_species"] != True]
        lines += [
            "\n## 5. Same-species vs different-species pairs\n",
            f"- Same-species pairs: {len(same):,} ({len(same)/len(syn)*100:.1f}%)\n",
            f"- Different-species pairs: {len(diff):,} ({len(diff)/len(syn)*100:.1f}%)\n",
            f"- Same-species median inverted_fraction: {same['syn2b_raw_inverted_fraction'].median():.4f}\n" if len(same) else "",
            f"- Different-species median inverted_fraction: {diff['syn2b_raw_inverted_fraction'].median():.4f}\n" if len(diff) else "",
        ]

    # Top rearranged pairs
    top = syn.sort_values(["syn2b_raw_inverted_fraction", "syn2b_breakpoints"], ascending=False).head(20)
    lines += [
        "\n## 6. Top 20 pairs by raw_inverted_fraction\n",
        "| pairid | inverted_fraction | breakpoints | scj_distance | shared_tags |\n",
        "|---|---:|---:|---:|---:|\n",
    ]
    for _, r in top.iterrows():
        lines.append(
            f"| {r['pairid']} | {r['syn2b_raw_inverted_fraction']:.4f} | "
            f"{int(r['syn2b_breakpoints'])} | {r['syn2b_scj_distance']:.2f} | {int(r['syn2b_shared_tags'])} |\n"
        )

    # Species identity for top pairs
    lines += [
        "\n## 7. Species identity for top pairs\n",
        "| pairid | same_species | q_species | r_species |\n",
        "|---|---|---|---|\n",
    ]
    for _, r in top.iterrows():
        same = r.get("same_species", "NA")
        lines.append(
            f"| {r['pairid']} | {same} | {r.get('q_species','NA')} | {r.get('r_species','NA')} |\n"
        )

    Path(args.out_md).write_text("".join(lines))
    print(f"wrote {args.out_md}")

    # Write full ranking TSV
    rank = syn[[
        "pairid",
        "syn2b_raw_inverted_fraction",
        "syn2b_inverted_fraction",
        "syn2b_breakpoints",
        "syn2b_breakpoint_density",
        "syn2b_scj_distance",
        "syn2b_scj_corrected",
        "syn2b_shared_tags",
        "syn2b_observable_fraction",
    ]].copy()
    rank = rank.merge(pairs[["pairid", "q_species", "r_species", "same_species"]], on="pairid", how="left")
    rank = rank.sort_values("syn2b_raw_inverted_fraction", ascending=False)
    rank.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
