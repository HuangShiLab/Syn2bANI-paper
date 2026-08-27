#!/usr/bin/env python3
"""Merge Syn2bANI, dnadiff, and minimap2 SV metrics and compute correlations.

Inputs:
  results/gtdb50k/s2b_50k.tsv
  results/gtdb50k/sv_truth_50k.tsv
  results/gtdb50k/sv_truth_50k_min{gap}.tsv  (optional, one or more)
  results/gtdb50k/minimap2_rows/slice_*.tsv

Outputs:
  results/gtdb50k/SV_COMPARISON_REPORT.md
  results/gtdb50k/sv_comparison_merged.tsv
"""
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SYN2BANI_ROOT", HERE.parent.parent))
RES = ROOT / "results" / "gtdb50k"


def load_minimap2():
    files = sorted(glob.glob(str(RES / "minimap2_rows" / "slice_*.tsv")))
    if not files:
        return None
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t", header=None,
                         names=["pairid", "mm2_blocks", "mm2_breakpoints",
                                "mm2_large_indels", "mm2_synteny_score", "mm2_status"])
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def safe_corr(x, y, method="pearson"):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if method == "pearson":
        return pearsonr(x[mask], y[mask])[0]
    return spearmanr(x[mask], y[mask])[0]


def main():
    s2b = pd.read_csv(RES / "s2b_50k.tsv", sep="\t")
    s2b = s2b[["pairid", "ani_gated", "af_query", "synteny_score", "synteny_blocks", "breakpoint_count"]].copy()

    dnadiff = pd.read_csv(RES / "sv_truth_50k.tsv", sep="\t")
    merged = s2b.merge(dnadiff, on="pairid", how="left")

    # optional filtered dnadiff files
    for gap_file in sorted(RES.glob("sv_truth_50k_min*.tsv")):
        gap_df = pd.read_csv(gap_file, sep="\t")
        merged = merged.merge(gap_df, on="pairid", how="left")

    mm2 = load_minimap2()
    if mm2 is not None:
        merged = merged.merge(mm2, on="pairid", how="left")

    merged.to_csv(RES / "sv_comparison_merged.tsv", sep="\t", index=False)

    L = ["# SV method comparison on GTDB 50k held-out pairs",
         f"Pairs with Syn2bANI output: {len(merged):,}",
         ""]
    cols = {
        "dnadiff_breakpoints": "dnadiff (all gaps)",
        "dnadiff_synteny_score": "dnadiff synteny (all gaps)",
        "mm2_breakpoints": "minimap2",
        "mm2_synteny_score": "minimap2 synteny",
    }
    # add filtered dnadiff columns if present
    for c in merged.columns:
        if c.startswith("dnadiff_breakpoints_min"):
            cols[c] = f"dnadiff min-gap {c.replace('dnadiff_breakpoints_min', '')} bp"
        if c.startswith("dnadiff_large_indels_min"):
            cols[c] = f"dnadiff large indels min-gap {c.replace('dnadiff_large_indels_min', '')} bp"

    L.append("## Breakpoint count correlations with Syn2bANI breakpoint_count")
    for col, label in cols.items():
        if col not in merged.columns or merged[col].isna().all():
            continue
        r = safe_corr(merged["breakpoint_count"].to_numpy(float),
                      merged[col].to_numpy(float), "pearson")
        rs = safe_corr(merged["breakpoint_count"].to_numpy(float),
                       merged[col].to_numpy(float), "spearman")
        mae = np.nanmean(np.abs(merged["breakpoint_count"] - merged[col]))
        L.append(f"- {label}: Pearson r={r:.4f}, Spearman r={rs:.4f}, MAE={mae:.1f}")

    L.append("")
    L.append("## Correlations with alignment-based synteny/coverage scores")
    syn_metrics = {
        "synteny_score": "Syn2bANI synteny_score (anchor-adjacency conservation)",
        "af_query": "Syn2bANI af_query (base-pair chain coverage)",
        "synteny_blocks": "Syn2bANI synteny_blocks",
    }
    syn_cols = [c for c, _ in cols.items() if "synteny" in c]
    for s2b_col, s2b_label in syn_metrics.items():
        L.append(f"\n### {s2b_label}")
        for col in syn_cols:
            if col not in merged.columns or merged[col].isna().all():
                continue
            label = cols[col]
            r = safe_corr(merged[s2b_col].to_numpy(float),
                          merged[col].to_numpy(float), "pearson")
            rs = safe_corr(merged[s2b_col].to_numpy(float),
                           merged[col].to_numpy(float), "spearman")
            L.append(f"- {label}: Pearson r={r:.4f}, Spearman r={rs:.4f}")

    L.append("")
    L.append("## Summary statistics")
    for col in ["breakpoint_count", "dnadiff_breakpoints", "mm2_breakpoints"]:
        if col in merged.columns and not merged[col].isna().all():
            L.append(f"- {col}: mean={merged[col].mean():.1f}, median={merged[col].median():.1f}")

    with open(RES / "SV_COMPARISON_REPORT.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
