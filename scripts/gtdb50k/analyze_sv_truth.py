#!/usr/bin/env python3
"""Analyze SV truth vs syn2bani predictions for the 43,334 held-out pairs."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")

def main():
    s2b = pd.read_csv(os.path.join(RES, "s2b_50k.tsv"), sep="\t")[["pairid", "breakpoint_count", "synteny_score"]].copy()
    sv = pd.read_csv(os.path.join(RES, "sv_truth_50k.tsv"), sep="\t")
    merged = s2b.merge(sv, on="pairid", how="inner")
    print(f"merged pairs: {len(merged):,}")

    bp_mae = np.mean(np.abs(merged["breakpoint_count"] - merged["dnadiff_breakpoints"]))
    bp_r = float(np.corrcoef(merged["breakpoint_count"], merged["dnadiff_breakpoints"])[0, 1])
    syn_r = float(np.corrcoef(merged["synteny_score"], merged["dnadiff_synteny_score"])[0, 1])

    truth_has = (merged["dnadiff_breakpoints"] > 0).astype(int)
    pred_has = (merged["breakpoint_count"] > 0).astype(int)
    tp = ((truth_has == 1) & (pred_has == 1)).sum()
    fp = ((truth_has == 0) & (pred_has == 1)).sum()
    fn = ((truth_has == 1) & (pred_has == 0)).sum()
    tn = ((truth_has == 0) & (pred_has == 0)).sum()
    precision = tp / (tp + fp) if (tp + fp) else np.nan
    recall = tp / (tp + fn) if (tp + fn) else np.nan
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan

    pairs = pd.read_csv(os.path.join(RES, "pairs_50k.tsv"), sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    merged = merged.merge(pairs[["pairid", "band"]], on="pairid", how="left")

    L = ["# SV truth evaluation against dnadiff 1-to-1 coords (GTDB 50k held-out)",
         f"Pairs with parseable dd.1coords: {len(sv):,}",
         f"Pairs with syn2bani output: {len(merged):,}",
         "",
         "## Overall metrics",
         f"- breakpoint_count MAE vs dnadiff: {bp_mae:.3f}",
         f"- breakpoint_count Pearson r: {bp_r:.3f}",
         f"- synteny_score Pearson r: {syn_r:.3f}",
         f"- Rearrangement detection (truth > 0 vs pred > 0): precision={precision:.3f}, recall={recall:.3f}, F1={f1:.3f}, specificity={specificity:.3f}",
         f"  TP={tp}, FP={fp}, FN={fn}, TN={tn}",
         "",
         "## Per-band breakpoint_count MAE"]
    for band in ["80-85", "85-90", "90-95", "95-100"]:
        sub = merged[merged["band"] == band]
        if len(sub):
            mae = np.mean(np.abs(sub["breakpoint_count"] - sub["dnadiff_breakpoints"]))
            L.append(f"- {band}: {mae:.3f} (n={len(sub)})")

    with open(os.path.join(RES, "SV_EVALUATION_REPORT.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()
