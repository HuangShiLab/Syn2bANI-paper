#!/usr/bin/env python3
"""Derive structural-variation truth from dnadiff 1-to-1 coords files.

For each pair in the GTDB 50k held-out benchmark, parse dd.1coords and
count:
  - dnadiff_blocks: number of 1-to-1 aligned blocks
  - dnadiff_breakpoints: number of collinearity breaks between consecutive
    query-ordered blocks (strand inversion or non-collinear reference jump)
  - dnadiff_large_indels: gaps between consecutive query-ordered blocks
    >= 1000 bp on either genome
  - dnadiff_anchor_adjacency: 1 - breakpoints/(blocks-1) when blocks>1, else 1.0

Compares these to syn2bani's breakpoint_count and anchor_adjacency.

Inputs:
  results/gtdb50k/out/{pairid}/dd.1coords
  results/gtdb50k/s2b_50k.tsv

Outputs:
  results/gtdb50k/sv_truth_50k.tsv
  results/gtdb50k/SV_EVALUATION_REPORT.md
"""
import os
import re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")
OUTDIR = os.path.join(RES, "out")
MIN_INDEL = 1000


def parse_one(pairid):
    path = os.path.join(OUTDIR, pairid, "dd.1coords")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("="):
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            try:
                s1, e1, s2, e2, l1, l2, ident, lenq, lenr, covq, covr = map(float, parts[:11])
            except ValueError:
                continue
            tagq, tagr = parts[11], parts[12]
            rows.append((int(s1), int(e1), int(s2), int(e2), tagq, tagr))
    if not rows:
        return None

    # sort by query coordinate
    rows.sort(key=lambda x: x[0])
    n = len(rows)
    if n == 1:
        return dict(blocks=1, breakpoints=0, large_indels=0, anchor_adjacency=1.0)

    breakpoints = 0
    large_indels = 0
    for i in range(n - 1):
        qs1, qe1, rs1, re1, _, _ = rows[i]
        qs2, qe2, rs2, re2, _, _ = rows[i + 1]
        qgap = qs2 - qe1 - 1
        rgap = 0
        strand1 = 1 if re1 >= rs1 else -1
        strand2 = 1 if re2 >= rs2 else -1
        if strand1 != strand2:
            breakpoints += 1
        else:
            if strand1 == 1:
                rgap = rs2 - re1 - 1
                if rgap < -1000:  # non-collinear: reference goes backward
                    breakpoints += 1
            else:
                rgap = rs1 - re2 - 1
                if rgap < -1000:
                    breakpoints += 1
        if qgap >= MIN_INDEL or abs(rgap) >= MIN_INDEL:
            large_indels += 1

    score = max(0.0, 1.0 - breakpoints / (n - 1))
    return dict(blocks=n, breakpoints=breakpoints, large_indels=large_indels,
                anchor_adjacency=score)


def main():
    s2b = pd.read_csv(os.path.join(RES, "s2b_50k.tsv"), sep="\t")
    s2b["pairid"] = s2b["pairid"] if "pairid" in s2b.columns else None
    # s2b_50k has pairid column from merge
    s2b = s2b[["pairid", "breakpoint_count", "anchor_adjacency"]].copy()

    pair_dirs = [d for d in os.listdir(OUTDIR) if os.path.isdir(os.path.join(OUTDIR, d))]
    print(f"found {len(pair_dirs)} pair directories")

    records = []
    for i, pid in enumerate(sorted(pair_dirs)):
        if i % 5000 == 0 and i > 0:
            print(f"  processed {i}")
        rec = parse_one(pid)
        if rec is None:
            continue
        rec["pairid"] = pid
        records.append(rec)

    sv = pd.DataFrame(records)
    sv = sv.rename(columns={"blocks": "dnadiff_blocks", "breakpoints": "dnadiff_breakpoints",
                            "large_indels": "dnadiff_large_indels", "anchor_adjacency": "dnadiff_anchor_adjacency"})
    sv.to_csv(os.path.join(RES, "sv_truth_50k.tsv"), sep="\t", index=False)

    merged = s2b.merge(sv, on="pairid", how="inner")
    print(f"merged with s2b: {len(merged)} pairs")

    # metrics
    bp_mae = np.mean(np.abs(merged["breakpoint_count"] - merged["dnadiff_breakpoints"]))
    bp_r = float(np.corrcoef(merged["breakpoint_count"], merged["dnadiff_breakpoints"])[0, 1])
    syn_r = float(np.corrcoef(merged["anchor_adjacency"], merged["dnadiff_anchor_adjacency"])[0, 1])
    # classification: any rearrangement
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

    # by band
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
         f"- anchor_adjacency Pearson r: {syn_r:.3f}",
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
