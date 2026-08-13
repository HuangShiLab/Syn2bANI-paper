#!/usr/bin/env python3
"""Evaluate Syn2bANI v8 / skani / FastANI against ANIm (dnadiff) truth on the
15 mid-ANI (85-95%) oral/gut validation pairs.

Inputs (copied back from HPC /lustre1/g/aos_shihuang/data/validation_mid_ani):
  --anim      anim/anim_truth.tsv        (query, reference, anim_ani, ...)
  --s2b       anim/syn2bani_v8.tsv       (syn2bani ani --verbose output)
  --skani     anim/skani.tsv             (skani dist on same genome set)
  --fastani   mid_ani_matrix_tools.tsv   (existing FastANI reference values)
  --pairs     mid_ani_pairs_85_95.tsv    (the 15 evaluation pairs)

Writes <outdir>/anim_midani_evaluation.tsv and prints a metrics table.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def metrics(pred, truth):
    pred, truth = np.asarray(pred, float), np.asarray(truth, float)
    d = pred - truth
    return {
        "n": len(d),
        "MAE": round(float(np.mean(np.abs(d))), 4),
        "RMSE": round(float(np.sqrt(np.mean(d ** 2))), 4),
        "bias": round(float(np.mean(d)), 4),
        "r": round(float(np.corrcoef(pred, truth)[0, 1]), 4) if len(d) > 2 else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anim", required=True)
    ap.add_argument("--s2b", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--fastani", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--map", required=True, help="accession<->first-seqid TSV")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    pairs["key"] = list(zip(pairs["query"], pairs["reference"]))

    anim = pd.read_csv(args.anim, sep="\t")
    anim["key"] = list(zip(anim["query"], anim["reference"]))

    df = pairs[["key", "query", "reference", "q_species", "r_species"]].merge(
        anim[["key", "anim_ani", "anim_aligned_query", "anim_aligned_ref"]], on="key", how="inner")

    # Syn2bANI v8: identify genomes by first-contig seqid -> accession map
    seq2acc = pd.read_csv(args.map, sep="\t", names=["accession", "seqid"])
    seq2acc = dict(zip(seq2acc["seqid"], seq2acc["accession"]))
    s2b = pd.read_csv(args.s2b, sep="\t")
    s2b["query"] = s2b["query"].map(seq2acc)
    s2b["reference"] = s2b["reference"].map(seq2acc)
    s2b["key"] = list(zip(s2b["query"], s2b["reference"]))
    keep = ["key", "ani", "ani_uniform", "af_query", "af_reference", "std_err"]
    for c in ("retention", "ani_from_loss", "ani_from_hist", "n_anchors", "n_chains", "n_tags", "flag"):
        if c in s2b.columns:
            keep.append(c)
    s2b = s2b[keep].rename(columns={"ani": "s2b_ani"})
    df = df.merge(s2b, on="key", how="left")

    # skani (fresh run on the same genome set); use file-path stems = accessions
    ska = pd.read_csv(args.skani, sep="\t")
    ska["query"] = ska["Query_file"].astype(str).str.split("/").str[-1].str.replace(r"\.fna$", "", regex=True)
    ska["reference"] = ska["Ref_file"].astype(str).str.split("/").str[-1].str.replace(r"\.fna$", "", regex=True)
    ska["key"] = list(zip(ska["query"], ska["reference"]))
    df = df.merge(ska[["key", "ANI"]].rename(columns={"ANI": "skani_ani"}), on="key", how="left")

    # FastANI reference (existing matrix; fraction -> percent)
    fa = pd.read_csv(args.fastani, sep="\t")
    fa_cols = {c.lower(): c for c in fa.columns}
    qcol, rcol = fa_cols.get("query"), fa_cols.get("reference")
    facol = next((c for c in fa.columns if "fastani" in c.lower()), None)
    if facol:
        fa["key"] = list(zip(fa[qcol], fa[rcol]))
        df = df.merge(fa[["key", facol]].rename(columns={facol: "fastani_ani"}), on="key", how="left")
        if df["fastani_ani"].max() <= 1.5:
            df["fastani_ani"] *= 100.0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "anim_midani_evaluation.tsv"
    df.drop(columns=["key"]).to_csv(out_tsv, sep="\t", index=False)

    methods = [c for c in ("s2b_ani", "ani_uniform", "skani_ani", "fastani_ani") if c in df.columns]
    rows = []
    for m in methods:
        sub = df.dropna(subset=[m, "anim_ani"])
        if len(sub):
            rows.append({"method": m, **metrics(sub[m], sub["anim_ani"])})
    res = pd.DataFrame(rows)
    res.to_csv(outdir / "anim_midani_metrics.tsv", sep="\t", index=False)
    print(res.to_string(index=False))
    print(f"\nwrote {out_tsv} ({len(df)} pairs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
