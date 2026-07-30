#!/usr/bin/env python3
"""Evaluate Syn2bANI output on the oral/gut validation set.

Joins the Syn2bANI ani (with optional ani_cal) output to the validation TSV
by mapping first sequence IDs back to assembly accessions, then reports
MAE/RMSE/bias/R2 against FastANI and skani references.

Usage:
    python3 scripts/evaluate_oral_gut.py \
        --s2b-ani results/oral_gut_11enzyme_cal.tsv \
        --validation data/oral_gut_validation_merged_v8.tsv \
        --genome-dir /lustre1/g/aos_shihuang/data/validation_oral_gut/genomes \
        --outdir results/oral_gut_eval
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def build_seqid_to_asm(genome_dir, suffix=".fna"):
    """Map first sequence ID of each genome FASTA to assembly accession."""
    genome_dir = Path(genome_dir)
    seqid_to_asm = {}
    for fna in sorted(genome_dir.glob(f"*{suffix}")):
        asm = fna.name[: -len(suffix)]
        with open(fna) as fh:
            first = fh.readline()
        if not first.startswith(">"):
            continue
        seqid = first[1:].split()[0]
        seqid_to_asm[seqid] = asm
    return seqid_to_asm


def read_s2b(path, seqid_to_asm):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"query": "query_seqid", "reference": "ref_seqid"})
    df["query_asm"] = df["query_seqid"].map(seqid_to_asm)
    df["ref_asm"] = df["ref_seqid"].map(seqid_to_asm)
    df = df.dropna(subset=["query_asm", "ref_asm"])
    cols = ["query_asm", "ref_asm", "ani"]
    if "ani_cal" in df.columns:
        cols.append("ani_cal")
    for c in ["synteny_blocks", "synteny_score", "breakpoint_count",
              "max_block_anchors", "mean_block_anchors"]:
        if c in df.columns:
            cols.append(c)
    return df[cols]


def metrics(pred, truth):
    pred = np.asarray(pred)
    truth = np.asarray(truth)
    mask = np.isfinite(pred) & np.isfinite(truth)
    pred, truth = pred[mask], truth[mask]
    if len(pred) == 0:
        return {"n": 0, "mae": None, "rmse": None, "bias": None, "r2": None}
    diff = pred - truth
    return {
        "n": int(len(diff)),
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "bias": float(np.mean(diff)),
        "r2": float(1 - np.sum(diff ** 2) / np.sum((truth - np.mean(truth)) ** 2))
        if np.std(truth) > 0
        else None,
    }


def evaluate(df, truth_col, pred_col):
    rows = []
    for band, sub in df.groupby("band"):
        if len(sub) < 2:
            continue
        rows.append({"band": band, **metrics(sub[pred_col].values, sub[truth_col].values)})
    rows.append({"band": "all", **metrics(df[pred_col].values, df[truth_col].values)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s2b-ani", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--outdir", default="results/oral_gut_eval")
    args = ap.parse_args()

    print("building seqid -> assembly map...")
    seqid_to_asm = build_seqid_to_asm(args.genome_dir)
    print(f"mapped {len(seqid_to_asm)} genomes")

    print("reading Syn2bANI output...")
    s2b = read_s2b(args.s2b_ani, seqid_to_asm)
    print(f"{len(s2b)} pairs with mapping")

    print("reading validation TSV...")
    val = pd.read_csv(args.validation, sep="\t")
    val = val.rename(columns={"query": "query_asm", "reference": "ref_asm"})
    # Use the ANI band label stored in the validation sheet, or fall back to binning
    if "label" in val.columns:
        val["band"] = val["label"]
    else:
        val["band"] = pd.cut(val["skani_ani"], bins=[0, 0.85, 0.90, 0.95, 1.0],
                             labels=["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"])

    merged = val[["query_asm", "ref_asm", "band", "fastani_ani", "skani_ani",
                  "q_species", "r_species", "q_genus", "r_genus"]].merge(
        s2b, on=["query_asm", "ref_asm"], how="inner"
    )
    print(f"merged {len(merged)} pairs")
    if len(merged) == 0:
        print("no merged pairs; aborting", file=sys.stderr)
        return 1

    # Reference tools report fractions; Syn2bANI reports percentages
    merged["fastani_ani_pct"] = merged["fastani_ani"] * 100.0
    merged["skani_ani_pct"] = merged["skani_ani"] * 100.0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(outdir / "eval_pairs.tsv", sep="\t", index=False)

    summary = {"n_pairs": len(merged)}
    methods = [("ani", "Syn2bANI_default")]
    if "ani_cal" in merged.columns:
        methods.append(("ani_cal", "Syn2bANI_cal"))

    # vs FastANI
    fastani_rows = []
    for col, name in methods:
        ev = evaluate(merged, "fastani_ani_pct", col)
        for _, row in ev.iterrows():
            fastani_rows.append({"method": name, **row.to_dict()})
    summary["vs_fastani"] = fastani_rows

    # vs skani
    skani_rows = []
    for col, name in methods:
        ev = evaluate(merged, "skani_ani_pct", col)
        for _, row in ev.iterrows():
            skani_rows.append({"method": name, **row.to_dict()})
    summary["vs_skani"] = skani_rows

    with open(outdir / "evaluation_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {outdir / 'evaluation_summary.json'}")

    print("\n=== MAE (%) vs FastANI by band ===")
    print(pd.DataFrame(fastani_rows).pivot(index="band", columns="method", values="mae").to_string())
    print("\n=== MAE (%) vs skani by band ===")
    print(pd.DataFrame(skani_rows).pivot(index="band", columns="method", values="mae").to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
