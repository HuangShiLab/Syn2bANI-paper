#!/usr/bin/env python3
"""Evaluate Syn2bANI default output and skani against ANIm truth.

Uses the stratified sample file to get genome paths and bands, maps assembly
accessions to first sequence IDs (so it can join the syn2bani ani output),
and writes a JSON summary plus a per-pair TSV.

Usage:
    python3 scripts/evaluate_vs_anim.py \
        --sample results/sample_anim_truth.tsv \
        --anim-dir anim_results \
        --s2b-ani results/gtdb_r207_100k_11enzyme.tsv \
        --matrix results/matrix_gtdb_r207_100k_v8_final.tsv \
        --outdir results/panel_by_band
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np


def build_seqid_map(paths):
    """Map assembly accession to first sequence ID in FASTA."""
    seqid_map = {}
    for asm, fna in sorted(paths.items()):
        if not Path(fna).exists():
            continue
        with open(fna) as fh:
            first = fh.readline()
        if not first.startswith(">"):
            continue
        seqid_map[asm] = first[1:].split()[0]
    return seqid_map


def read_truth(anim_dir):
    truth = {}
    for p in Path(anim_dir).glob("anim_*.tsv"):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                q, r, ident = f[0], f[1], f[2]
                try:
                    v = float(ident)
                except ValueError:
                    continue
                if not v > 0:
                    continue
                truth[(q, r)] = v
    return truth


def read_sample(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "query": "query_asm",
        "reference": "ref_asm",
        "band": "band",
        "group": "group",
    })
    paths = {}
    for _, row in df.iterrows():
        paths[row["query_asm"]] = row["query_path"]
        paths[row["ref_asm"]] = row["ref_path"]
    return df, paths


def read_s2b_ani(path, seqid_to_asm):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "query": "query_seqid",
        "reference": "ref_seqid",
        "ani": "s2b_ani",
    })
    df["query_asm"] = df["query_seqid"].map(seqid_to_asm)
    df["ref_asm"] = df["ref_seqid"].map(seqid_to_asm)
    df = df.dropna(subset=["query_asm", "ref_asm"])
    return df[["query_asm", "ref_asm", "s2b_ani"]]


def metrics(pred, truth):
    pred = np.asarray(pred)
    truth = np.asarray(truth)
    diff = pred - truth
    return {
        "n": int(len(diff)),
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "bias": float(np.mean(diff)),
        "r2": float(1 - np.sum(diff ** 2) / np.sum((truth - np.mean(truth)) ** 2)) if np.std(truth) > 0 else None,
    }


def evaluate(df):
    rows = []
    for band, sub in df.groupby("band"):
        if len(sub) < 2:
            continue
        rows.append({
            "band": band,
            **metrics(sub["s2b_ani"].values, sub["anim_ani"].values),
            "method": "Syn2bANI_default",
        })
        rows.append({
            "band": band,
            **metrics(sub["skani_ani"].values, sub["anim_ani"].values),
            "method": "skani",
        })
    rows.append({
        "band": "all",
        **metrics(df["s2b_ani"].values, df["anim_ani"].values),
        "method": "Syn2bANI_default",
    })
    rows.append({
        "band": "all",
        **metrics(df["skani_ani"].values, df["anim_ani"].values),
        "method": "skani",
    })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--anim-dir", default="anim_results")
    ap.add_argument("--s2b-ani", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--outdir", default="results/panel_by_band")
    args = ap.parse_args()

    sample_df, paths = read_sample(args.sample)
    print(f"sample: {len(sample_df)} pairs, {len(paths)} genomes")

    print("building assembly -> seqid map...")
    seqid_map = build_seqid_map(paths)
    print(f"mapped {len(seqid_map)} assemblies")
    seqid_to_asm = {v: k for k, v in seqid_map.items()}

    print("reading ANIm truth...")
    truth = read_truth(args.anim_dir)
    print(f"truth pairs: {len(truth)}")

    print("reading Syn2bANI ani output...")
    s2b_df = read_s2b_ani(args.s2b_ani, seqid_to_asm)
    print(f"s2b pairs with mapping: {len(s2b_df)}")

    print("reading matrix...")
    mat = pd.read_csv(args.matrix, sep="\t")
    mat = mat.rename(columns={
        "query": "query_asm",
        "reference": "ref_asm",
        "skani_ani": "skani_ani",
        "q_phylum": "q_phylum",
        "r_phylum": "r_phylum",
    })

    # Build evaluation frame
    truth_rows = []
    for (q, r), v in truth.items():
        truth_rows.append({"query_asm": q, "ref_asm": r, "anim_ani": v})
    truth_df = pd.DataFrame(truth_rows)

    merged = sample_df[["query_asm", "ref_asm", "band", "group"]].merge(
        truth_df, on=["query_asm", "ref_asm"], how="inner"
    ).merge(
        s2b_df, on=["query_asm", "ref_asm"], how="inner"
    ).merge(
        mat[["query_asm", "ref_asm", "skani_ani", "q_phylum", "r_phylum"]],
        on=["query_asm", "ref_asm"], how="inner",
    )
    # skani reports ANI as a fraction; Syn2bANI and ANIm report percentages
    merged["skani_ani"] = merged["skani_ani"] * 100.0
    print(f"merged pairs: {len(merged)}")
    if len(merged) == 0:
        print("no merged pairs; aborting", file=sys.stderr)
        return 1

    merged["phylum_pair"] = merged["q_phylum"] + " vs " + merged["r_phylum"]
    merged.to_csv(Path(args.outdir) / "eval_pairs.tsv", sep="\t", index=False)

    summary = {
        "n_pairs": len(merged),
        "by_band": evaluate(merged).to_dict(orient="records"),
    }

    # Per-phylum summary for the most common phyla
    phylum_counts = merged["q_phylum"].value_counts()
    phylum_rows = []
    for phylum in phylum_counts.head(15).index:
        sub = merged[merged["q_phylum"] == phylum]
        phylum_rows.append({
            "phylum": phylum,
            **metrics(sub["s2b_ani"].values, sub["anim_ani"].values),
            "method": "Syn2bANI_default",
        })
        phylum_rows.append({
            "phylum": phylum,
            **metrics(sub["skani_ani"].values, sub["anim_ani"].values),
            "method": "skani",
        })
    summary["by_phylum"] = phylum_rows

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "evaluation_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {outdir / 'evaluation_summary.json'}")

    # Print readable table
    print("\n=== MAE (%) vs ANIm by band ===")
    band_df = evaluate(merged)
    print(band_df.pivot(index="band", columns="method", values="mae").to_string())

    return 0


if __name__ == "__main__":
    sys.exit(main())
