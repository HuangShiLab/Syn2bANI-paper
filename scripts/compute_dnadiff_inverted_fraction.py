#!/usr/bin/env python3
"""Compute inverted aligned fraction from dnadiff dd.1coords files (parallel).

For each pair, parse dd.1coords and compute:
    inverted_bp  = sum |E2 - S2| over blocks where query is reversed (S2 > E2)
    aligned_bp   = sum |E2 - S2| over all blocks
    inverted_fraction = inverted_bp / aligned_bp

This is a length-weighted ratio, invariant to alignment fragmentation.

Inputs:
    results/gtdb50k/out/<pairid>/dd.1coords
    results/gtdb50k/pairs_50k.tsv

Outputs:
    results/gtdb50k/dnadiff_inverted_fraction.tsv
"""

import argparse
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "results" / "gtdb50k"


def parse_one(args):
    pairid, outdir = args
    path = Path(outdir) / pairid / "dd.1coords"
    if not path.exists():
        return None

    inverted_bp = 0
    aligned_bp = 0
    n_blocks = 0
    n_inverted = 0

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("="):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            try:
                s1, e1, s2, e2, l1, l2, ident, lenq, lenr, covq, covr = map(float, parts[:11])
            except ValueError:
                continue

            bp = abs(int(e2) - int(s2)) + 1
            aligned_bp += bp
            n_blocks += 1

            if s2 > e2:
                inverted_bp += bp
                n_inverted += 1

    if aligned_bp == 0:
        return None

    return {
        "pairid": pairid,
        "dnadiff_blocks": n_blocks,
        "dnadiff_inverted_blocks": n_inverted,
        "dnadiff_aligned_bp": aligned_bp,
        "dnadiff_inverted_bp": inverted_bp,
        "dnadiff_inverted_fraction": inverted_bp / aligned_bp,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", default=str(RES / "pairs_50k.tsv"))
    p.add_argument("--outdir", default=str(RES / "out"))
    p.add_argument("--outfile", default=str(RES / "dnadiff_inverted_fraction.tsv"))
    p.add_argument("--workers", type=int, default=min(32, cpu_count() or 1))
    args = p.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    if "pairid" not in pairs.columns:
        pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    pairids = pairs["pairid"].tolist()

    tasks = [(pid, args.outdir) for pid in pairids]

    print(f"parsing {len(tasks)} pairs from {args.outdir} with {args.workers} workers ...", flush=True)

    records = []
    with Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(parse_one, tasks, chunksize=100)):
            if i % 5000 == 0 and i > 0:
                print(f"  processed {i}/{len(tasks)}", flush=True)
            if rec is not None:
                records.append(rec)

    df = pd.DataFrame(records)
    keep_cols = [c for c in ["pairid", "band", "phylum"] if c in pairs.columns]
    df = pairs[keep_cols].merge(df, on="pairid", how="left")
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"wrote {args.outfile}: {len(df)} pairs, {df['dnadiff_inverted_fraction'].notna().sum()} parseable")

    if "band" in df.columns:
        print("\nSummary by band:")
        print(df.groupby("band").agg(
            n=("dnadiff_inverted_fraction", "size"),
            parseable=("dnadiff_inverted_fraction", lambda x: x.notna().sum()),
            mean_inverted_fraction=("dnadiff_inverted_fraction", "mean"),
            median_inverted_fraction=("dnadiff_inverted_fraction", "median"),
        ).round(5))


if __name__ == "__main__":
    main()
