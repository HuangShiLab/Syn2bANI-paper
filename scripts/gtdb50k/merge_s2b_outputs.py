#!/usr/bin/env python3
"""Merge per-pair syn2bani --verbose TSV outputs into one table.

Input:  results/gtdb50k/s2b_out/*.tsv
Output: results/gtdb50k/s2b_50k.tsv
"""
import argparse
import csv
import glob
import os
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", default="results/gtdb50k/s2b_out")
    p.add_argument("--out", default="results/gtdb50k/s2b_50k.tsv")
    p.add_argument("--expected", type=int, default=None,
                   help="expected number of rows (for a sanity check)")
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.indir, "*.tsv")))
    if not files:
        print(f"no .tsv files found in {args.indir}", file=sys.stderr)
        sys.exit(1)

    written_header = False
    n_rows = 0
    with open(args.out, "w", newline="") as outfh:
        writer = None
        for path in files:
            try:
                with open(path) as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    if not written_header:
                        writer = csv.DictWriter(outfh, fieldnames=reader.fieldnames, delimiter="\t")
                        writer.writeheader()
                        written_header = True
                    for row in reader:
                        writer.writerow(row)
                        n_rows += 1
            except Exception as e:
                print(f"skip {path}: {e}", file=sys.stderr)

    print(f"merged {n_rows} rows from {len(files)} files -> {args.out}")
    if args.expected is not None:
        if n_rows != args.expected:
            print(f"WARNING: expected {args.expected} rows, got {n_rows}", file=sys.stderr)
        else:
            print("row count matches expected")


if __name__ == "__main__":
    main()
