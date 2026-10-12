#!/usr/bin/env python3
"""Diagnose block-pair coverage for one cluster in a normalized table.

This is intentionally based on the final canonical pair table rather than
worker checkpoints, because old and current task layouts can coexist in the
raw output directory. The final table is already globally deduplicated.
"""

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--block-size", type=int, default=500)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    safe = args.cluster.replace("/", "_")
    accessions = Path(args.workdir, "effective_cluster_accessions", safe + ".txt").read_text().split()
    index = {acc: i for i, acc in enumerate(accessions)}
    n = len(accessions)
    nblocks = (n + args.block_size - 1) // args.block_size
    expected = {}
    for bi in range(nblocks):
        for bj in range(bi, nblocks):
            ni = min(args.block_size, n - bi * args.block_size)
            nj = min(args.block_size, n - bj * args.block_size)
            expected[(bi, bj)] = ni * (ni - 1) // 2 if bi == bj else ni * nj

    observed = defaultdict(int)
    total = in_cluster = invalid_id = 0
    opener = gzip.open if str(args.table).endswith(".gz") else open
    with opener(args.table, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ia, ib = header.index("genome_A"), header.index("genome_B")
        for line in fh:
            total += 1
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(ia, ib):
                continue
            a, b = row[ia], row[ib]
            x, y = index.get(a), index.get(b)
            if x is None or y is None:
                if (a in index) != (b in index):
                    invalid_id += 1
                continue
            if x == y:
                invalid_id += 1
                continue
            in_cluster += 1
            lo, hi = sorted((x, y))
            observed[(lo // args.block_size, hi // args.block_size)] += 1

    fields = [
        "block_i", "block_j", "expected_pairs", "observed_pairs",
        "missing_pairs", "coverage", "status",
    ]
    tmp = Path(str(args.output) + ".tmp")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(expected):
            exp, obs = expected[key], observed[key]
            missing = max(0, exp - obs)
            status = "complete" if missing == 0 and obs == exp else (
                "partial" if obs else "missing")
            if obs > exp:
                status = "unexpected_extra_rows"
            writer.writerow({
                "block_i": key[0], "block_j": key[1],
                "expected_pairs": exp, "observed_pairs": obs,
                "missing_pairs": missing,
                "coverage": obs / exp if exp else "",
                "status": status,
            })
    tmp.replace(args.output)

    complete = sum(observed.get(k, 0) == v for k, v in expected.items())
    partial = sum(0 < observed.get(k, 0) < v for k, v in expected.items())
    missing = sum(observed.get(k, 0) == 0 for k, v in expected.items())
    extra = sum(observed.get(k, 0) > v for k, v in expected.items())
    report = {
        "cluster": args.cluster,
        "block_size": args.block_size,
        "genomes": n,
        "blocks": nblocks,
        "expected_pairs": sum(expected.values()),
        "observed_pairs": in_cluster,
        "missing_pairs_vs_expected": max(0, sum(expected.values()) - in_cluster),
        "block_pairs_expected": len(expected),
        "block_pairs_complete": complete,
        "block_pairs_partial": partial,
        "block_pairs_missing": missing,
        "block_pairs_extra": extra,
        "table_rows_scanned": total,
        "cross_cluster_or_unmapped_rows_ignored": total - in_cluster,
        "invalid_in_cluster_ids": invalid_id,
    }
    report_tmp = Path(str(args.report) + ".tmp")
    report_tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report_tmp.replace(args.report)
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
