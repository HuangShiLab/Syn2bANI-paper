#!/usr/bin/env python3
"""Merge one-to-all benchmark row files into a single TSV + summary.

Usage: python3 merge_one_to_all_rows.py <bench_base> <out_tsv>

Reads <bench_base>/results/rows_ani/*.tsv and rows_sv/*.tsv (one data row each,
headerless; columns per bench_gtdb_one_to_all.py HEADER), writes one TSV with
header, and prints a per-mode summary (median wall time, max RSS, n queries).
"""
import os
import sys
import glob
import statistics

HEADER = [
    "stage", "mode", "query_id", "ref_id", "n_refs", "threads", "tool", "tool_version",
    "exit_code", "wall_s", "timev_elapsed_s", "max_rss_mb_timev", "tree_rss_mb_polled",
    "n_rows_out", "kept_output", "host", "slurm_job", "slurm_array_task", "start_ts",
]


def main():
    base = sys.argv[1]
    out_tsv = sys.argv[2]
    rows = []
    for d in ("rows_ani", "rows_sv"):
        for path in sorted(glob.glob(os.path.join(base, "results", d, "*.tsv"))):
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == len(HEADER):
                        rows.append(parts)
                    else:
                        print(f"WARN: skipping malformed row in {path}: {len(parts)} fields")

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"wrote {len(rows)} rows -> {out_tsv}")

    by_mode = {}
    for r in rows:
        by_mode.setdefault(r[1], []).append(r)

    print(f"\n{'mode':16s} {'n':>3s} {'wall_s(med)':>12s} {'wall_s(range)':>24s} {'maxRSS_MB':>10s} {'rc!=0':>6s}")
    for mode in sorted(by_mode):
        rs = by_mode[mode]
        walls = sorted(float(r[9]) for r in rs)
        rss = []
        for r in rs:
            vals = []
            for col in (11, 12):
                try:
                    vals.append(float(r[col]))
                except ValueError:
                    pass
            rss.append(max(vals) if vals else 0.0)
        bad = sum(1 for r in rs if r[8] != "0")
        med = statistics.median(walls)
        print(f"{mode:16s} {len(rs):3d} {med:12.2f} {walls[0]:10.2f}..{walls[-1]:10.2f} {max(rss):10.1f} {bad:6d}")


if __name__ == "__main__":
    main()
