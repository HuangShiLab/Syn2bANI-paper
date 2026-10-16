#!/usr/bin/env python3
"""Canonicalize a completed skani triangle raw table without recomputation."""

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--pigz-threads", type=int, default=4)
    ap.add_argument("--delete-raw", action="store_true")
    args = ap.parse_args()

    accs = Path(args.accessions).read_text().split()
    expected = len(accs) * (len(accs) - 1) // 2
    acc_set = set(accs)
    tmp = Path(str(args.output) + ".tmp")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    pigz = subprocess.Popen(
        ["pigz", "-p", str(args.pigz_threads)], stdin=subprocess.PIPE,
        stdout=tmp.open("wb"), text=True,
    )
    rows = invalid = self_pairs = unexpected_accessions = 0
    with open(args.raw) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"Ref_file", "Query_file", "ANI"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"unexpected triangle columns: {reader.fieldnames}")
        pigz.stdin.write("genome_A\tgenome_B\tskani_ani\n")
        for row in reader:
            aa = Path(row["Ref_file"]).stem
            bb = Path(row["Query_file"]).stem
            try:
                ani = float(row["ANI"])
            except (TypeError, ValueError):
                invalid += 1
                continue
            if not math.isfinite(ani):
                invalid += 1
                continue
            if aa == bb:
                self_pairs += 1
                continue
            if aa not in acc_set or bb not in acc_set:
                unexpected_accessions += 1
                continue
            if bb < aa:
                aa, bb = bb, aa
            pigz.stdin.write(f"{aa}\t{bb}\t{row['ANI']}\n")
            rows += 1
    pigz.stdin.close()
    pigz.wait()
    if pigz.returncode:
        raise SystemExit(f"pigz failed with code {pigz.returncode}")
    tmp.replace(args.output)

    report = {
        "cluster": args.cluster,
        "raw": args.raw,
        "output": args.output,
        "input_accessions": len(accs),
        "expected_pairs": expected,
        "reported_pairs": rows,
        "missing_pairs": max(0, expected - rows),
        "invalid_ani": invalid,
        "self_pairs": self_pairs,
        "unexpected_accessions": unexpected_accessions,
        "complete": rows == expected,
    }
    report_tmp = Path(str(args.report) + ".tmp")
    report_tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report_tmp.replace(args.report)
    if args.delete_raw and report["complete"]:
        Path(args.raw).unlink()
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
