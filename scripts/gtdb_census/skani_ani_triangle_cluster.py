#!/usr/bin/env python3
"""Run a complete skani triangle ANI pass for one GTDB species cluster."""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--pigz-threads", type=int, default=4)
    ap.add_argument("--min-af", type=float, default=0.0)
    ap.add_argument("--screen", type=float, default=0.0)
    args = ap.parse_args()

    root = Path(args.workdir)
    safe = args.cluster.replace("/", "_")
    accs = Path(root, "effective_cluster_accessions", safe + ".txt").read_text().split()
    manifest = json.loads((root / "manifest.json").read_text())
    missing_fasta = [acc for acc in accs if not manifest.get(acc) or not Path(manifest[acc]).is_file()]
    paths = [manifest[acc] for acc in accs if acc not in set(missing_fasta)]
    if len(paths) < 2:
        raise SystemExit(f"cluster has fewer than two usable FASTAs: {args.cluster}")

    work = Path(args.tmpdir)
    work.mkdir(parents=True, exist_ok=True)
    list_path = work / f"{safe}.paths.txt"
    list_path.write_text("\n".join(paths) + "\n")
    raw = work / f"{safe}.triangle.raw"
    err = work / f"{safe}.triangle.err"
    cmd = [
        args.skani, "triangle", "-E", "-t", str(args.threads),
        "--min-af", str(args.min_af), "-s", str(args.screen),
        "-l", str(list_path), "-o", str(raw),
    ]
    with err.open("w") as ef:
        subprocess.run(cmd, stdout=ef, stderr=ef, check=True)

    tmp_out = Path(str(args.output) + ".tmp")
    tmp_out.parent.mkdir(parents=True, exist_ok=True)
    pigz = subprocess.Popen(
        ["pigz", "-p", str(args.pigz_threads)], stdin=subprocess.PIPE,
        stdout=tmp_out.open("wb"), text=True,
    )
    expected = len(accs) * (len(accs) - 1) // 2
    rows = invalid = self_pairs = 0
    with raw.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"Ref_file", "Query_file", "ANI"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"unexpected triangle columns: {reader.fieldnames}")
        pigz.write("genome_A\tgenome_B\tskani_ani\n")
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
            if bb < aa:
                aa, bb = bb, aa
            pigz.write(f"{aa}\t{bb}\t{row['ANI']}\n")
            rows += 1
    pigz.close()
    if pigz.returncode:
        raise SystemExit(f"pigz failed with code {pigz.returncode}")
    tmp_out.replace(args.output)
    raw.unlink(missing_ok=True)
    list_path.unlink(missing_ok=True)

    report = {
        "cluster": args.cluster,
        "expected_pairs": expected,
        "reported_pairs": rows,
        "missing_pairs": max(0, expected - rows),
        "invalid_ani": invalid,
        "self_pairs": self_pairs,
        "input_accessions": len(accs),
        "missing_or_absent_fasta": len(missing_fasta),
        "command": " ".join(cmd),
        "output": str(args.output),
    }
    report_path = Path(str(args.output) + ".report.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True), flush=True)
    if rows != expected:
        print(
            f"WARNING: {args.cluster} ANI coverage {rows}/{expected}",
            file=sys.stderr, flush=True,
        )


if __name__ == "__main__":
    main()
