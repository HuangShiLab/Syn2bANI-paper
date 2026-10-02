#!/usr/bin/env python3
"""Restrict GTDB scope to genomes with non-empty digest TGTs.

Scanning the TGT directory once is much faster on Lustre than stat-ing two
files for every planned accession.
"""
import argparse
import os
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()
    root = Path(a.workdir)
    tgt = root / "tgt"
    valid = {}
    counts = Counter()
    with os.scandir(tgt) as it:
        for entry in it:
            if not entry.name.endswith(".tgt") or not entry.is_file(follow_symlinks=False):
                continue
            acc = entry.name[:-4]
            try:
                if entry.stat(follow_symlinks=False).st_size > 0:
                    valid[acc] = entry.path
                    counts["included"] += 1
            except FileNotFoundError:
                pass
    outdir = root / "effective_cluster_accessions"
    outdir.mkdir(exist_ok=True)
    expected = 0
    clusters = 0
    for list_path in (root / "cluster_accessions").glob("*.txt"):
        cluster = list_path.name[:-4]
        included = [x for x in list_path.read_text().split() if x in valid]
        if len(included) >= 2:
            clusters += 1
            expected += len(included) * (len(included) - 1) // 2
        (outdir / f"{cluster}.txt").write_text(
            "".join(x + "\n" for x in included))
    with open(root / "effective_scope_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"effective_genomes\t{len(valid)}\n")
        fh.write(f"effective_multi_genome_clusters\t{clusters}\n")
        fh.write(f"effective_unique_pairs\t{expected}\n")
    print(f"effective_genomes={len(valid)} clusters={clusters} pairs={expected}")


if __name__ == "__main__":
    main()
