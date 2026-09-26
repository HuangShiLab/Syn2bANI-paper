#!/usr/bin/env python3
"""Merge census task outputs into the final per-pair table + QC + species summary.

  python3 merge_census.py --workdir ... --out results/census

Outputs:
  gtdb_r207_within_species_sv.tsv.gz   the census table
  census_qc.md                         expected-vs-actual rows per cluster
  species_summary.tsv                  per-cluster junction statistics
"""
import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

NUM = re.compile(r"^-?\d+(\.\d+)?([eE][-+]?\d+)?$")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--output-name",
                    default="gtdb_r207_within_species_sv.tsv.gz",
                    help="final compressed table filename")
    a = ap.parse_args()
    root, out = Path(a.workdir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    final = out / a.output_name
    stats = defaultdict(lambda: [0, [], []])   # rows, junctions, inverted

    header = None
    dup_guard = defaultdict(int)
    cluster_sizes = {}
    expected = {}
    for f in (root / "cluster_accessions").glob("*.txt"):
        cl = f.name[:-4]
        n = sum(1 for _ in open(f))
        cluster_sizes[cl] = n
        expected[cl] = n * (n - 1) // 2
    with gzip.open(final, "wt") as fout:
        for gz in sorted((root / "outputs").glob("*.tsv.gz")):
            with gzip.open(gz, "rt") as fin:
                first = fin.readline().rstrip("\n")
                if header is None:
                    header = first.split("\t")
                    fout.write(first + "\n")
                    cols = {name: i for i, name in enumerate(header)}
                for line in fin:
                    f = line.rstrip("\n").split("\t")
                    a = f[cols["genome_A"]]
                    b = f[cols["genome_B"]]
                    key = (a, b)
                    dup_guard[key] += 1
                    if dup_guard[key] > 1:
                        continue
                    fout.write(line)
                    pair_cluster = gz.name.split("|", 1)[0]
                    s = stats[pair_cluster]
                    s[0] += 1
                    for col, slot in (("breakpoints", 1), ("raw_inverted_fraction", 2)):
                        v = f[cols[col]]
                        if v and re.fullmatch(r"-?\d+(\.\d+)?([eE][-+]?\d+)?", v):
                            s[slot].append(float(v))

    with open(out / "species_summary.tsv", "w") as fh:
        fh.write("cluster\tgenomes\texpected_pairs\trows\tmedian_junctions\t"
                 "pct_pairs_ge2_junctions\tmedian_raw_inverted_fraction\n")
        for cl in sorted(stats):
            rows, js, iv = stats[cl]
            js.sort()
            fh.write(f"{cl}\t{cluster_sizes.get(cl, 0)}\t{expected.get(cl, 0)}\t{rows}\t"
                     f"{js[len(js)//2] if js else ''}\t"
                     f"{100*sum(1 for j in js if j >= 2)/len(js) if js else ''}\t"
                     f"{sorted(iv)[len(iv)//2] if iv else ''}\n")

    dup = sum(1 for v in dup_guard.values() if v > 1)
    short = [(c, expected[c], stats[c][0]) for c in expected
             if c in stats and stats[c][0] < expected[c]]
    with open(out / "census_qc.md", "w") as fh:
        fh.write("# Census QC\n\n")
        fh.write(f"- unique pairs written: {len(dup_guard):,}\n")
        fh.write(f"- duplicate rows suppressed: {dup:,}\n")
        fh.write(f"- clusters with fewer rows than C(n,2) (missing genomes): "
                 f"{len(short):,}\n")
        for c, e, r in sorted(short, key=lambda x: x[1]-x[2], reverse=True)[:20]:
            fh.write(f"  - {c}: expected {e:,} got {r:,}\n")
    print(f"wrote {final} ({len(dup_guard):,} unique pairs); QC in census_qc.md")


if __name__ == "__main__":
    main()
