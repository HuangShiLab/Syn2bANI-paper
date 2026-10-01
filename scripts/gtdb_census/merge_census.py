#!/usr/bin/env python3
"""Merge census task outputs into a per-pair table and low-memory QC tables.

  python3 merge_census.py --workdir ... --out results/census
"""
import argparse
import gzip
import random
import re
from collections import defaultdict
from pathlib import Path

NUM = re.compile(r"^-?\d+(\.\d+)?([eE][-+]?\d+)?$")


def median(values):
    values = sorted(values)
    if not values:
        return ""
    n = len(values)
    mid = n // 2
    return values[mid] if n % 2 else (values[mid - 1] + values[mid]) / 2.0


def add_sample(state, key, value, cap):
    state["seen"] += 1
    sample = state[key]
    if len(sample) < cap:
        sample.append(value)
    else:
        j = state["rng"].randrange(state["seen"])
        if j < cap:
            sample[j] = value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--output-name",
                    default="gtdb_r207_within_species_sv.tsv.gz",
                    help="final compressed table filename")
    ap.add_argument("--summary-sample-size", type=int, default=1000,
                    help="per-cluster reservoir size for medians")
    a = ap.parse_args()

    root, out = Path(a.workdir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    final = out / a.output_name

    # GTDB has ~711M pairs.  Keep exact counts and exact P(>=2 breakpoints),
    # but use reservoirs for medians so the merger stays low-memory.
    cap = max(1, a.summary_sample_size)
    empty = lambda: {
        "rows": 0, "ge2": 0, "bp": [], "inv": [], "seen": 0,
        "rng": random.Random(0x534e32),
    }
    stats = defaultdict(empty)
    total_rows = 0
    duplicate_local = 0

    cluster_sizes = {}
    expected = {}
    for f in (root / "cluster_accessions").glob("*.txt"):
        cl = f.name[:-4]
        n = sum(1 for _ in open(f))
        cluster_sizes[cl] = n
        expected[cl] = n * (n - 1) // 2

    header = None
    with gzip.open(final, "wt") as fout:
        for gz in sorted((root / "outputs").glob("*.tsv.gz")):
            pair_cluster = gz.name.split("|", 1)[0]
            with gzip.open(gz, "rt") as fin:
                first = fin.readline().rstrip("\n")
                if not first:
                    continue
                if header is None:
                    header = first.split("\t")
                    fout.write(first + "\n")
                    cols = {name: i for i, name in enumerate(header)}
                    for required in ("genome_A", "genome_B", "breakpoints",
                                     "raw_inverted_fraction"):
                        if required not in cols:
                            raise ValueError(f"missing column {required}")

                # Duplicate checks are local to an output: the worker design
                # partitions pairs by query batch, so a global pair dictionary
                # is unnecessary and would not fit in memory for 711M pairs.
                local_seen = set()
                for line in fin:
                    f = line.rstrip("\n").split("\t")
                    a, b = f[cols["genome_A"]], f[cols["genome_B"]]
                    key = (a, b)
                    if key in local_seen:
                        duplicate_local += 1
                        continue
                    local_seen.add(key)
                    fout.write(line)
                    total_rows += 1

                    s = stats[pair_cluster]
                    s["rows"] += 1
                    try:
                        bp = float(f[cols["breakpoints"]])
                    except ValueError:
                        bp = None
                    try:
                        inv = float(f[cols["raw_inverted_fraction"]])
                    except ValueError:
                        inv = None
                    if bp is not None:
                        s["ge2"] += int(bp >= 2)
                        add_sample(s, "bp", bp, cap)
                    if inv is not None:
                        add_sample(s, "inv", inv, cap)

    with open(out / "species_summary.tsv", "w") as fh:
        fh.write("cluster\tgenomes\texpected_pairs\trows\tmedian_junctions\t"
                 "pct_pairs_ge2_junctions\tmedian_raw_inverted_fraction\n")
        for cl in sorted(stats):
            s = stats[cl]
            rows = s["rows"]
            pct = 100.0 * s["ge2"] / rows if rows else ""
            fh.write(f"{cl}\t{cluster_sizes.get(cl, 0)}\t{expected.get(cl, 0)}\t{rows}\t"
                     f"{median(s['bp'])}\t{pct}\t{median(s['inv'])}\n")

    expected_total = sum(expected.values())
    short = [(c, expected[c], stats[c]["rows"]) for c in expected
             if c in stats and stats[c]["rows"] < expected[c]]
    with open(out / "census_qc.md", "w") as fh:
        fh.write("# Census QC\n\n")
        fh.write(f"- unique pairs written: {total_rows:,}\n")
        fh.write(f"- expected pairs: {expected_total:,}\n")
        fh.write(f"- local duplicate rows suppressed: {duplicate_local:,}\n")
        fh.write(f"- clusters with fewer rows than C(n,2): {len(short):,}\n")
        for c, e, r in sorted(short, key=lambda x: x[1] - x[2], reverse=True)[:20]:
            fh.write(f"  - {c}: expected {e:,} got {r:,}\n")
    print(f"wrote {final} ({total_rows:,} pairs); QC in census_qc.md")


if __name__ == "__main__":
    main()
