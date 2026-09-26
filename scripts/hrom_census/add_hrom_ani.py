#!/usr/bin/env python3
"""Attach cluster-wise skani ANI to the HROM structural census table.

The structural merge writes rows grouped by cluster because outputs are
processed in sorted filename order. This script exploits that ordering and
loads one cluster's ANI triangle at a time, so it does not keep all 62.9M
pairwise ANI values in memory.
"""
import argparse
import csv
import gzip
import json
import random
from collections import Counter
from pathlib import Path

BINS = [
    ("below_95", lambda x: x < 95.0),
    ("95_97", lambda x: 95.0 <= x < 97.0),
    ("97_98", lambda x: 97.0 <= x < 98.0),
    ("98_99", lambda x: 98.0 <= x < 99.0),
    ("99_99.5", lambda x: 99.0 <= x < 99.5),
    ("ge_99.5", lambda x: x >= 99.5),
]


def median(values):
    values = sorted(values)
    if not values:
        return ""
    n = len(values)
    mid = n // 2
    return values[mid] if n % 2 else (values[mid - 1] + values[mid]) / 2.0


def load_cluster_ani(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[(row["genome_A"], row["genome_B"])] = float(row["skani_ani"])
    return out


class Reservoir:
    def __init__(self, size=100000, seed=42):
        self.size = size
        self.values = []
        self.n = 0
        self.rng = random.Random(seed)

    def add(self, value):
        self.n += 1
        if len(self.values) < self.size:
            self.values.append(value)
        else:
            j = self.rng.randrange(self.n)
            if j < self.size:
                self.values[j] = value

    def median(self):
        return median(self.values)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--metadata", required=True,
                    help="hrom_genome_metadata.tsv from prepare_hrom_census.py")
    ap.add_argument("--structural", required=True,
                    help="merged hrom_within_species_sv.tsv.gz")
    ap.add_argument("--out", required=True,
                    help="output gz with appended skani_ani column")
    ap.add_argument("--summary", required=True, help="ANI-bin summary TSV")
    ap.add_argument("--qc", required=True, help="merge QC report")
    ap.add_argument("--max-missing", type=int, default=20,
                    help="number of missing-ANI examples to record")
    args = ap.parse_args()

    work = Path(args.workdir)
    genome_cluster = {}
    with open(args.metadata) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_cluster[row["genome"]] = row["cluster"]

    current_cluster = None
    ani = {}
    missing = []
    matched = 0
    total = 0
    bin_stats = {
        name: {"n": 0, "ge2": 0, "junction_res": Reservoir(),
               "inverted_res": Reservoir(), "sum_junctions": 0.0}
        for name, _ in BINS
    }

    with gzip.open(args.structural, "rt") as fin, \
         gzip.open(args.out, "wt") as fout:
        header = fin.readline().rstrip("\n")
        if not header:
            raise ValueError("structural table is empty")
        cols = header.split("\t")
        if "genome_A" not in cols or "genome_B" not in cols:
            raise ValueError("structural table lacks genome_A/genome_B")
        if "breakpoints" not in cols or "raw_inverted_fraction" not in cols:
            raise ValueError("structural table lacks breakpoint columns")
        out_header = cols + ["skani_ani"]
        fout.write("\t".join(out_header) + "\n")
        ia, ib = cols.index("genome_A"), cols.index("genome_B")
        ibp = cols.index("breakpoints")
        iinv = cols.index("raw_inverted_fraction")

        for line in fin:
            total += 1
            fields = line.rstrip("\n").split("\t")
            ga, gb = fields[ia], fields[ib]
            cl = genome_cluster.get(ga)
            if cl is None or genome_cluster.get(gb) != cl:
                missing.append(f"{ga}|{gb}: non-cluster pair")
                fields.append("NA")
            else:
                if cl != current_cluster:
                    ani = load_cluster_ani(work / "ani" / f"{cl}.ani.tsv.gz")
                    current_cluster = cl
                akey = (ga, gb) if (ga, gb) in ani else (gb, ga)
                value = ani.get(akey)
                if value is None:
                    missing.append(f"{ga}|{gb}: missing ANI")
                    fields.append("NA")
                else:
                    fields.append(f"{value:.6f}")
                    matched += 1
                    try:
                        junctions = float(fields[ibp])
                        inverted = float(fields[iinv])
                    except ValueError:
                        junctions = inverted = None
                    if junctions is not None:
                        for name, predicate in BINS:
                            if predicate(value):
                                s = bin_stats[name]
                                s["n"] += 1
                                s["ge2"] += junctions >= 2
                                s["sum_junctions"] += junctions
                                s["junction_res"].add(junctions)
                                s["inverted_res"].add(inverted)
                                break
            if len(missing) <= args.max_missing:
                pass
            fout.write("\t".join(fields) + "\n")

    with open(args.summary, "w") as fh:
        fh.write("ani_bin\tn\tpct_ge2_junctions\tmedian_junctions\t"
                 "mean_junctions\tmedian_raw_inverted_fraction\n")
        for name, _ in BINS:
            s = bin_stats[name]
            if s["n"]:
                pct = 100.0 * s["ge2"] / s["n"]
                mean = s["sum_junctions"] / s["n"]
                fh.write(f"{name}\t{s['n']}\t{pct:.6g}\t"
                         f"{s['junction_res'].median():.6g}\t{mean:.6g}\t"
                         f"{s['inverted_res'].median():.6g}\n")
            else:
                fh.write(f"{name}\t0\tNA\tNA\tNA\tNA\n")

    with open(args.qc, "w") as fh:
        fh.write("# HROM structural x ANI merge QC\n\n")
        fh.write(f"- structural rows read: {total:,}\n")
        fh.write(f"- rows with skani ANI: {matched:,}\n")
        fh.write(f"- rows without skani ANI: {total - matched:,}\n")
        if missing:
            fh.write("\nExamples:\n")
            for item in missing[:args.max_missing]:
                fh.write(f"- {item}\n")
    print(f"rows={total:,} matched_ani={matched:,} missing={total - matched:,}")


if __name__ == "__main__":
    main()
