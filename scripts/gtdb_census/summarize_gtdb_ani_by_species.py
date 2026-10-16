#!/usr/bin/env python3
"""Validate the merged GTDB ANI table and summarize it by species.

The merged table is expected to be globally sorted on the canonical
``genome_A``/``genome_B`` key and to contain one row for each effective
within-species pair. The script never holds pair-level data in memory.
"""

import argparse
import csv
import gzip
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani", required=True, help="merged sorted ANI TSV.gz")
    ap.add_argument("--workdir", required=True, help="GTDB census workdir")
    ap.add_argument("--quality", default="", help="optional genome quality table")
    ap.add_argument("--output", required=True, help="species summary TSV")
    ap.add_argument("--report", required=True, help="validation JSON")
    args = ap.parse_args()

    acc_species = {}
    expected = Counter()
    genomes = Counter()
    for path in sorted(Path(args.workdir, "effective_cluster_accessions").glob("*.txt")):
        species = path.name[:-4]
        accessions = path.read_text().split()
        expected[species] = len(accessions) * (len(accessions) - 1) // 2
        genomes[species] = len(accessions)
        for accession in accessions:
            acc_species[accession] = species

    taxonomy = {}
    if args.quality:
        with open(args.quality, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                taxonomy.setdefault(row["species_cluster"], row.get("gtdb_taxonomy", ""))

    stats = defaultdict(lambda: {
        "rows": 0, "ani_n": 0, "sum": 0.0, "min": math.inf,
        "max": -math.inf, "ge95": 0, "ge99": 0, "ge999": 0,
    })
    checks = {
        "header_ok": True,
        "sorted_canonical": True,
        "duplicate_keys": 0,
        "self_pairs": 0,
        "unknown_genomes": 0,
        "cross_species_pairs": 0,
        "invalid_ani": 0,
    }
    previous = None
    total_rows = 0
    opener = gzip.open if str(args.ani).endswith(".gz") else open
    with opener(args.ani, "rt") as fh:
        expected_header = ["genome_A", "genome_B", "skani_ani"]
        first_line = fh.readline().rstrip("\n")
        first = first_line.split("\t")
        if first == expected_header:
            checks["header_ok"] = True
            ix = {name: i for i, name in enumerate(expected_header)}
            lines = fh
        else:
            # Some sorted streams place the lowercase header after accession
            # keys; process the first line and continue scanning for fields.
            checks["header_ok"] = False
            ix = {name: i for i, name in enumerate(expected_header)}
            lines = itertools.chain((first_line,), fh)
        for line in lines:
            row = line.rstrip("\n").split("\t")
            if row[:3] == expected_header:
                continue
            if len(row) < 3:
                continue
            total_rows += 1
            a, b = row[ix["genome_A"]], row[ix["genome_B"]]
            value = fnum(row[ix["skani_ani"]])
            if a == b:
                checks["self_pairs"] += 1
            if a not in acc_species or b not in acc_species:
                checks["unknown_genomes"] += 1
            sa, sb = acc_species.get(a, ""), acc_species.get(b, "")
            if sa and sb and sa != sb:
                checks["cross_species_pairs"] += 1
            key = (a, b)
            if previous is not None and key < previous:
                checks["sorted_canonical"] = False
            elif previous == key:
                checks["duplicate_keys"] += 1
            previous = key
            if value is None or not math.isfinite(value):
                checks["invalid_ani"] += 1
                continue
            if sa:
                st = stats[sa]
                st["rows"] += 1
                st["ani_n"] += 1
                st["sum"] += value
                st["min"] = min(st["min"], value)
                st["max"] = max(st["max"], value)
                st["ge95"] += value >= 95
                st["ge99"] += value >= 99
                st["ge999"] += value >= 99.9

    fields = [
        "species_cluster", "taxonomy", "n_genomes", "expected_pairs",
        "ani_pairs", "ani_coverage", "mean_skani_ani", "min_skani_ani",
        "max_skani_ani", "frac_ani_ge95", "frac_ani_ge99", "frac_ani_ge999",
        "ani_layer",
    ]
    tmp = Path(str(args.output) + ".tmp")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for species in sorted(set(expected) | set(stats)):
            st = stats[species]
            n = st["ani_n"]
            exp = expected[species]
            writer.writerow({
                "species_cluster": species,
                "taxonomy": taxonomy.get(species, ""),
                "n_genomes": genomes[species],
                "expected_pairs": exp,
                "ani_pairs": st["rows"],
                "ani_coverage": st["rows"] / exp if exp else "",
                "mean_skani_ani": st["sum"] / n if n else "",
                "min_skani_ani": st["min"] if n else "",
                "max_skani_ani": st["max"] if n else "",
                "frac_ani_ge95": st["ge95"] / n if n else "",
                "frac_ani_ge99": st["ge99"] / n if n else "",
                "frac_ani_ge999": st["ge999"] / n if n else "",
                "ani_layer": "complete" if exp and st["rows"] == exp else "partial",
            })
    tmp.replace(args.output)

    report = {
        "ani_table": str(args.ani),
        "total_rows": total_rows,
        "effective_genomes": len(acc_species),
        "effective_species": len(expected),
        "expected_unique_pairs": sum(expected.values()),
        "reported_unique_pairs": total_rows,
        "reported_within_expected_species": sum(st["rows"] for st in stats.values()),
        "missing_pairs_vs_effective_scope": sum(expected.values()) - total_rows,
        **checks,
    }
    report_tmp = Path(str(args.report) + ".tmp")
    report_tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report_tmp.replace(args.report)
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
