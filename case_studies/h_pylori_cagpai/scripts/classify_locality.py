#!/usr/bin/env python3
"""Split `complete_rearranged` by where the structural call actually breaks.

The extended state calls a genome `complete_rearranged` when any inversion or
translocation call overlaps the cagPAI window (NC_000915.1: 547,327-583,481,
+/- 2 kb). Overlap is satisfied by a 100-600 kb inversion that contains the
whole island; in that case the island is intact and merely lies on the other
strand of a larger block. This script separates three situations:

  island_internal   at least one call lies entirely inside the window
  island_boundary   at least one call has an endpoint inside the window
                    (the island itself is broken) but none lies entirely inside
  island_spanned    every overlapping call spans the whole window; the island
                    is intact inside a larger inverted/translocated block

Input : results/cagpai_states_extended_filtered.tsv
Output: results/cagpai_states_locality.tsv, counts on stdout
"""
import csv
import re
from collections import Counter
from pathlib import Path

WORK = Path(__file__).resolve().parent.parent
INP = WORK / "results" / "cagpai_states_extended_filtered.tsv"
OUT = WORK / "results" / "cagpai_states_locality.tsv"
LO, HI = 547_327 - 2_000, 583_481 + 2_000


def spans(s):
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+)-(\d+)", s or "")]


def locality(sv_list):
    calls = spans(sv_list)
    if any(LO <= a and b <= HI for a, b in calls):
        return "island_internal"
    if any(LO <= a <= HI or LO <= b <= HI for a, b in calls):
        return "island_boundary"
    return "island_spanned"


def main():
    rows = list(csv.DictReader(open(INP), delimiter="\t"))
    counts = Counter()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "status_extended", "locality", "cag_sv_list"])
        for r in rows:
            st = r["status_extended"]
            loc = locality(r["cag_sv_list"]) if st == "complete_rearranged" else ""
            counts[(st, loc)] += 1
            w.writerow([r["genome"], st, loc, r["cag_sv_list"]])
    for (st, loc), n in sorted(counts.items()):
        print(f"{st:22s} {loc:16s} {n}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
