#!/usr/bin/env python3
"""parse_repsearch.py <skani_search.tsv> > map.tsv
Pick nearest GTDB rep per bin: among top hits, prefer highest ANI with
Align_fraction_query >= 10%; if none qualify, fall back to the best-AF hit.
Output: bin \t rep_path \t skani_ani \t skani_af_query
skani 0.3.1 columns: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query [Ref_name Query_name]
"""
import sys

rows = {}
with open(sys.argv[1]) as f:
    header = f.readline().rstrip("\n").split("\t")
    for line in f:
        d = dict(zip(header, line.rstrip("\n").split("\t")))
        q = d.get("Query_file", "")
        if not q:
            continue
        binid = q.rsplit("/", 1)[-1][:-3] if q.endswith(".fa") else q
        ref = d["Ref_file"]
        ani = float(d["ANI"])
        afq = float(d.get("Align_fraction_query", 0))
        cur = rows.setdefault(binid, {"best": None, "fallback": None})
        if afq >= 10.0 and (cur["best"] is None or ani > cur["best"][1]):
            cur["best"] = (ref, ani, afq)
        if cur["fallback"] is None or afq > cur["fallback"][2]:
            cur["fallback"] = (ref, ani, afq)

sys.stdout.write("bin\trep_path\tskani_ani\tskani_af_query\n")
for binid in sorted(rows):
    pick = rows[binid]["best"] or rows[binid]["fallback"]
    if pick:
        sys.stdout.write(f"{binid}\t{pick[0]}\t{pick[1]:.3f}\t{pick[2]:.2f}\n")
