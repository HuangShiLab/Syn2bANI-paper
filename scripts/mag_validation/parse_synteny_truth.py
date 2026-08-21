#!/usr/bin/env python3
"""parse_synteny_truth.py <binid> <dd.1coords> -> one TSV row on stdout:
bin \t n_aligns \t nucmer_breakpoints

Derives a rearrangement-breakpoint count from dnadiff 1-to-1 alignments.
Alignments are sorted by reference start; walking the sorted list, a break
is counted when two CONSECUTIVE alignments that land on the SAME query
contig violate collinearity (orientation flip, or query coordinate not
increasing for + orientation / not decreasing for -). Query-contig switches
are not breaks: on a multi-contig bin, contig order is arbitrary.

dnadiff .1coords layout (show-coords -THrcl: no header, whitespace columns):
rS rE qS qE len1 len2 %id lenR lenQ covR covQ tagR tagQ  (13 fields)
"""
import sys

binid, path = sys.argv[1], sys.argv[2]

rows = []
for line in open(path):
    f = line.split()
    if len(f) != 13 or not f[0].isdigit():
        continue
    try:
        rs, re_ = int(f[0]), int(f[1])
        qs, qe = int(f[2]), int(f[3])
    except ValueError:
        continue
    rows.append((min(rs, re_), qs, qe, f[12]))

rows.sort(key=lambda t: t[0])

breaks = 0
prev = None  # (qctg, orient, leading_qpos)
for rs, qs, qe, qctg in rows:
    orient = 1 if qe >= qs else -1
    qpos = qe if orient == 1 else qs  # leading edge along query
    if prev is not None and prev[0] == qctg:
        if orient != prev[1]:
            breaks += 1
        elif orient == 1 and qpos < prev[2]:
            breaks += 1
        elif orient == -1 and qpos > prev[2]:
            breaks += 1
    prev = (qctg, orient, qpos)

print(f"{binid}\t{len(rows)}\t{breaks}")
