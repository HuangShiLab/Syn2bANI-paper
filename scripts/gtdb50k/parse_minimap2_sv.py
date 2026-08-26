#!/usr/bin/env python3
"""Parse a minimap2 PAF file and emit SV metrics for one pair.

Usage:
    python3 parse_minimap2_sv.py <pairid> <mm2.paf>

Output: one TSV line
    pairid, mm2_blocks, mm2_breakpoints, mm2_large_indels, mm2_synteny_score, status
"""
import sys
from pathlib import Path

MIN_INDEL = 1000


def parse_paf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            try:
                qs, qe = int(parts[2]), int(parts[3])
                strand = parts[4]
                rs, re = int(parts[7]), int(parts[8])
            except (ValueError, IndexError):
                continue
            rows.append((qs, qe, rs, re, strand))
    return rows


def metrics(rows):
    if not rows:
        return dict(blocks=0, breakpoints=0, large_indels=0, synteny_score=1.0)
    rows = sorted(rows, key=lambda x: x[0])
    n = len(rows)
    if n == 1:
        return dict(blocks=1, breakpoints=0, large_indels=0, synteny_score=1.0)
    breakpoints = 0
    large_indels = 0
    for i in range(n - 1):
        qs1, qe1, rs1, re1, st1 = rows[i]
        qs2, qe2, rs2, re2, st2 = rows[i + 1]
        qgap = qs2 - qe1 - 1
        if st1 != st2:
            breakpoints += 1
        else:
            if st1 == "+":
                rgap = rs2 - re1 - 1
                if rgap < -MIN_INDEL:
                    breakpoints += 1
            else:
                rgap = rs1 - re2 - 1
                if rgap < -MIN_INDEL:
                    breakpoints += 1
        if qgap >= MIN_INDEL or abs(rgap) >= MIN_INDEL:
            large_indels += 1
    score = max(0.0, 1.0 - breakpoints / (n - 1)) if n > 1 else 1.0
    return dict(blocks=n, breakpoints=breakpoints, large_indels=large_indels,
                synteny_score=score)


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 1
    pairid, path = sys.argv[1], sys.argv[2]
    if not Path(path).exists():
        print(f"{pairid}\tNA\tNA\tNA\tmissing_paf")
        return 0
    rows = parse_paf(path)
    m = metrics(rows)
    print(f"{pairid}\t{m['blocks']}\t{m['breakpoints']}\t{m['large_indels']}\t{m['synteny_score']:.6f}\tok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
