#!/usr/bin/env python3
"""Re-parse dnadiff 1-to-1 coords with a minimum gap threshold so that
small alignment fragmentation is ignored and only large rearrangements are
counted.

Usage:
    python3 parse_dnadiff_sv_filtered.py <min_gap_bp> <out_tsv>

Example:
    python3 parse_dnadiff_sv_filtered.py 5000 sv_truth_50k_min5k.tsv
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SYN2BANI_ROOT", HERE.parent.parent))
RES = ROOT / "results" / "gtdb50k"
OUTDIR = RES / "out"


def parse_one(pairid, min_gap):
    path = OUTDIR / pairid / "dd.1coords"
    if not path.exists():
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("="):
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            try:
                s1, e1, s2, e2 = map(int, parts[:4])
            except ValueError:
                continue
            rows.append((s1, e1, s2, e2))
    if len(rows) < 2:
        return dict(blocks=len(rows), breakpoints=0, large_indels=0)
    rows.sort(key=lambda x: x[0])
    breakpoints = 0
    large_indels = 0
    for i in range(len(rows) - 1):
        qs1, qe1, rs1, re1 = rows[i]
        qs2, qe2, rs2, re2 = rows[i + 1]
        qgap = qs2 - qe1 - 1
        strand1 = 1 if re1 >= rs1 else -1
        strand2 = 1 if re2 >= rs2 else -1
        if strand1 != strand2:
            breakpoints += 1
            # rgap is ill-defined across an inversion junction; only the
            # query-side gap can host an indel here.
            if qgap >= min_gap:
                large_indels += 1
        else:
            if strand1 == 1:
                rgap = rs2 - re1 - 1
            else:
                rgap = rs1 - re2 - 1
            if rgap < -min_gap:
                breakpoints += 1
            if qgap >= min_gap or abs(rgap) >= min_gap:
                large_indels += 1
    return dict(blocks=len(rows), breakpoints=breakpoints, large_indels=large_indels)


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 1
    min_gap = int(sys.argv[1])
    out_tsv = sys.argv[2]

    pair_dirs = [d for d in os.listdir(OUTDIR) if (OUTDIR / d).is_dir()]
    print(f"found {len(pair_dirs)} pair directories")
    records = []
    for i, pid in enumerate(sorted(pair_dirs)):
        if i % 5000 == 0 and i > 0:
            print(f"  processed {i}")
        rec = parse_one(pid, min_gap)
        if rec is None:
            continue
        rec["pairid"] = pid
        records.append(rec)

    df = pd.DataFrame(records)
    df = df.rename(columns={
        "blocks": f"dnadiff_blocks_min{min_gap}",
        "breakpoints": f"dnadiff_breakpoints_min{min_gap}",
        "large_indels": f"dnadiff_large_indels_min{min_gap}",
    })
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv} ({len(df)} pairs)")


if __name__ == "__main__":
    sys.exit(main())
