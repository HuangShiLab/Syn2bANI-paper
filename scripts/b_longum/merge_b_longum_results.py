#!/usr/bin/env python3
"""Merge B. longum Syn2bANI per-pair TSVs into one matrix."""
import glob
import os
from pathlib import Path

import pandas as pd

WORK = Path(os.environ.get("B_LONGUM_WORK", "results/b_longum_abfA"))
OUT = WORK / "s2b_out"


def main():
    files = sorted(glob.glob(str(OUT / "*.tsv")))
    print(f"found {len(files)} per-pair files")
    rows = []
    for f in files:
        pid = Path(f).stem
        qa, ra = pid.split("__", 1)
        df = pd.read_csv(f, sep="\t")
        if len(df) == 0:
            continue
        r = df.iloc[0].to_dict()
        r["query"] = qa
        r["reference"] = ra
        r["pairid"] = pid
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(WORK / "b_longum_s2b_matrix.tsv", sep="\t", index=False)
    print(f"wrote {WORK / 'b_longum_s2b_matrix.tsv'} ({len(out)} pairs)")


if __name__ == "__main__":
    main()
