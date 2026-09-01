#!/usr/bin/env python3
"""Extract dnadiff's classified rearrangement events from existing dd.report files.

`compute_dnadiff_inverted_fraction.py` reduces each pair to one number, which
validates only the magnitude of Syn2b's orientation channel. dnadiff's report
already separates the events by *kind*:

    [Feature Estimates]
    Breakpoints        <ref>  <qry>
    Relocations        <ref>  <qry>     within-replicon movement
    Translocations     <ref>  <qry>     between-replicon movement
    Inversions         <ref>  <qry>     orientation flips
    Insertions         <ref>  <qry>     indels, NOT rearrangements
    ...

That separation is what the junction/count channel has to be tested against, and
it also decides whether the +290 intercept in `dnadiff_breakpoints ~
breakpoint_count` (SV_REANALYSIS.md) is rearrangement signal or indel bookkeeping.

No new alignment is needed: dd.report is already on disk for every pair, written
by scripts/gtdb50k/run_dnadiff_slice.sh. This is a parse pass, minutes not hours.

Usage:
    python3 scripts/compute_dnadiff_events.py \
        --pairs results/gtdb50k/pairs_50k.tsv \
        --outdir results/gtdb50k/out \
        --outfile results/gtdb50k/dnadiff_events_50k.tsv
"""

import argparse
import re
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "results" / "gtdb50k"

FEATURES = [
    "Breakpoints",
    "Relocations",
    "Translocations",
    "Inversions",
    "Insertions",
    "InsertionSum",
    "TandemIns",
    "TandemInsSum",
]
# dnadiff prints the ref column first, then the qry column.
ROW = re.compile(r"^(\w+)\s+(\S+)\s+(\S+)\s*$", re.M)


def parse_one(args):
    pairid, outdir = args
    path = Path(outdir) / pairid / "dd.report"
    if not path.exists():
        return None
    txt = path.read_text()

    block = txt.split("[Feature Estimates]", 1)
    if len(block) < 2:
        return None

    rec = {"pairid": pairid}
    for name, ref, qry in ROW.findall(block[1]):
        if name not in FEATURES:
            continue
        for side, raw in (("ref", ref), ("qry", qry)):
            try:
                rec[f"dd_{name.lower()}_{side}"] = float(raw)
            except ValueError:
                pass
    return rec if len(rec) > 1 else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", default=str(RES / "pairs_50k.tsv"))
    p.add_argument("--outdir", default=str(RES / "out"))
    p.add_argument("--outfile", default=str(RES / "dnadiff_events_50k.tsv"))
    p.add_argument("--workers", type=int, default=min(32, cpu_count() or 1))
    args = p.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    if "pairid" not in pairs.columns:
        pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]

    tasks = [(pid, args.outdir) for pid in pairs["pairid"].tolist()]
    print(
        f"parsing {len(tasks)} dd.report files from {args.outdir} "
        f"with {args.workers} workers ...",
        flush=True,
    )

    records = []
    with Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(parse_one, tasks, chunksize=200)):
            if i % 5000 == 0 and i > 0:
                print(f"  processed {i}/{len(tasks)}", flush=True)
            if rec is not None:
                records.append(rec)

    df = pd.DataFrame(records)
    keep = [c for c in ["pairid", "band", "phylum"] if c in pairs.columns]
    df = pairs[keep].merge(df, on="pairid", how="left")
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"wrote {args.outfile}: {len(df)} pairs, {df.filter(like='dd_').notna().all(axis=1).sum()} complete")

    cols = [c for c in df.columns if c.startswith("dd_")]
    if cols:
        print("\nmedian per feature:")
        print(df[cols].median().round(2).to_string())


if __name__ == "__main__":
    main()
