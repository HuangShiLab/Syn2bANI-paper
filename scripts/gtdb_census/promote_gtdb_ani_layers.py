#!/usr/bin/env python3
"""Promote completed GTDB ANI chunk norms/triangle outputs for resume."""

import argparse
import gzip
import json
import os
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    root = Path(args.workdir)
    normdir = root / "ani_chunk_norm"
    triangledir = root / "ani_triangle"
    anidir = root / "ani"
    anidir.mkdir(parents=True, exist_ok=True)

    promoted = []
    skipped_existing = []
    empty_norms = []
    # Link complete triangle outputs first. This makes the four large clusters
    # visible before legacy chunk norms are inspected.
    linked = []
    for src in sorted(triangledir.glob("*.ani.tsv.gz")):
        final = anidir / src.name
        if final.exists():
            skipped_existing.append(src.name)
            continue
        try:
            os.link(src, final)
        except OSError:
            final.write_bytes(src.read_bytes())
        linked.append({"source": str(src), "output": str(final)})

    invalid_norms = []
    for src in sorted(normdir.glob("*.norm")):
        final = anidir / (src.name[:-5] + ".ani.tsv.gz")
        if final.exists():
            skipped_existing.append(src.name)
            continue
        with src.open("rb") as probe:
            magic = probe.read(2)
        if src.stat().st_size == 0:
            empty_norms.append(src.name)
            continue
        if magic == b"\x1f\x8b":
            inp_ctx = gzip.open(src, "rt")
        elif magic == b"ge":
            # Earlier chunk workers wrote uncompressed TSV but kept `.norm`.
            inp_ctx = src.open("rt")
        else:
            invalid_norms.append({"source": str(src), "magic": magic.hex()})
            continue
        tmp = Path(str(final) + ".tmp")
        rows = 0
        with inp_ctx as inp, gzip.open(tmp, "wt") as out:
            out.write("genome_A\tgenome_B\tskani_ani\n")
            for line in inp:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                aa, bb, ani = fields[:3]
                if aa == bb:
                    continue
                if bb < aa:
                    aa, bb = bb, aa
                out.write(f"{aa}\t{bb}\t{ani}\n")
                rows += 1
        tmp.replace(final)
        promoted.append({"source": str(src), "output": str(final), "rows": rows})

    report = {
        "norm_promoted": len(promoted),
        "norm_promoted_rows": sum(x["rows"] for x in promoted),
        "triangle_linked": len(linked),
        "skipped_existing": len(skipped_existing),
        "empty_norm_left_for_retry": empty_norms,
        "invalid_norms_left_unpromoted": invalid_norms,
        "promoted": promoted,
        "linked": linked,
    }
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    summary = {
            key: report[key] for key in (
            "norm_promoted", "norm_promoted_rows", "triangle_linked",
            "skipped_existing", "empty_norm_left_for_retry",
            "invalid_norms_left_unpromoted",
        )
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
