#!/usr/bin/env python3
"""Collect ANIm truth and run syn2bani panel --greedy per ANI band.

Usage:
    python3 run_panel_per_band.py \
        --sample results/sample_anim_truth.tsv \
        --anim-dir anim_results \
        --strata results/gtdb_r207_100k_strata.tsv \
        --syn2bani /lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani \
        --outdir results/panel_by_band
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def read_sample(path):
    """Return dict (query,ref) -> band."""
    bands = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, r, band = f[0], f[1], f[2]
            bands[(q, r)] = band
    return bands


def read_anim_truth(anim_dir):
    """Collect dnadiff results, skipping NA/failed."""
    truth = {}
    for p in Path(anim_dir).glob("anim_*.tsv"):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                q, r, ident = f[0], f[1], f[2]
                try:
                    v = float(ident)
                except ValueError:
                    continue
                if not v > 0:
                    continue
                truth[(q, r)] = v
    return truth


def split_truth_by_band(truth, bands, outdir):
    """Write per-band truth files plus a pooled one."""
    per_band = defaultdict(list)
    all_lines = []
    for (q, r), v in truth.items():
        band = bands.get((q, r))
        if band is None:
            continue
        line = f"{q}\t{r}\t{v:.4f}\n"
        all_lines.append(line)
        per_band[band].append(line)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pooled = outdir / "truth_all.tsv"
    with open(pooled, "w") as fh:
        fh.writelines(all_lines)
    print(f"wrote {pooled}: {len(all_lines)} pairs")

    band_paths = {}
    for band, lines in sorted(per_band.items()):
        p = outdir / f"truth_{band.replace('.', '_')}.tsv"
        with open(p, "w") as fh:
            fh.writelines(lines)
        print(f"wrote {p}: {len(lines)} pairs")
        band_paths[band] = p
    return pooled, band_paths


def split_strata_by_band(strata_path, bands, outdir):
    """Filter strata file by band using query/reference keys."""
    outdir = Path(outdir)
    band_paths = {}
    band_handles = {}
    with open(strata_path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            q, r = f[0], f[1]
            band = bands.get((q, r))
            if band is None:
                continue
            if band not in band_handles:
                p = outdir / f"strata_{band.replace('.', '_')}.tsv"
                band_paths[band] = p
                band_handles[band] = open(p, "w")
                band_handles[band].write(header)
            band_handles[band].write(line)
    for h in band_handles.values():
        h.close()
    for band, p in sorted(band_paths.items()):
        print(f"wrote {p}")
    return band_paths


def run_panel(syn2bani, strata_path, truth_path, band, outdir):
    print(f"\n=== panel {band} ===")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_file = outdir / f"panel_{band.replace('.', '_')}.txt"
    cmd = [
        syn2bani, "panel",
        "--strata", str(strata_path),
        "--truth", str(truth_path),
        "--greedy",
    ]
    with open(out_file, "w") as fh:
        fh.write(f"# {' '.join(cmd)}\n")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        fh.write(proc.stdout)
    print(f"wrote {out_file} (exit {proc.returncode})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True, help="sample_anim_truth.tsv")
    ap.add_argument("--anim-dir", default="anim_results", help="directory with anim_*.tsv")
    ap.add_argument("--strata", required=True, help="strata.tsv from syn2bani ani --strata-out")
    ap.add_argument("--syn2bani", default="/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani")
    ap.add_argument("--outdir", default="results/panel_by_band")
    args = ap.parse_args()

    if not Path(args.strata).exists():
        print(f"error: strata file not found: {args.strata}", file=sys.stderr)
        return 1

    bands = read_sample(args.sample)
    truth = read_anim_truth(args.anim_dir)
    print(f"collected truth for {len(truth)} pairs")

    pooled_truth, truth_by_band = split_truth_by_band(truth, bands, args.outdir)
    strata_by_band = split_strata_by_band(args.strata, bands, args.outdir)

    # Run pooled panel
    run_panel(args.syn2bani, args.strata, pooled_truth, "all", args.outdir)

    # Run per-band panels
    for band in sorted(set(bands.values())):
        st = strata_by_band.get(band)
        tr = truth_by_band.get(band)
        if st is None or tr is None:
            print(f"skipping {band}: missing strata or truth")
            continue
        run_panel(args.syn2bani, st, tr, band, args.outdir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
