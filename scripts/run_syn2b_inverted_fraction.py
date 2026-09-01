#!/usr/bin/env python3
"""Compute Syn2b inverted_fraction for a list of genome pairs.

Workflow:
  1. Digest every unique genome with the 4-enzyme panel (BcgI,AlfI,AloI,FalI).
  2. For each pair, run syn2b synteny on the two TGTs and parse inverted_fraction.

Inputs:
    --pairs      TSV with columns pairid, q_acc, r_acc
    --genome-dir directory containing {acc}.fna files
    --syn2b      path to syn2b binary
    --out        output TSV

Outputs:
    TSV with pairid, syn2b_breakpoints, syn2b_inverted_fraction, syn2b_observable_fraction, ...
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd


ENZYMES = "BcgI,AlfI,AloI,FalI"


def run(cmd, **kw):
    env = os.environ.copy()
    env.setdefault("RAYON_NUM_THREADS", "1")
    kw.setdefault("env", env)
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\nstderr: {r.stderr}\nstdout: {r.stdout}")
    return r


def digest_genome(args):
    acc, genome_path, tgt_dir, syn2b = args
    out_tgt = tgt_dir / f"{acc}.tgt"
    if out_tgt.exists() and out_tgt.stat().st_size > 0:
        return acc, str(out_tgt)
    try:
        run([syn2b, "digest", "-i", str(genome_path), "-o", str(out_tgt),
             "-e", ENZYMES, "-f", "text"], timeout=120)
        return acc, str(out_tgt)
    except Exception as e:
        return acc, f"ERROR: {e}"


def run_pair(args):
    pairid, q_acc, r_acc, tgt_dir, syn2b = args
    q_tgt = tgt_dir / f"{q_acc}.tgt"
    r_tgt = tgt_dir / f"{r_acc}.tgt"
    if not q_tgt.exists() or not r_tgt.exists():
        return {"pairid": pairid, "status": "missing_tgt"}

    tmpdir = tgt_dir / "tmp_pairs" / pairid
    tmpdir.mkdir(parents=True, exist_ok=True)
    q_link = tmpdir / f"{q_acc}.tgt"
    r_link = tmpdir / f"{r_acc}.tgt"
    if not q_link.is_symlink():
        q_link.symlink_to(q_tgt.resolve())
    if not r_link.is_symlink():
        r_link.symlink_to(r_tgt.resolve())

    out_csv = tmpdir / "synteny.csv"
    try:
        run([syn2b, "synteny", "-i", str(tmpdir), "-o", str(out_csv)], timeout=120)
    except Exception as e:
        err_file = tmpdir / "synteny.err"
        err_file.write_text(str(e))
        return {"pairid": pairid, "status": f"synteny_error: {e}"}

    try:
        with open(out_csv) as fh:
            # skip comment lines starting with '#'
            lines = [line for line in fh if not line.startswith("#")]
            if len(lines) < 2:
                return {"pairid": pairid, "status": "no_data"}
            reader = csv.DictReader(lines)
            row = next(reader)
            return {
                "pairid": pairid,
                "status": "ok",
                "syn2b_breakpoints": row["breakpoints"],
                "syn2b_scj_distance": row["scj_distance"],
                "syn2b_breakpoint_density": row["breakpoint_density"],
                "syn2b_inverted_fraction": row["inverted_fraction"],
                "syn2b_orientation_mismatches": row["orientation_mismatches"],
                "syn2b_orientation_uninformative": row["orientation_uninformative"],
                "syn2b_observable_fraction": row["observable_fraction"],
                "syn2b_observable_adjacencies": row["observable_adjacencies"],
                "syn2b_structural": row["structural"],
                "syn2b_shared_tags": row["shared_tags"],
                "syn2b_repeats_dropped": row["repeats_dropped"],
                "syn2b_landmarks_collapsed": row["landmarks_collapsed"],
                "syn2b_circular": row["circular"],
            }
    except Exception as e:
        return {"pairid": pairid, "status": f"parse_error: {e}"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", required=True)
    p.add_argument("--genome-dir", required=True)
    p.add_argument("--syn2b", default="/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b")
    p.add_argument("--tgt-dir", default=None, help="cache directory for .tgt files")
    p.add_argument("--out", required=True)
    p.add_argument("--workers", type=int, default=min(16, cpu_count() or 1))
    args = p.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    if "pairid" not in pairs.columns:
        pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]

    # genome path lookup
    genome_dir = Path(args.genome_dir)
    accs = set(pairs["q_acc"]) | set(pairs["r_acc"])
    print(f"{len(pairs)} pairs, {len(accs)} unique genomes", flush=True)

    # tgt cache
    if args.tgt_dir:
        tgt_dir = Path(args.tgt_dir)
    else:
        tgt_dir = Path(tempfile.mkdtemp(prefix="syn2b_tgts_"))
    tgt_dir.mkdir(parents=True, exist_ok=True)
    print(f"tgt cache: {tgt_dir}", flush=True)

    # digest
    digest_tasks = []
    missing = []
    for acc in accs:
        gp = genome_dir / f"{acc}.fna"
        if not gp.exists():
            # try .fasta
            gp = genome_dir / f"{acc}.fasta"
        if not gp.exists():
            missing.append(acc)
            continue
        digest_tasks.append((acc, str(gp), tgt_dir, args.syn2b))

    if missing:
        print(f"WARNING: {len(missing)} genomes missing, e.g. {missing[:5]}", flush=True)

    print(f"digesting {len(digest_tasks)} genomes with {args.workers} workers ...", flush=True)
    with Pool(args.workers) as pool:
        digest_results = pool.imap_unordered(digest_genome, digest_tasks, chunksize=10)
        for i, (acc, result) in enumerate(digest_results):
            if i % 500 == 0 and i > 0:
                print(f"  digested {i}/{len(digest_tasks)}", flush=True)
            if result.startswith("ERROR"):
                print(f"  digest failed for {acc}: {result}", flush=True)

    # synteny per pair
    synteny_tasks = []
    for _, row in pairs.iterrows():
        synteny_tasks.append((row["pairid"], row["q_acc"], row["r_acc"], tgt_dir, args.syn2b))

    print(f"running synteny for {len(synteny_tasks)} pairs ...", flush=True)
    records = []
    with Pool(max(1, args.workers // 2)) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_pair, synteny_tasks, chunksize=10)):
            if i % 500 == 0 and i > 0:
                print(f"  processed {i}/{len(synteny_tasks)}", flush=True)
            records.append(rec)

    out_df = pd.DataFrame(records)
    out_df = pairs[["pairid", "q_acc", "r_acc"]].merge(out_df, on="pairid", how="left")
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}: {len(out_df)} rows, {(out_df['status'] == 'ok').sum()} ok")


if __name__ == "__main__":
    main()
