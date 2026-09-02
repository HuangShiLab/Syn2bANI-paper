#!/usr/bin/env python3
"""Compute Syn2b inverted_fraction for a list of genome pairs.

Workflow:
  1. Extract landmarks from every unique genome, either by restriction digest
     (--mode 2brad, default, with --enzymes) or by FracMinHash (--mode fracminhash,
     with --kmer and --scale).
  2. For each pair, run syn2b synteny on the two TGTs and parse inverted_fraction.

The two landmark sources are not interchangeable inside one run: syn2b refuses a
comparison across them. Give each mode its own --tgt-dir, since the cache keys on
accession alone. FracMinHash k=31 --scale 750 gives ~6,030 landmarks on a 4.54 Mb
genome, matching the four-enzyme panel's ~6,080, which is the density-matched
comparison point.

Inputs:
    --pairs      TSV with columns pairid, q_acc, r_acc
    --genome-dir directory containing {acc}.fna files
    --syn2b      path to syn2b binary
    --out        output TSV

Outputs:
    One TSV row per pair carrying Syn2b's whole structural channel, not only
    inverted_fraction -- the script name is a leftover from what it was first
    written for. Columns:

      syn2b_*                breakpoints, scj_distance, breakpoint_density, both
                             inverted_fraction variants, the three orientation
                             counters, observable_fraction/adjacencies, structural,
                             shared_tags, repeats_dropped, landmarks_collapsed,
                             circular, legacy_adjacency
      syn2b_junctions        breakpoint coordinates in the reference frame,
                             comma-separated
      syn2b_rev_*            the same pair with the roles swapped (--reverse,
                             on by default), including query-frame junctions
      syn2b_scj_corrected    scj_distance - hidden_a - hidden_b

`scj_distance` is the one metric here that is NOT fragmentation-immune -- it is the
raw symmetric difference of the adjacency sets, so each contig break genuinely
removes an adjacency and adds to it. `syn2b_scj_corrected` subtracts that term using
Syn2b's own output, and it needs both directions because `observable_fraction` is
defined on genome_A's adjacencies alone. Measured on E. coli K-12 shattered
independently on both sides with truth SCJ 0, the one-sided correction leaves
+4 / +9 / +18 / +39 / +77 / +141 at K = 5 / 10 / 20 / 40 / 80 / 160 -- no correction
at all -- while the two-sided one leaves +0.2 / +0.1 / +0.1 / +0.2 / +0.1 / -4.0.
With a real 200 kb inversion underneath (truth SCJ 4) it reads 4.2 / 4.1 / 3.8 at
K = 5 / 20 / 80.
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


DEFAULT_ENZYMES = "BcgI,AlfI,AloI,FalI"


def run(cmd, **kw):
    env = os.environ.copy()
    env.setdefault("RAYON_NUM_THREADS", "1")
    kw.setdefault("env", env)
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\nstderr: {r.stderr}\nstdout: {r.stdout}")
    return r


def digest_genome(args):
    acc, genome_path, tgt_dir, syn2b, enzymes, mode, kmer, scale = args
    out_tgt = tgt_dir / f"{acc}.tgt"
    if out_tgt.exists() and out_tgt.stat().st_size > 0:
        return acc, str(out_tgt)
    if mode == "fracminhash":
        cmd = [syn2b, "digest", "-i", str(genome_path), "-o", str(out_tgt),
               "--mode", "fracminhash", "--kmer", str(kmer), "--scale", str(scale),
               "-f", "text"]
    else:
        cmd = [syn2b, "digest", "-i", str(genome_path), "-o", str(out_tgt),
               "-e", enzymes, "-f", "text"]
    try:
        run(cmd, timeout=300)
        return acc, str(out_tgt)
    except Exception as e:
        return acc, f"ERROR: {e}"


SYN2B_COLS = ("breakpoints", "scj_distance", "breakpoint_density",
              "inverted_fraction", "raw_inverted_fraction",
              "orientation_mismatches", "orientation_mismatches_raw",
              "orientation_uninformative", "observable_fraction",
              "observable_adjacencies", "structural", "shared_tags",
              "repeats_dropped", "landmarks_collapsed", "circular",
              # The metric the current one replaced. Kept because the paper's
              # argument is the *contrast*: legacy_adjacency correlates r=+0.982
              # with SynTracker's APSS while responding to divergence rather than
              # order, and the new metric is its mirror image. Dropping it leaves
              # that figure unplottable (MATH_REVIEW.md section 6).
              "legacy_adjacency")


def hidden_adjacencies(rec, prefix=""):
    """Adjacencies of genome_A that genome_B is in no position to judge.

    `observable_fraction = observable_adjacencies / |adj_A|`, so
    `|adj_A| - observable = observable * (1/f - 1)`. This is the quantity
    `scj_distance` is inflated by, one contig break at a time.
    """
    try:
        f = float(rec[f"syn2b_{prefix}observable_fraction"])
        obs = float(rec[f"syn2b_{prefix}observable_adjacencies"])
    except (KeyError, TypeError, ValueError):
        return None
    if not (0.0 < f <= 1.0):
        return None
    return obs * (1.0 / f - 1.0)


def _synteny_once(pairid, a_acc, b_acc, tgt_dir, syn2b, subdir, prefix):
    """Run one direction. `a_acc` becomes genome_A, so all coordinates and the
    fixed-reference raw_inverted_fraction are in its frame."""
    a_tgt = tgt_dir / f"{a_acc}.tgt"
    b_tgt = tgt_dir / f"{b_acc}.tgt"
    if not a_tgt.exists() or not b_tgt.exists():
        return None, "missing_tgt"

    tmpdir = tgt_dir / subdir / pairid
    tmpdir.mkdir(parents=True, exist_ok=True)
    # Remove stale symlinks from previous naming conventions so only the two
    # current TGTs are present per pair.
    for existing in tmpdir.glob("*.tgt"):
        existing.unlink()
    # Name files so a_acc sorts first and becomes genome_A.
    (tmpdir / f"a_ref_{a_acc}.tgt").symlink_to(a_tgt.resolve())
    (tmpdir / f"b_qry_{b_acc}.tgt").symlink_to(b_tgt.resolve())

    out_csv = tmpdir / "synteny.csv"
    try:
        run([syn2b, "synteny", "-i", str(tmpdir), "-o", str(out_csv)], timeout=120)
    except Exception as e:
        (tmpdir / "synteny.err").write_text(str(e))
        return None, f"synteny_error: {e}"

    try:
        with open(out_csv) as fh:
            lines = [line for line in fh if not line.startswith("#")]
        if len(lines) < 2:
            return None, "no_data"
        row = next(csv.DictReader(lines))
    except Exception as e:
        return None, f"parse_error: {e}"

    out = {f"syn2b_{prefix}{c}": row.get(c, "NA") for c in SYN2B_COLS}
    # `syn2b synteny` writes the breakpoint coordinates next to the CSV and this
    # runner used to read only the CSV, so they were produced 43k times and never
    # collected. They are what turns "the two methods agree on counts" into "the
    # two methods agree on positions".
    jpath = out_csv.with_suffix(".junctions.tsv")
    pos = []
    if jpath.exists():
        try:
            with open(jpath) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    pos.append(int(r["junction_pos_in_A"]))
        except Exception:
            pos = []
    out[f"syn2b_{prefix}junctions"] = ",".join(map(str, sorted(pos)))
    return out, "ok"


def run_pair(args):
    pairid, q_acc, r_acc, tgt_dir, syn2b, reverse = args
    # r_acc (reference, matching `dnadiff <ref> <qry>`) is genome_A in the forward
    # direction, so Syn2b's junction coordinates and dnadiff's are in the same frame.
    rec, status = _synteny_once(pairid, r_acc, q_acc, tgt_dir, syn2b, "tmp_pairs", "")
    if rec is None:
        return {"pairid": pairid, "status": status}
    rec["pairid"] = pairid
    rec["status"] = status

    if reverse:
        # The same pair with the roles swapped. Needed because `observable_fraction`
        # is defined on genome_A's adjacencies only, so it accounts for exactly one
        # genome's contig breaks -- and `scj_distance` is inflated by both.
        #
        # Measured on E. coli K-12 shattered independently on both sides, truth
        # SCJ 0: subtracting one side leaves +4 / +9 / +18 / +39 / +77 / +141 at
        # K = 5 / 10 / 20 / 40 / 80 / 160, i.e. no correction at all. Subtracting
        # both leaves +0.2 / +0.1 / +0.1 / +0.2 / +0.1 / -4.0. With a real 200 kb
        # inversion underneath (truth SCJ 4) the two-sided value reads 4.2 / 4.1 /
        # 3.8 at K = 5 / 20 / 80.
        #
        # It also puts the query-side breakpoint coordinates in hand, which is the
        # frame dnadiff's discarded dd.qdiff was in.
        rev, rstatus = _synteny_once(pairid, q_acc, r_acc, tgt_dir, syn2b,
                                     "tmp_pairs_rev", "rev_")
        if rev is not None:
            rec.update(rev)
        rec["reverse_status"] = rstatus

        hA = hidden_adjacencies(rec)
        hB = hidden_adjacencies(rec, "rev_")
        try:
            scj = float(rec["syn2b_scj_distance"])
        except (KeyError, TypeError, ValueError):
            scj = None
        if scj is not None and hA is not None and hB is not None:
            rec["syn2b_scj_corrected"] = f"{scj - hA - hB:.2f}"
            rec["syn2b_hidden_a"] = f"{hA:.2f}"
            rec["syn2b_hidden_b"] = f"{hB:.2f}"
    return rec




def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", required=True)
    p.add_argument("--genome-dir", required=True)
    p.add_argument("--syn2b", default="/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b")
    p.add_argument("--tgt-dir", default=None, help="cache directory for .tgt files")
    p.add_argument("--out", required=True)
    p.add_argument("--workers", type=int, default=min(16, cpu_count() or 1))
    p.add_argument("--enzymes", default=DEFAULT_ENZYMES,
                   help="comma-separated enzyme panel passed to `syn2b digest`. "
                        "Use a --tgt-dir per panel: the cache keys on accession "
                        "only, so reusing one across panels silently mixes them.")
    p.add_argument("--mode", default="2brad", choices=["2brad", "fracminhash"],
                   help="landmark source. `2brad` digests with --enzymes; "
                        "`fracminhash` selects k-mers by hash and uses --kmer/--scale. "
                        "Needs its own --tgt-dir, for the same reason as --enzymes.")
    p.add_argument("--kmer", type=int, default=31,
                   help="FracMinHash k-mer length, 1..32. Ignored in 2brad mode.")
    p.add_argument("--reverse", dest="reverse", action="store_true", default=True,
                   help="also run each pair with the roles swapped (default on). "
                        "Doubles the synteny step, which is cheap next to digestion, "
                        "and is what makes scj_distance usable: observable_fraction "
                        "is defined on genome_A only, so it accounts for one "
                        "genome's contig breaks while scj_distance is inflated by "
                        "both. It also yields query-side junction coordinates.")
    p.add_argument("--no-reverse", dest="reverse", action="store_false",
                   help="forward direction only. scj_corrected is then not emitted.")
    p.add_argument("--scale", type=int, default=750,
                   help="FracMinHash compression: expected landmark spacing in bp. "
                        "750 matches the four-enzyme panel's density on a 4.5 Mb "
                        "genome. Ignored in 2brad mode.")
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
        digest_tasks.append((acc, str(gp), tgt_dir, args.syn2b, args.enzymes,
                             args.mode, args.kmer, args.scale))

    if missing:
        print(f"WARNING: {len(missing)} genomes missing, e.g. {missing[:5]}", flush=True)

    what = (f"FracMinHash k={args.kmer} scale={args.scale}"
            if args.mode == "fracminhash" else f"enzymes {args.enzymes}")
    print(f"extracting landmarks from {len(digest_tasks)} genomes "
          f"[{what}] with {args.workers} workers ...", flush=True)
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
        synteny_tasks.append((row["pairid"], row["q_acc"], row["r_acc"], tgt_dir,
                              args.syn2b, args.reverse))

    print(f"running synteny for {len(synteny_tasks)} pairs "
          f"({'both directions' if args.reverse else 'forward only'}) ...", flush=True)
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
    if args.reverse and "syn2b_scj_corrected" in out_df.columns:
        n = out_df["syn2b_scj_corrected"].notna().sum()
        print(f"  scj_corrected written for {n} pairs "
              f"(= scj_distance - hidden_a - hidden_b; see run_pair)")
    if "syn2b_junctions" in out_df.columns:
        n = (out_df["syn2b_junctions"].fillna("") != "").sum()
        print(f"  junction coordinates collected for {n} pairs")


if __name__ == "__main__":
    main()
