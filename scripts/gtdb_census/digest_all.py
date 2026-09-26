#!/usr/bin/env python3
"""One-time pre-digestion of every census genome (thundering-herd fix).

census_worker.py digests missing TGTs per task, which races badly on mega
clusters (every task attempt re-digests the same genomes). This pass digests
each genome exactly once, with a per-accession mkdir claim, so the census
workers later find everything present.

Idempotent and multi-node safe: run as a SLURM array over shards, or several
single-node jobs; already-digested accessions are skipped.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def sanitize_fasta(src, tmpdir, genome_id=None):
    """Copy a FASTA, stripping ENA's leading 'ENA|' header field so genome and
    contig IDs are the unique INSDC accessions.  For catalogs whose FASTA
    headers are contig IDs (HROM), ``genome_id`` prefixes every header so all
    contigs in one file become one TGT genome."""
    dst = Path(tmpdir) / (Path(src).name + ".sanitized.fna")
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue
            header = line[1:].rstrip("\n")
            if genome_id is not None:
                if not header.startswith(genome_id + "|"):
                    header = f"{genome_id}|{header}"
            elif header.startswith("ENA|"):
                header = header[4:]
            fout.write(f">{header}\n")
    return dst


def digest_one(syn2b, enzymes, fasta, tgt_out, tmpdir, genome_id=None):
    fasta = sanitize_fasta(fasta, tmpdir, genome_id)
    try:
        subprocess.run([syn2b, "digest", "-i", str(fasta), "-o", str(tgt_out),
                        "-e", enzymes], check=True, capture_output=True)
    finally:
        fasta.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--syn2b", required=True)
    ap.add_argument("--enzymes", default="BcgI,AlfI,AloI,FalI")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--purge-broken", action="store_true",
                    help="first delete any of this shard's TGTs whose genome "
                         "id is the broken literal 'ENA' token (leftover from "
                         "unsanitized ENA digests)")
    ap.add_argument("--ensure-filename-genome-id", action="store_true",
                    help="prefix each sanitized FASTA header with the manifest "
                         "accession; required for HROM contig headers")
    ap.add_argument("--max-per-core", type=int, default=0,
                    help="stop after N digests (0 = unlimited); use for chunked runs")
    args = ap.parse_args()

    root = Path(args.workdir)
    tgt = root / "tgt"
    locks = root / "locks"
    tgt.mkdir(parents=True, exist_ok=True)
    locks.mkdir(parents=True, exist_ok=True)
    manifest = __import__("json").loads((root / "manifest.json").read_text())

    items = sorted(manifest.items())
    mine = items[args.shard::args.nshards]
    if args.purge_broken:
        purged = 0
        for acc, _ in mine:
            t = tgt / f"{acc}.tgt"
            if t.exists() and t.stat().st_size == 0:
                t.unlink(); purged += 1; continue
            try:
                with t.open() as fh:
                    if fh.readline().startswith(">ENA|"):
                        t.unlink(); purged += 1
            except OSError:
                pass
        if purged:
            print(f"shard {args.shard}: purged {purged} broken TGTs", flush=True)
    done = 0
    t0 = time.time()
    for acc, fasta in mine:
        if args.max_per_core and done >= args.max_per_core:
            break
        if (tgt / f"{acc}.tgt").exists():
            continue
        try:
            (locks / f"digest.{acc}").mkdir()
        except FileExistsError:
            continue
        try:
            if not os.path.exists(fasta):
                print(f"MISSING-FASTA {acc}", flush=True)
                (locks / f"digest.{acc}").rmdir()
                continue
            try:
                genome_id = acc if args.ensure_filename_genome_id else None
                digest_one(args.syn2b, args.enzymes, fasta,
                           tgt / f"{acc}.tgt", tempfile.gettempdir(),
                           genome_id)
                done += 1
                if done % 500 == 0:
                    rate = done / (time.time() - t0)
                    print(f"shard {args.shard}: {done} digested "
                          f"({rate:.0f}/s)", flush=True)
            finally:
                try:
                    (locks / f"digest.{acc}").rmdir()
                except OSError:
                    pass
        except Exception as e:
            print(f"FAIL {acc}: {str(e)[:200]}", flush=True)
            # leave the tgt absent and the lock released so a later pass can retry
            tgt_file = tgt / f"{acc}.tgt"
            if tgt_file.exists() and tgt_file.stat().st_size == 0:
                tgt_file.unlink()
            try:
                (locks / f"digest.{acc}").rmdir()
            except OSError:
                pass
    print(f"shard {args.shard}: DONE {done} new digests", flush=True)


if __name__ == "__main__":
    main()
