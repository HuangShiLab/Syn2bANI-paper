#!/usr/bin/env python3
"""Benchmark one-to-all ANI search for one query against the full GTDB-R207 database.

Usage:
    python3 bench_gtdb_one_to_all.py <mode> <query_fasta> <out_tsv>

mode: syn2b_fasta, syn2b_search, skani_fasta, skani_sketch, fastani

Expects environment:
  GTDB_REF_LIST  - file with all reference FASTA paths (for syn2b_fasta, skani_fasta, fastani)
  GTDB_SKANI_DB  - skani sketch database directory
  GTDB_S2B_DB    - Syn2bANI sketch database directory
"""
import sys
import os
import time
import subprocess
import tempfile
import shutil

S2B = "/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani"
SKANI = "/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani"
FASTANI = "/group/aos_shihuang/conda/envs/fastani_env/bin/fastani"
THREADS = 16


def run_cmd(cmd, stdout=None, stderr=None):
    t0 = time.time()
    subprocess.run(cmd, stdout=stdout, stderr=stderr, check=True)
    return time.time() - t0


def write_list(path, items):
    with open(path, "w") as fh:
        for item in items:
            fh.write(item + "\n")


def main():
    mode = sys.argv[1]
    query_fasta = sys.argv[2]
    out_tsv = sys.argv[3]

    base = "/lustre1/g/aos_shihuang/Syn2bANI-paper-bench"
    ref_list = os.environ.get("GTDB_REF_LIST", os.path.join(base, "gtdb_r207_references.txt"))
    skani_db = os.environ.get("GTDB_SKANI_DB", os.path.join(base, "gtdb_r207_skani_sketches"))
    s2b_db = os.environ.get("GTDB_S2B_DB", os.path.join(base, "gtdb_r207_s2b_sketches"))

    workdir = tempfile.mkdtemp(prefix=f"gtdb_1toall_{mode}_")
    query_name = os.path.basename(query_fasta).replace(".fna", "")

    wall_s = None
    peak_rss_mb = "NA"

    if mode == "syn2b_fasta":
        cmd = [S2B, "ani", "--ql", query_fasta, "--rl", ref_list,
               "-t", str(THREADS), "-o", os.path.join(workdir, "out.tsv")]
        wall_s = run_cmd(cmd, subprocess.DEVNULL, subprocess.DEVNULL)

    elif mode == "syn2b_search":
        ql = os.path.join(workdir, "query.txt")
        write_list(ql, [query_fasta])
        cmd = [S2B, "search", "--ql", ql, s2b_db,
               "-t", str(THREADS), "-o", os.path.join(workdir, "out.tsv")]
        wall_s = run_cmd(cmd, subprocess.DEVNULL, subprocess.DEVNULL)

    elif mode == "skani_fasta":
        cmd = [SKANI, "dist", "-t", str(THREADS), "--ql", query_fasta, "--rl", ref_list,
               "-o", os.path.join(workdir, "out.tsv")]
        wall_s = run_cmd(cmd, subprocess.DEVNULL, subprocess.DEVNULL)

    elif mode == "skani_sketch":
        ql = os.path.join(workdir, "query.txt")
        write_list(ql, [query_fasta])
        # Build list of reference sketch files from DB directory
        sk_list = os.path.join(workdir, "ref_sketches.txt")
        write_list(sk_list, [os.path.join(skani_db, f) for f in os.listdir(skani_db)
                             if f.endswith(".sketch")])
        cmd = [SKANI, "dist", "-t", str(THREADS), "--ql", ql, "--rl", sk_list,
               "-o", os.path.join(workdir, "out.tsv")]
        wall_s = run_cmd(cmd, subprocess.DEVNULL, subprocess.DEVNULL)

    elif mode == "fastani":
        cmd = [FASTANI, "--ql", query_fasta, "--rl", ref_list,
               "-t", str(THREADS), "-o", os.path.join(workdir, "out.tsv")]
        wall_s = run_cmd(cmd, subprocess.DEVNULL, subprocess.DEVNULL)

    else:
        raise ValueError(f"Unknown mode: {mode}")

    shutil.rmtree(workdir, ignore_errors=True)

    with open(out_tsv, "a") as fh:
        fh.write(f"{mode}\t{query_name}\t{wall_s:.3f}\t{peak_rss_mb}\n")

    print(f"{mode} {query_name}: {wall_s:.1f}s")


if __name__ == "__main__":
    main()
