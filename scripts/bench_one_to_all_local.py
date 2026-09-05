#!/usr/bin/env python3
"""Benchmark one-to-all ANI search locally (Mac Studio).

Usage:
    python3 bench_one_to_all_local.py <fasta_list> <n> <rep> <out_tsv>

fasta_list: file with one FASTA path per line (first line is the query)
n: number of references to use from the list (excluding query)
rep: repetition number
out_tsv: append row to this file
"""
import sys
import os
import time
import subprocess
import tempfile
import shutil

SYN = "/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani"
SKANI = "/Users/macstudio/.cargo/bin/skani"
FASTANI = "/opt/homebrew/bin/fastani"
THREADS = 16


def run_and_log(tool, mode, n, rep, out_tsv, cmd):
    """Run cmd and capture wall time; append to TSV."""
    t0 = time.time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    wall = time.time() - t0
    with open(out_tsv, "a") as fh:
        fh.write(f"{tool}\t{mode}\t{n}\t{rep}\t{wall:.3f}\tNA\n")
    print(f"{tool} {mode} n={n} rep={rep}: wall={wall:.2f}s")


def main():
    list_file = sys.argv[1]
    n = int(sys.argv[2])
    rep = int(sys.argv[3])
    out_tsv = sys.argv[4]

    with open(list_file) as fh:
        fastas = [line.strip() for line in fh if line.strip()]

    query = fastas[0]
    refs = fastas[1:n+1]

    workdir = tempfile.mkdtemp(prefix=f"one_to_all_n{n}_r{rep}_")

    query_list = os.path.join(workdir, "query_list.txt")
    with open(query_list, "w") as fh:
        fh.write(query + "\n")
    ref_list = os.path.join(workdir, "ref_list.txt")
    with open(ref_list, "w") as fh:
        for f in refs:
            fh.write(f + "\n")

    # syn2bani FASTA mode
    run_and_log("syn2bani", "one_to_all_fasta", n, rep, out_tsv,
                [SYN, "ani", "--ql", query_list, "--rl", ref_list,
                 "-t", "0", "-o", os.path.join(workdir, "s2b_fasta.tsv")])

    # syn2bani sketch-reuse mode
    s2b_db = os.path.join(workdir, "s2b_db")
    os.makedirs(s2b_db, exist_ok=True)
    subprocess.run([SYN, "sketch", "--enzymes", "BcgI,AlfI,AloI,FalI", "-t", "0",
                    "-o", s2b_db] + refs + [query],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    query_s2ba = os.path.join(s2b_db, os.path.basename(query).replace(".fasta", ".s2ba"))
    ref_s2ba_list = os.path.join(workdir, "ref_s2ba_list.txt")
    query_s2ba_list = os.path.join(workdir, "query_s2ba_list.txt")
    with open(query_s2ba_list, "w") as fh:
        fh.write(query_s2ba + "\n")
    with open(ref_s2ba_list, "w") as fh:
        for f in refs:
            fh.write(os.path.join(s2b_db, os.path.basename(f).replace(".fasta", ".s2ba")) + "\n")
    run_and_log("syn2bani", "one_to_all_sketches", n, rep, out_tsv,
                [SYN, "ani", "--ql", query_s2ba_list, "--rl", ref_s2ba_list,
                 "-t", "0", "-o", os.path.join(workdir, "s2b_sk.tsv")])

    # skani sketch + dist
    sk_db = os.path.join(workdir, "sk_db")
    if os.path.exists(sk_db):
        shutil.rmtree(sk_db)
    subprocess.run([SKANI, "sketch", "-t", str(THREADS), "-o", sk_db] + refs + [query],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    query_sk = os.path.join(sk_db, os.path.basename(query) + ".sketch")
    ref_sk_list = os.path.join(workdir, "ref_sk_list.txt")
    query_sk_list = os.path.join(workdir, "query_sk_list.txt")
    with open(query_sk_list, "w") as fh:
        fh.write(query_sk + "\n")
    with open(ref_sk_list, "w") as fh:
        for f in refs:
            fh.write(os.path.join(sk_db, os.path.basename(f) + ".sketch") + "\n")
    run_and_log("skani", "one_to_all_dist", n, rep, out_tsv,
                [SKANI, "dist", "-t", str(THREADS),
                 "--ql", query_sk_list, "--rl", ref_sk_list,
                 "-o", os.path.join(workdir, "sk.tsv")])

    # FastANI
    run_and_log("fastani", "one_to_all", n, rep, out_tsv,
                [FASTANI, "--ql", query_list, "--rl", ref_list,
                 "-t", str(THREADS), "-o", os.path.join(workdir, "fa.tsv")])

    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
