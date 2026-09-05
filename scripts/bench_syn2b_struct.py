#!/usr/bin/env python3
"""Benchmark Syn2bANI struct on all-vs-all pairs.

Usage:
    python3 bench_syn2b_struct.py <fasta_list> <n> <rep> <out_tsv>

fasta_list: file with one FASTA path per line
n: number of genomes (used for output)
rep: repetition number
out_tsv: append row to this file
"""
import sys
import os
import time
import subprocess
import tempfile
import shutil
from itertools import product
from multiprocessing import Pool

SYN2B = "/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani"
THREADS = 16


def run_struct_pair(args):
    """Run syn2bani struct for one query-reference pair."""
    qfile, rfile, outdir = args
    qname = os.path.splitext(os.path.basename(qfile))[0]
    rname = os.path.splitext(os.path.basename(rfile))[0]
    outfile = os.path.join(outdir, f"{qname}__{rname}.tsv")
    cmd = [
        SYN2B, "struct",
        "--enzymes", "BcgI,AlfI,AloI,FalI",
        "--rearrangement", "--indel",
        qfile, rfile,
        "-o", outfile,
    ]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return (qfile, rfile, "ok")
    except subprocess.CalledProcessError:
        return (qfile, rfile, "fail")


def main():
    list_file = sys.argv[1]
    n = int(sys.argv[2])
    rep = int(sys.argv[3])
    out_tsv = sys.argv[4]

    with open(list_file) as fh:
        fastas = [line.strip() for line in fh if line.strip()]

    workdir = tempfile.mkdtemp(prefix=f"s2b_struct_n{n}_r{rep}_")
    pairs = [(q, r, workdir) for q, r in product(fastas, fastas)]

    t0 = time.time()
    with Pool(THREADS) as pool:
        results = pool.map(run_struct_pair, pairs)
    t1 = time.time()
    struct_time = t1 - t0

    n_pairs = len(pairs)
    n_ok = sum(1 for _, _, status in results if status == "ok")

    with open(out_tsv, "a") as fh:
        fh.write(f"syn2b_struct\t{n}\t{n_pairs}\t{rep}\t{struct_time:.3f}\t{n_ok}\n")

    shutil.rmtree(workdir, ignore_errors=True)
    print(f"syn2b_struct n={n} rep={rep}: struct={struct_time:.1f}s ok={n_ok}/{n_pairs}")


if __name__ == "__main__":
    main()
