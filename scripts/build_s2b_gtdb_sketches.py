#!/usr/bin/env python3
"""Build Syn2bANI sketch database for all GTDB-R207 representatives in batches.

Usage:
    python3 build_s2b_gtdb_sketches.py <ref_list> <out_dir> <batch_size> <threads>
"""
import sys
import os
import subprocess
import math
import time

S2B = "/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani"
ENZYMES = "BcgI,AlfI,AloI,FalI"


def main():
    ref_list = sys.argv[1]
    out_dir = sys.argv[2]
    batch_size = int(sys.argv[3])
    threads = int(sys.argv[4])

    os.makedirs(out_dir, exist_ok=True)

    with open(ref_list) as fh:
        refs = [line.strip() for line in fh if line.strip()]

    total = len(refs)
    n_batches = math.ceil(total / batch_size)

    t0 = time.time()
    for i in range(n_batches):
        batch = refs[i * batch_size:(i + 1) * batch_size]
        bt0 = time.time()
        subprocess.run([S2B, "sketch", "--enzymes", ENZYMES, "-t", str(threads),
                        "-o", out_dir] + batch,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        bt1 = time.time()
        print(f"batch {i+1}/{n_batches}: {len(batch)} genomes in {bt1-bt0:.1f}s")

    t1 = time.time()
    print(f"Total: {total} genomes in {t1-t0:.1f}s ({(t1-t0)/total:.2f}s per genome)")


if __name__ == "__main__":
    main()
