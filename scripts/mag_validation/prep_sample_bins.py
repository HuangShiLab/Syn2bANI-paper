#!/usr/bin/env python3
"""prep_sample_bins.py <bins_dir> <dataset> <sample> <out.fa> <min_bp>
Concatenate one sample's MetaBAT2 bins (>= min_bp) into a single fasta with
headers {binid}|{contig} for minimap2 contig-truth assignment.
"""
import os
import sys

bins_dir, dataset, sample, out_fa, min_bp = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])

def fasta_lens(path):
    lens, cur = [], 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lens.append(cur)
    return lens

kept = 0
with open(out_fa, "w") as out:
    for fn in sorted(os.listdir(bins_dir)):
        if not (fn.startswith("bin.") and fn.endswith(".fa")):
            continue
        src = os.path.join(bins_dir, fn)
        if sum(fasta_lens(src)) < min_bp:
            continue
        binid = f"{dataset}__{sample}__{fn[:-3]}"
        kept += 1
        with open(src) as fi:
            for l in fi:
                if l.startswith(">"):
                    out.write(f">{binid}|{l[1:]}")
                else:
                    out.write(l)
print(f"[prep_sample_bins] {dataset}__{sample}: kept {kept} bins >= {min_bp} bp")
