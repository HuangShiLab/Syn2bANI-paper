#!/usr/bin/env python3
"""stage_bins.py <bins_root> <bins_all_dir> <bin_stats.tsv> <min_bp>
Collect MetaBAT2 bins from bins_root/{dataset}__{sample}/bin.*.fa, keep bins
>= min_bp total, rewrite fasta headers to {binid}|{orig} (syn2bani reports the
first record id, so the binid must lead), and emit per-bin stats.
binid = {dataset}__{sample}__{bin_stem}   (bin_stem e.g. bin.3)
"""
import os
import sys

bins_root, out_dir, stats_out, min_bp = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
os.makedirs(out_dir, exist_ok=True)

def fasta_lengths(path):
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

def n50(lens):
    s = sorted(lens, reverse=True)
    half = sum(s) / 2.0
    acc = 0
    for L in s:
        acc += L
        if acc >= half:
            return L
    return 0

rows = []
for tag in sorted(os.listdir(bins_root)):
    d = os.path.join(bins_root, tag)
    if not os.path.isdir(d):
        continue
    dataset, sample = tag.split("__", 1)
    for fn in sorted(os.listdir(d)):
        if not (fn.startswith("bin.") and fn.endswith(".fa")):
            continue
        src = os.path.join(d, fn)
        lens = fasta_lengths(src)
        size = sum(lens)
        if size < min_bp:
            continue
        binid = f"{dataset}__{sample}__{fn[:-3]}"
        dst = os.path.join(out_dir, binid + ".fa")
        with open(src) as fi, open(dst, "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    fo.write(f">{binid}|{line[1:]}")
                else:
                    fo.write(line)
        rows.append((binid, dataset, sample, size, n50(lens), len(lens), dst))

with open(stats_out, "w") as f:
    f.write("bin\tdataset\tsample\tsize_bp\tn50\tn_contigs\tpath\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print(f"[stage_bins] staged {len(rows)} bins >= {min_bp} bp into {out_dir}")
