#!/usr/bin/env python3
"""select_hi95_pairs.py

Select dnadiff candidates from the skani batch screening output.

Keeps unordered same-genus / different-species pairs with skani ANI in
[94.5, 99.8) and min aligned fraction >= MIN_AF, excluding pairs already in
the existing ANIm-truth set (either orientation). Caps pairs per GTDB
species and flattens the ANI distribution by round-robin across 0.5%-wide
skani-ANI bins.

Usage: select_hi95_pairs.py [--n-target 650]
"""
import argparse
import csv
import os
from collections import defaultdict

WORK = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work"
SKANI_GLOB_DIR = f"{WORK}/skani_out"
MIN_ANI, MAX_ANI = 94.5, 99.8
MIN_AF = 15.0
MAX_PER_SPECIES = 5

ap = argparse.ArgumentParser()
ap.add_argument("--n-target", type=int, default=650)
args = ap.parse_args()

tax = {}
with open(f"{WORK}/acc_tax.tsv") as f:
    for line in f:
        a, p, g, s = line.rstrip("\n").split("\t")
        tax[a] = (p, g, s)

excl = set()
with open(f"{WORK}/existing_pairs.txt") as f:
    for line in f:
        q, r = line.rstrip("\n").split("\t")
        excl.add((q, r))

def acc_of(path):
    return os.path.basename(path)[:-4]  # strip .fna

pairs = {}
n_lines = 0
for fn in sorted(os.listdir(SKANI_GLOB_DIR)):
    if not fn.endswith(".tsv"):
        continue
    with open(f"{SKANI_GLOB_DIR}/{fn}") as f:
        hdr = f.readline()
        for line in f:
            n_lines += 1
            c = line.rstrip("\n").split("\t")
            # skani 0.3.x: Ref_file Query_file ANI AF_ref AF_query ...
            r, q = acc_of(c[0]), acc_of(c[1])
            if r == q:
                continue
            ani, af_r, af_q = float(c[2]), float(c[3]), float(c[4])
            if not (MIN_ANI <= ani < MAX_ANI):
                continue
            if min(af_r, af_q) < MIN_AF:
                continue
            if r not in tax or q not in tax:
                continue
            key = tuple(sorted((q, r)))
            if key in pairs:
                continue
            pairs[key] = (q, r, ani, af_r, af_q)
print(f"skani lines: {n_lines}, unique pairs in ANI/AF window: {len(pairs)}")

cand = []
n_dup_excl = n_same_sp = 0
for (q, r), v in pairs.items():
    if (q, r) in excl:
        n_dup_excl += 1
        continue
    pq, gq, sq = tax[q]
    pr, gr, sr = tax[r]
    if gq != gr or pq != pr:
        continue  # cross-genus leakage from batching; keep same-genus only
    if sq == sr:
        n_same_sp += 1
        continue
    cand.append((q, r, v[2], v[3], v[4], pq, gq, sq, sr))
print(f"excluded as existing-truth dup: {n_dup_excl}, same-species: {n_same_sp}")
print(f"candidates after genus/species filter: {len(cand)}")

# Round-robin across 0.5%-wide ANI bins to flatten the distribution.
bins = defaultdict(list)
for row in cand:
    bins[int((row[2] - MIN_ANI) / 0.5)].append(row)
for b in bins.values():
    b.sort()  # deterministic

picked = []
per_species = defaultdict(int)
bin_ids = sorted(bins)
idx = {b: 0 for b in bin_ids}
while len(picked) < args.n_target:
    progress = False
    for b in bin_ids:
        while idx[b] < len(bins[b]):
            row = bins[b][idx[b]]
            idx[b] += 1
            if per_species[row[7]] >= MAX_PER_SPECIES or \
               per_species[row[8]] >= MAX_PER_SPECIES:
                continue
            picked.append(row)
            per_species[row[7]] += 1
            per_species[row[8]] += 1
            progress = True
            break
        if len(picked) >= args.n_target:
            break
    if not progress:
        break

with open(f"{WORK}/pairs_dnadiff.tsv", "w") as out:
    out.write("query\treference\tskani_ani\tskani_af_ref\tskani_af_query\t"
              "phylum\tgenus\tq_species\tr_species\n")
    for row in picked:
        out.write("\t".join(str(x) for x in row) + "\n")

n_sp = len({s for row in picked for s in (row[7], row[8])})
n_ph = len({row[5] for row in picked})
print(f"selected for dnadiff: {len(picked)} "
      f"(distinct species: {n_sp}, phyla: {n_ph})")
from collections import Counter
hist = Counter(int((row[2] - MIN_ANI) / 0.5) for row in picked)
for b in sorted(hist):
    lo = MIN_ANI + b * 0.5
    print(f"  skani ANI {lo:.1f}-{lo + 0.5:.1f}: {hist[b]}")
