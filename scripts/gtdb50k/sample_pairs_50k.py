#!/usr/bin/env python3
"""sample_pairs_50k.py — band-stratified 50k-pair sample from the GTDB-R207
representative frame, prescreened by the existing same-genus skani screen.

Frame: hi95_work/skani_out/batch_*.tsv (628,617 non-self same-genus pairs
with skani ANI + aligned fractions). Excludes every pair already in the
calibration training sets (2,074-pair benchmark + 467 hi95 pairs, both
orientations) so the new set is a pure held-out test set.

Sampling: per skani-ANI band targets {80-85: 16k, 85-90: 16k, 90-95: 16k,
95+: all available}; min(skani AF) >= 15; round-robin across 0.25-wide ANI
bins (flattens within-band distribution); <= 5 pairs per genome; <= 30% of
a band's target from any single phylum.

Output: $WORK/pairs_50k.tsv  (q_acc, r_acc, skani_ani, skani_af_min, band,
phylum) plus stdout summary.
"""
import glob
import os
import sys

import pandas as pd

HI95 = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work"
WORK = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k"
EXCL_2074 = f"{WORK}/excl_2074_pairs.tsv"   # query_asm, ref_asm (local eval_pairs.tsv cols)
EXCL_HI95 = f"{WORK}/excl_hi95_pairs.tsv"   # query, reference

TARGETS = {"80-85": 16000, "85-90": 16000, "90-95": 16000, "95-100": 10**9}
BIN_W = 0.25
MAX_PER_GENOME = 5
PHYLUM_CAP_FRAC = 0.30
MIN_AF = 15.0

os.makedirs(WORK, exist_ok=True)

# --- load skani screen ---
cols = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"]
frames = []
for f in sorted(glob.glob(f"{HI95}/skani_out/batch_*.tsv")):
    frames.append(pd.read_csv(f, sep="\t", usecols=cols,
                              dtype={"ANI": float,
                                     "Align_fraction_ref": float,
                                     "Align_fraction_query": float}))
sc = pd.concat(frames, ignore_index=True)
sc = sc[sc.Ref_file != sc.Query_file]
sc["a1"] = sc.Ref_file.str.rsplit("/", n=1).str[-1].str.replace(".fna", "", regex=False)
sc["a2"] = sc.Query_file.str.rsplit("/", n=1).str[-1].str.replace(".fna", "", regex=False)
# single orientation per unordered pair
key = sc[["a1", "a2"]].apply(lambda r: "__".join(sorted(r)), axis=1)
sc["pairkey"] = key
sc = sc.drop_duplicates("pairkey")
sc["af_min"] = sc[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
sc = sc[sc.af_min >= MIN_AF]
sc = sc[(sc.ANI >= 80) & (sc.ANI <= 100)]
print("frame after AF/dedupe filters:", len(sc))

# --- taxonomy ---
tax = {}
with open(f"{HI95}/acc_tax.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        tax[p[0]] = (p[1], p[2], p[3])  # phylum, genus, species

# --- exclusions (training pairs, both orientations) ---
excl = set()
for path, c1, c2 in [(EXCL_2074, "query_asm", "ref_asm"),
                     (EXCL_HI95, "query", "reference")]:
    d = pd.read_csv(path, sep="\t")
    for a, b in zip(d[c1], d[c2]):
        excl.add("__".join(sorted([a, b])))
sc = sc[~sc.pairkey.isin(excl)]
print("frame after excluding training pairs:", len(sc))

# --- bands ---
def band(a):
    if a < 85: return "80-85"
    if a < 90: return "85-90"
    if a < 95: return "90-95"
    return "95-100"
sc["band"] = sc.ANI.map(band)
sc["phylum"] = sc.a1.map(lambda a: tax.get(a, ("?",))[0])

# --- stratified sampling ---
rng = None  # deterministic: sort by ANI, round-robin across fine bins
out_rows = []
for bname, target in TARGETS.items():
    sub = sc[sc.band == bname].copy()
    sub["finebin"] = (sub.ANI // BIN_W).astype(int)
    sub = sub.sort_values(["finebin", "ANI", "pairkey"])
    # round-robin across fine bins
    order = (sub.groupby("finebin", sort=True)
                .apply(lambda g: g.reset_index(drop=True), include_groups=False))
    per_genome = {}
    per_phylum = {}
    phycap = int(target * PHYLUM_CAP_FRAC) if target < 10**8 else 10**9
    taken = 0
    max_len = sub.groupby("finebin").size().max()
    by_bin = {i: g for i, g in sub.groupby("finebin", sort=True)}
    idx = {i: 0 for i in by_bin}
    while taken < min(target, len(sub)):
        progressed = False
        for fb in sorted(by_bin):
            g = by_bin[fb]
            i = idx[fb]
            while i < len(g):
                row = g.iloc[i]
                i += 1
                if per_genome.get(row.a1, 0) >= MAX_PER_GENOME or \
                   per_genome.get(row.a2, 0) >= MAX_PER_GENOME:
                    continue
                if per_phylum.get(row.phylum, 0) >= phycap:
                    continue
                per_genome[row.a1] = per_genome.get(row.a1, 0) + 1
                per_genome[row.a2] = per_genome.get(row.a2, 0) + 1
                per_phylum[row.phylum] = per_phylum.get(row.phylum, 0) + 1
                out_rows.append((row.a1, row.a2, row.ANI, row.af_min, bname,
                                 row.phylum.replace("p__", "")))
                taken += 1
                progressed = True
                break
            idx[fb] = i
            if taken >= min(target, len(sub)):
                break
        if not progressed:
            break
    print(f"band {bname}: took {taken} (available {len(sub)})")

res = pd.DataFrame(out_rows, columns=["q_acc", "r_acc", "skani_ani",
                                      "skani_af_min", "band", "phylum"])
res = res.sample(frac=1.0, random_state=42).reset_index(drop=True)  # shuffle slices
res.to_csv(f"{WORK}/pairs_50k.tsv", sep="\t", index=False,
           float_format="%.2f")
print("wrote", len(res), "pairs ->", f"{WORK}/pairs_50k.tsv")
print(res.band.value_counts())
print("unique genomes:", len(set(res.q_acc) | set(res.r_acc)))
print("phyla:", res.phylum.nunique())
