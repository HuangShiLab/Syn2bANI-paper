#!/usr/bin/env python3
"""assemble_hi95_gated.py

Assemble the per-pair `syn2bani ani --verbose` fragments (commit 15da386,
default 4-enzyme panel BcgI,AlfI,AloI,FalI) into
results/anim_truth_hi95_gated.tsv with the exact column layout of
results/anim_truth_2074_gated.tsv. query/reference are re-keyed to the
assembly accessions (the binary prints FASTA genome ids).

QA: row count, field-count check, flag/gate distributions, NaN census.
"""
import csv
import os
import sys
from collections import Counter

WORK = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work"
RESULTS = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results"
PAIRS = f"{RESULTS}/anim_truth_hi95.tsv"
OUT = f"{RESULTS}/anim_truth_hi95_gated.tsv"

HEADER = ("query reference ani ani_uniform af_query af_reference std_err "
          "synteny_blocks anchor_adjacency breakpoint_count het_shape retention "
          "ani_from_loss ani_from_hist enzyme_spread enzyme_chi2 per_enzyme "
          "n_anchors n_chains n_tags max_block_anchors mean_block_anchors "
          "flag ani_gated gate").split()

pairs = []
with open(PAIRS) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        pairs.append((row["query"], row["reference"]))

rows, missing, bad = [], [], []
for q, ref in pairs:
    frag = f"{WORK}/s2b_out/{q}__{ref}.tsv"
    if not (os.path.exists(frag) and os.path.getsize(frag) > 0):
        missing.append((q, ref))
        continue
    lines = [l for l in open(frag).read().splitlines() if l.strip()]
    if len(lines) != 1:
        bad.append((q, ref, f"{len(lines)} lines"))
        continue
    fields = lines[0].split("\t")
    if len(fields) != len(HEADER):
        bad.append((q, ref, f"{len(fields)} fields"))
        continue
    fields[0], fields[1] = q, ref  # re-key to accessions
    rows.append(fields)

with open(OUT, "w") as out:
    out.write("\t".join(HEADER) + "\n")
    for f_ in rows:
        out.write("\t".join(f_) + "\n")

print(f"pairs: {len(pairs)}, assembled: {len(rows)}, missing: {len(missing)}, "
      f"bad: {len(bad)}")
if missing:
    print("missing:", missing[:10])
if bad:
    print("bad:", bad[:10])

flags = Counter(f_[22] for f_ in rows)
gates = Counter(f_[24] for f_ in rows)
print("flag distribution:", dict(flags))
print("gate distribution:", dict(gates))
nan_gated = sum(1 for f_ in rows if f_[23] in ("NaN", "", "nan"))
print(f"ani_gated NaN (legit below-detection/none): {nan_gated}")
# every NaN ani_gated row must have gate == 'none'
bad_nan = [f_[:2] for f_ in rows
           if f_[23] in ("NaN", "", "nan") and f_[24] != "none"]
print(f"NaN ani_gated with gate != none: {len(bad_nan)}", bad_nan[:5])
