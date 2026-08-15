#!/usr/bin/env python3
"""finalize_hi95.py

Parse dnadiff .report files for the hi95 candidate pairs, keep pairs whose
1-to-1 ANIm ANI falls in [95, 99.5], and write the final deliverables:

  results/anim_truth_hi95.tsv        query, reference, anim_ani,
                                     anim_aligned_query, anim_aligned_ref
  results/anim_truth_hi95_pairs.tsv  pair metadata (taxonomy + skani + ANIm)

Also runs QA: duplicate check against the existing 2,074-pair truth set
(both orientations), ANI-band histogram, species/phylum counts.
"""
import csv
import re
import sys
from collections import Counter

WORK = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work"
RESULTS = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results"
EXISTING = f"{RESULTS}/sample_anim_truth.tsv"
LO_ANI, HI_ANI = 95.0, 99.5

def parse_report(path):
    txt = open(path).read()
    def grab(label):
        m = re.search(rf"^{label}\s+(.*)$", txt, re.M)
        return m.group(1).split() if m else []
    ab = grab("AlignedBases")   # ['123(45.67%)', '124(46.78%)'] (1-to-1 block)
    ai = grab("AvgIdentity")    # ['97.12', '97.23']
    aq = float(re.search(r"\(([\d.]+)%\)", ab[0]).group(1))
    ar = float(re.search(r"\(([\d.]+)%\)", ab[1]).group(1))
    ident = (float(ai[0]) + float(ai[1])) / 2.0
    return ident, aq, ar

meta = {}
with open(f"{WORK}/pairs_dnadiff.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        meta[(row["query"], row["reference"])] = row

rows, failed, out_of_band = [], [], 0
for (q, r), m in sorted(meta.items()):
    rep = f"{WORK}/dnadiff/{q}__{r}.report"
    try:
        ident, aq, ar = parse_report(rep)
    except Exception:
        failed.append((q, r))
        continue
    if not (LO_ANI <= ident <= HI_ANI):
        out_of_band += 1
        continue
    rows.append((q, r, ident, aq, ar, m))

print(f"dnadiff pairs attempted: {len(meta)}")
print(f"parse failures: {len(failed)}, ANIm out of [95,99.5]: {out_of_band}")
print(f"kept: {len(rows)}")

with open(f"{RESULTS}/anim_truth_hi95.tsv", "w") as out:
    out.write("query\treference\tanim_ani\tanim_aligned_query\tanim_aligned_ref\n")
    for q, r, ident, aq, ar, _ in rows:
        out.write(f"{q}\t{r}\t{ident:.4f}\t{aq:.2f}\t{ar:.2f}\n")

with open(f"{RESULTS}/anim_truth_hi95_pairs.tsv", "w") as out:
    out.write("query\treference\tphylum\tgenus\tq_species\tr_species\t"
              "skani_ani\tskani_af_ref\tskani_af_query\tanim_ani\n")
    for q, r, ident, aq, ar, m in rows:
        out.write(f"{q}\t{r}\t{m['phylum']}\t{m['genus']}\t{m['q_species']}\t"
                  f"{m['r_species']}\t{m['skani_ani']}\t{m['skani_af_ref']}\t"
                  f"{m['skani_af_query']}\t{ident:.4f}\n")

# --- QA ----------------------------------------------------------------------
excl = set()
with open(EXISTING) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    qi, ri = hdr.index("query"), hdr.index("reference")
    for line in f:
        c = line.rstrip("\n").split("\t")
        excl.add((c[qi], c[ri]))
        excl.add((c[ri], c[qi]))
dups = [(q, r) for q, r, *_ in rows if (q, r) in excl or (r, q) in excl]
print(f"QA duplicates vs existing truth set: {len(dups)}")
if dups:
    print("  DUP:", dups[:5])

hist = Counter(int((ident - LO_ANI) / 0.5) for _, _, ident, *_ in rows)
print("QA ANIm band histogram:")
for b in sorted(hist):
    lo = LO_ANI + b * 0.5
    print(f"  {lo:.1f}-{lo + 0.5:.1f}: {hist[b]}")
n_sp = len({s for *_, m in rows for s in (m["q_species"], m["r_species"])})
phy = Counter(m["phylum"] for *_, m in rows)
print(f"QA distinct species: {n_sp}, phyla: {len(phy)}")
print("QA top phyla:", phy.most_common(10))
if failed:
    print("failed pairs (first 10):", failed[:10])
