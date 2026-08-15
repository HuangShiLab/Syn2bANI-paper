#!/usr/bin/env python3
"""make_hi95_batches.py

Build skani screening batches for the hi95 (95-99.5% ANI) ANIm-truth expansion.

genomes_all contains one representative per GTDB species, so 95-99.5% ANI
pairs can only come from cross-species pairs within the same genus. We pack
all multi-species genera into ~500-genome batches for efficient skani
all-vs-all screening, and dump the taxonomy lookup + existing-pair exclusion
set needed downstream.

Outputs (under WORK):
  batches/batch_NNN.list   genome fasta paths per skani batch
  acc_tax.tsv              accession, phylum, genus, species
  existing_pairs.txt       existing truth pairs, both orientations
"""
import csv
import os
from collections import defaultdict

GDIR = "/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all"
META = "/lustre1/g/aos_shihuang/data/gtdb-r207/metadata"
WORK = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work"
EXISTING = "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/sample_anim_truth.tsv"
BATCH_SIZE = 500  # genomes per skani batch

have = {f[:-4] for f in os.listdir(GDIR) if f.endswith(".fna")}

tax = {}                    # acc -> (phylum, genus, species)
genera = defaultdict(list)  # (phylum, genus) -> [acc]
for mf in ("bac120_metadata_r207.tsv", "ar53_metadata_r207.tsv"):
    with open(f"{META}/{mf}") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = row["accession"].replace("GB_", "").replace("RS_", "")
            if acc not in have:
                continue
            t = row["gtdb_taxonomy"].split(";")
            tax[acc] = (t[1], t[5], t[6])
            genera[(t[1], t[5])].append(acc)

multi = {g: sorted(a) for g, a in genera.items()
         if len({tax[x][2] for x in a}) >= 2}
n_gen = sum(len(v) for v in multi.values())
n_pairs = sum(len(v) * (len(v) - 1) // 2 for v in multi.values())
print(f"multi-species genera: {len(multi)}, genomes: {n_gen}, "
      f"within-genus pairs: {n_pairs}")

# Greedy bin-pack whole genera into batches (large genera first).
items = sorted(multi.items(), key=lambda kv: -len(kv[1]))
batches, sizes = [], []
for g, accs in items:
    for i, s in enumerate(sizes):
        if s + len(accs) <= BATCH_SIZE:
            batches[i].append(g)
            sizes[i] += len(accs)
            break
    else:
        batches.append([g])
        sizes.append(len(accs))

os.makedirs(f"{WORK}/batches", exist_ok=True)
for i, gs in enumerate(batches):
    with open(f"{WORK}/batches/batch_{i:03d}.list", "w") as bl:
        for g in gs:
            for a in multi[g]:
                bl.write(f"{GDIR}/{a}.fna\n")
print(f"batches: {len(batches)}")

with open(f"{WORK}/acc_tax.tsv", "w") as out:
    for a, (p, g, s) in sorted(tax.items()):
        out.write(f"{a}\t{p}\t{g}\t{s}\n")

excl = set()
with open(EXISTING) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    qi, ri = hdr.index("query"), hdr.index("reference")
    for line in f:
        c = line.rstrip("\n").split("\t")
        excl.add((c[qi], c[ri]))
        excl.add((c[ri], c[qi]))
with open(f"{WORK}/existing_pairs.txt", "w") as out:
    for q, r in sorted(excl):
        out.write(f"{q}\t{r}\n")
print(f"existing pairs excluded (both orientations): {len(excl)}")
