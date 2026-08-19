#!/usr/bin/env python3
"""build_cohort.py <work>
Merge bin_stats + CheckM2 + per-sample assignment tables into:
  collect/bins.tsv        one row per bin (cohort table)
  pairs/pairs.tsv         bin, ref_id, ref_path, role   (anchor + top-3 same-species alternates)
  pairs/pairs_anchor.tsv  bin, anchor_path              (dnadiff truth list)
Alternates: same-species-cluster source genomes (excluding anchor), ranked by
assigned bp desc; ties/shortfall filled by skani ANI to the majority genome.
"""
import glob
import os
import sys
from collections import defaultdict

WORK = sys.argv[1]

def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(header, line.rstrip("\n").split("\t")))

# --- bin stats ---
bins = {}
for r in read_tsv(f"{WORK}/bins_all/bin_stats.tsv"):
    bins[r["bin"]] = r

# --- checkm2 ---
cm2 = {}
qrep = f"{WORK}/checkm2/quality_report.tsv"
if os.path.exists(qrep):
    for r in read_tsv(qrep):
        cm2[r["Name"]] = r

# --- assignment tables ---
assign = {}
for p in glob.glob(f"{WORK}/assign/*/assign.tsv"):
    for r in read_tsv(p):
        assign[r["bin"]] = r

# --- species clusters + genome paths per dataset ---
cluster_of, gid_path, gid_ani = {}, {}, defaultdict(dict)  # gid_ani[ds][(a,b)] = ani
for ds in ("strain", "marine"):
    gpath = {}
    for r in read_tsv(f"{WORK}/refs/{ds}/genomes.tsv"):
        gpath[r["genome_id"]] = r["path"]
    gid_path[ds] = gpath
    for r in read_tsv(f"{WORK}/refs/{ds}/species_clusters.tsv"):
        cluster_of[r["genome_id"]] = r["cluster_id"]
    # skani triangle for alternate ranking
    tri = f"{WORK}/refs/{ds}/triangle.tsv"
    if os.path.exists(tri):
        with open(tri) as f:
            n = int(f.readline().strip())
            names = []
            for i in range(n):
                parts = f.readline().rstrip("\n").split("\t")
                names.append(parts[0])
                for j, v in enumerate(parts[1:]):
                    try:
                        ani = float(v)
                    except ValueError:
                        continue
                    if j < i and ani > 0:
                        a, b = names[j], parts[0]
                        gid_ani[ds][(a, b)] = gid_ani[ds][(b, a)] = ani
    # (triangle row/col names are full paths; gid_ani keys stay paths)

def tier(comp, cont):
    try:
        c, x = float(comp), float(cont)
    except (TypeError, ValueError):
        return "NA"
    if c >= 90 and x <= 5:
        return "HQ"
    if c >= 50 and x <= 10:
        return "MQ"
    return "LQ"

os.makedirs(f"{WORK}/collect", exist_ok=True)
os.makedirs(f"{WORK}/pairs", exist_ok=True)

n_pairs = 0
with open(f"{WORK}/collect/bins.tsv", "w") as bo, \
     open(f"{WORK}/pairs/pairs.tsv", "w") as po, \
     open(f"{WORK}/pairs/pairs_anchor.tsv", "w") as ao:
    bo.write("bin\tdataset\tsample\tsize_bp\tn50\tn_contigs\tcheckm2_comp\tcheckm2_cont\ttier\t"
             "maj_comp_est\tassigned_frac\tcontam_pct\tclass\tmajority_genome\tn_contigs_multihit\n")
    po.write("bin\tref_id\tref_path\trole\n")
    ao.write("bin\tanchor_path\n")
    for binid, st in sorted(bins.items()):
        ds = st["dataset"]
        a = assign.get(binid, {})
        c = cm2.get(binid, {})
        comp, cont = c.get("Completeness", "NA"), c.get("Contamination", "NA")
        t = tier(comp, cont)
        maj = a.get("majority_genome", "NA")
        total = float(st["size_bp"])
        assigned_frac = (float(a.get("assigned_bp", 0)) / total) if total else 0.0
        bo.write("\t".join(map(str, [binid, ds, st["sample"], st["size_bp"], st["n50"],
                                     st["n_contigs"], comp, cont, t,
                                     a.get("maj_comp_est", "0"), f"{assigned_frac:.4f}",
                                     a.get("contam_pct", "0"), a.get("class", "unassigned"),
                                     maj, a.get("n_contigs_multihit", "0")])) + "\n")
        if maj == "NA" or maj not in gid_path.get(ds, {}):
            continue
        anchor_path = gid_path[ds][maj]
        ao.write(f"{binid}\t{anchor_path}\n")
        po.write(f"{binid}\t{maj}\t{anchor_path}\tanchor\n")
        n_pairs += 1
        # alternates: same cluster, ranked by assigned bp, then skani ANI to majority
        members = [g for g, cl in cluster_of.items()
                   if cl == cluster_of.get(maj) and g != maj and g in gid_path.get(ds, {})]
        ani = gid_ani[ds]  # keys are full paths (skani triangle row/col names)
        def ani_to_maj(g):
            return ani.get((gid_path[ds][g], gid_path[ds][maj]), 0.0)
        members.sort(key=lambda g: -ani_to_maj(g))
        for k, g in enumerate(members[:3], 1):
            po.write(f"{binid}\t{g}\t{gid_path[ds][g]}\talt{k}\n")

print(f"[cohort] {len(bins)} bins, {n_pairs} with anchor")
