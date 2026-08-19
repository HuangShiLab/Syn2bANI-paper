#!/usr/bin/env python3
"""cluster_species.py triangle.tsv genomes.tsv ani_threshold > species_clusters.tsv
skani triangle output format (verified): line1 = N; then N rows, row i =
name \t ani_to_genome0 \t ... (lower triangle, self included as last value 100-ish).
genomes.tsv: genome_id \t path (order irrelevant; names in triangle are full paths).
Union-find at ANI >= threshold. Output: genome_id \t cluster_id
"""
import sys

tri_path, genomes_path, thr = sys.argv[1], sys.argv[2], float(sys.argv[3])

with open(genomes_path) as f:
    id_of = {path: gid for gid, path in
             (l.rstrip("\n").split("\t") for l in f if not l.startswith("genome_id\t"))}

with open(tri_path) as f:
    n = int(f.readline().strip())
    names = []
    edges = []
    for i in range(n):
        parts = f.readline().rstrip("\n").split("\t")
        names.append(parts[0])
        for j, v in enumerate(parts[1:]):
            try:
                ani = float(v)
            except ValueError:
                continue
            if j < i and ani >= thr:
                edges.append((i, j))

parent = list(range(n))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb

for a, b in edges:
    union(a, b)

clusters = {}
for i in range(n):
    clusters.setdefault(find(i), []).append(i)

sys.stdout.write("genome_id\tcluster_id\n")
for ci, (_, members) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
    for m in members:
        gid = id_of.get(names[m], names[m])
        sys.stdout.write(f"{gid}\tsp{ci:04d}\n")
