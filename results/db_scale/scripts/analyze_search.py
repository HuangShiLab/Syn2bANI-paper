#!/usr/bin/env python3
"""Search-quality comparison: syn2bani search vs skani search outputs.

Per query: top-hit identity + ANI; on the intersection of reported (query, ref)
pairs: Pearson/Spearman correlation and mean |delta ANI|.
Usage: analyze_search.py <s2b_search.tsv> <skani_search.tsv>
"""
import csv, sys, math, statistics

def load_s2b(path):
    # query_file ref_file ani af_q af_r query_name ref_name shared_tags sv_count
    d = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = row["query_name"]; r = row["ref_name"]
            try: ani = float(row["ani"])
            except ValueError: continue
            d[(q, r)] = ani
    return d

def load_skani(path):
    # Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
    d = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = row["Query_file"].split("/")[-1].replace(".fna", "")
            r = row["Ref_file"].split("/")[-1].replace(".fna", "").replace(".sketch", "")
            try: ani = float(row["ANI"]) / 100.0
            except ValueError: continue
            d[(q, r)] = ani
    return d

s2b = load_s2b(sys.argv[1])
skani = load_skani(sys.argv[2])
print(f"syn2bani reported pairs: {len(s2b)}  ({len(set(q for q,_ in s2b))} queries)")
print(f"skani reported pairs:    {len(skani)}  ({len(set(q for q,_ in skani))} queries)")

def top_hits(d):
    top = {}
    for (q, r), a in d.items():
        if q not in top or a > top[q][1]:
            top[q] = (r, a)
    return top

ts, tk = top_hits(s2b), top_hits(skani)
common_q = sorted(set(ts) & set(tk))
agree_id = sum(1 for q in common_q if ts[q][0] == tk[q][0])
print(f"queries with top hit in both: {len(common_q)}; same top-hit genome: {agree_id}")

inter = sorted(set(s2b) & set(skani))
print(f"intersection pairs: {len(inter)}")
if len(inter) >= 3:
    xs = [s2b[k] for k in inter]; ys = [skani[k] for k in inter]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    cov = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x-mx)**2 for x in xs)); sy = math.sqrt(sum((y-my)**2 for y in ys))
    r = cov/(sx*sy) if sx*sy > 0 else float("nan")
    md = statistics.mean(abs(x-y) for x, y in zip(xs, ys))
    print(f"ANI correlation (Pearson r): {r:.4f}")
    print(f"mean |dANI|: {md:.4f} (fraction units; x100 for ANI points)")
    big = [(k, s2b[k], skani[k]) for k in inter if abs(s2b[k]-skani[k]) > 0.05]
    print(f"pairs with |dANI|>5 points: {len(big)}")
    for k, a, b in big[:10]:
        print(f"  {k[0]} vs {k[1]}: s2b={a:.4f} skani={b:.4f}")
