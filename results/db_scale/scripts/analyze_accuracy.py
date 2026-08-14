#!/usr/bin/env python3
"""Accuracy cross-check stats: dist vs ani (MLE truth) from dist_vs_ani.tsv."""
import csv, sys, math, statistics

rows = []
with open(sys.argv[1]) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rows.append(row)

def f(x):
    try: return float(x)
    except (ValueError, TypeError): return None

n = len(rows)
scr = sum(1 for r in rows if r["dist_status"] == "screened_out")
print(f"pairs: {n}; dist screened out (min_af=0.1 default): {scr} ({100*scr/n:.1f}%)")

# sanity: ani rerun vs matrix truth
both = [(f(r["truth_ani_mle_v8"]), f(r["ani_mle_rerun"])) for r in rows]
both = [(a, b) for a, b in both if a is not None and b is not None]
if both:
    d = [abs(a-b) for a, b in both]
    print(f"ani rerun vs matrix truth: n={len(both)} mean|d|={statistics.mean(d):.4f} max|d|={max(d):.4f}")

for col, label in [("dist_chainedkmer", "dist default (chained-kmer)"), ("dist_gbrt_v7", "dist --mash-ani (GBRT v7)")]:
    pairs = [(f(r["truth_ani_mle_v8"]), f(r[col])) for r in rows]
    pairs = [(t, d) for t, d in pairs if t is not None and d is not None]
    zero = sum(1 for t, d in pairs if d == 0.0)
    print(f"\n{label}: reported n={len(pairs)}; ANI==0.0000 despite truth>=80: {zero} ({100*zero/max(len(pairs),1):.1f}%)")
    pos = [(t, d) for t, d in pairs if d > 0]
    if len(pos) >= 3:
        xs = [t for t, _ in pos]; ys = [d*100 for _, d in pos]  # dist ani is 0-1 fraction
        mx, my = statistics.mean(xs), statistics.mean(ys)
        cov = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
        sx = math.sqrt(sum((x-mx)**2 for x in xs)); sy = math.sqrt(sum((y-my)**2 for y in ys))
        r = cov/(sx*sy) if sx*sy > 0 else float("nan")
        md = statistics.mean(abs(x-y) for x, y in zip(xs, ys))
        print(f"  on ANI>0 subset n={len(pos)}: Pearson r={r:.4f} mean|dANI|={md:.2f} points")
        bands = {}
        for t, d in zip(xs, ys):
            b = "80-85" if t < 85 else "85-90" if t < 90 else "90-95" if t < 95 else "95-100"
            bands.setdefault(b, []).append(abs(t-d))
        for b in ["80-85", "85-90", "90-95", "95-100"]:
            if b in bands:
                print(f"    band {b}: n={len(bands[b])} mean|dANI|={statistics.mean(bands[b]):.2f}")
