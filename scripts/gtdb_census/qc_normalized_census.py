#!/usr/bin/env python3
"""Exact low-memory QC for normalized deduplicated census pair table."""
import argparse, csv, gzip, random, statistics
from collections import defaultdict
from pathlib import Path


def median(xs):
    if not xs: return ""
    xs=sorted(xs); n=len(xs); m=n//2
    return xs[m] if n%2 else (xs[m-1]+xs[m])/2


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--table",required=True)
    ap.add_argument("--workdir",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    root=Path(a.workdir)
    acc_cluster={}; expected={}
    for p in (root/"effective_cluster_accessions").glob("*.txt"):
        accs=p.read_text().split(); cl=p.name[:-4]
        expected[cl]=len(accs)*(len(accs)-1)//2
        for acc in accs: acc_cluster[acc]=cl
    stats=defaultdict(lambda: {"rows":0,"ge2":0,"bp":[]})
    total=0
    with gzip.open(a.table,"rt") as fh:
        header=next(fh).rstrip("\n").split("\t")
        ia=header.index("genome_A"); ib=header.index("genome_B"); ip=header.index("breakpoints")
        for line in fh:
            f=line.rstrip("\n").split("\t"); cl=acc_cluster.get(f[ia])
            if cl is None or acc_cluster.get(f[ib])!=cl: raise SystemExit(f"cross-cluster pair: {line}")
            total+=1; s=stats[cl]; s["rows"]+=1
            try: bp=float(f[ip])
            except ValueError: continue
            s["ge2"] += int(bp>=2); s["bp"].append(bp)
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    with open(out/"effective_species_summary.tsv","w") as fh:
        fh.write("cluster\teffective_genomes\teffective_pairs\trows\tpct_ge2\tmedian_breakpoints\n")
        for cl in sorted(expected):
            s=stats.get(cl)
            if not s: fh.write(f"{cl}\t0\t{expected[cl]}\t0\t0\t\n"); continue
            fh.write(f"{cl}\t{round((1+(1+8*s['rows'])**0.5)/2)}\t{expected[cl]}\t{s['rows']}\t{100*s['ge2']/s['rows']}\t{median(s['bp'])}\n")
    bad=[(c,expected[c],stats[c]["rows"] if c in stats else 0) for c in expected if stats.get(c,{}).get("rows",0)!=expected[c]]
    ge2=sum(s["ge2"] for s in stats.values())
    with open(out/"effective_census_qc.md","w") as fh:
        fh.write("# Effective GTDB structural QC\n\n")
        fh.write(f"- normalized deduplicated pairs: {total:,}\n")
        fh.write(f"- expected effective pairs: {sum(expected.values()):,}\n")
        fh.write(f"- pairs with breakpoints >= 2: {ge2:,}\n")
        fh.write(f"- fraction breakpoints >= 2: {ge2/total:.6%}\n")
        fh.write(f"- clusters with row-count mismatch: {len(bad):,}\n")
        for c,e,r in sorted(bad,key=lambda x:abs(x[1]-x[2]),reverse=True)[:50]: fh.write(f"  - {c}: expected {e:,}, got {r:,}\n")
    print(f"pairs={total} ge2={ge2} frac={ge2/total:.6%} mismatch_clusters={len(bad)}")

if __name__=="__main__": main()
