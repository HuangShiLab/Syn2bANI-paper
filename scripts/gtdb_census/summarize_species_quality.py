#!/usr/bin/env python3
"""Summarize natural within-species genome-quality heterogeneity in GTDB."""
import argparse, csv, math, statistics
from collections import Counter, defaultdict
from pathlib import Path


def quantile(xs, q):
    xs=sorted(xs)
    if not xs: return ""
    i=(len(xs)-1)*q; lo=math.floor(i); hi=math.ceil(i)
    return xs[lo] if lo==hi else xs[lo]*(hi-i)+xs[hi]*(i-lo)


def group(comp, cont):
    if comp is None: return "missing"
    if comp < 70: return "low_completeness"
    if cont is not None and cont > 5: return "high_contamination"
    if comp < 90: return "medium"
    return "high"


def num(x):
    try: return float(x)
    except (TypeError,ValueError): return None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--quality",default="/lustre1/g/aos_shihuang/data/syn2b_census/gtdb_genome_quality.tsv")
    ap.add_argument("--out",default="/lustre1/g/aos_shihuang/data/syn2b_census/gtdb_species_quality_summary.tsv")
    a=ap.parse_args()
    d=defaultdict(list)
    with open(a.quality,newline="") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            comp=num(r["checkm_completeness"]); cont=num(r["checkm_contamination"])
            n50=num(r["n50_contigs"]); length=num(r["longest_contig"])
            d[r["species_cluster"]].append((comp,cont,n50,length,group(comp,cont),r["gtdb_taxonomy"]))
    with open(a.out,"w",newline="") as out:
        w=csv.writer(out,delimiter="\t",lineterminator="\n")
        w.writerow(["species_cluster","taxonomy","genomes","high","medium","high_contamination","low_completeness","missing_quality",
                    "frac_high","n_quality_groups","quality_entropy","comp_min","comp_q1","comp_median","comp_q3","comp_max","comp_iqr",
                    "contam_median","contam_max","n50_median","n50_q1","n50_q3"])
        for cl,rows in sorted(d.items()):
            if len(rows)<2: continue
            comps=[x[0] for x in rows if x[0] is not None]
            conts=[x[1] for x in rows if x[1] is not None]
            n50s=[x[2] for x in rows if x[2] is not None]
            lens=[x[3] for x in rows if x[3] is not None]
            groups=Counter(x[4] for x in rows)
            probs=[v/len(rows) for v in groups.values() if v]
            entropy=-sum(p*math.log2(p) for p in probs)
            tax=next((x[5] for x in rows if x[5]),"")
            q1,med,q3=quantile(comps,.25),quantile(comps,.5),quantile(comps,.75)
            w.writerow([cl,tax,len(rows),groups["high"],groups["medium"],groups["high_contamination"],groups["low_completeness"],groups["missing"],
                        groups["high"]/len(rows),len(groups),entropy,min(comps,default=""),q1,med,q3,max(comps,default=""),
                        (q3-q1) if comps else "",quantile(conts,.5),max(conts,default=""),quantile(n50s,.5),
                        quantile(n50s,.25),quantile(n50s,.75)])
    # console overview
    multi=[rows for rows in d.values() if len(rows)>=2]
    allgroups=Counter(x[4] for rows in multi for x in rows)
    mixed=sum(len(set(x[4] for x in rows))>1 for rows in multi)
    print({"multi_species":len(multi),"genomes":sum(map(len,multi)),"species_mixed_quality_group":mixed,
           "frac_mixed":round(mixed/len(multi),4),"genome_quality_groups":dict(allgroups)},flush=True)

if __name__=="__main__": main()
