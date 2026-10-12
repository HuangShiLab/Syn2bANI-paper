#!/usr/bin/env python3
"""Quality-stratified sensitivity for effective GTDB structural pairs."""
import argparse,csv,gzip

def num(x):
    try:return float(x)
    except (TypeError,ValueError):return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--table",required=True)
    ap.add_argument("--quality",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    q={}
    with open(a.quality) as f:
        for r in csv.DictReader(f,delimiter="\t"):
            q[r["accession"]]={
                "comp":num(r.get("checkm_completeness")),
                "cont":num(r.get("checkm_contamination")),
                "n50":num(r.get("n50_contigs")),
            }
    strata=["all","min_completeness_ge90","min_completeness_lt70","min_completeness_70_90",
            "max_contamination_le5","max_contamination_gt5","HQ_comp90_contam5",
            "min_N50_ge10kb","min_N50_ge100kb","HQ_N50_ge10kb"]
    pairs={s:0 for s in strata}; ge2={s:0 for s in strata}
    n=0; missing_quality=0
    with gzip.open(a.table,"rt") as f:
        header=f.readline().rstrip("\n").split("\t")
        ix={x:i for i,x in enumerate(header)}
        for line in f:
            f=line.rstrip("\n").split("\t"); n+=1
            a1=f[ix["genome_A"]]; b1=f[ix["genome_B"]]
            qa=q.get(a1); qb=q.get(b1)
            if qa is None or qb is None: missing_quality+=1; continue
            mc=min(qa["comp"],qb["comp"]); xc=max(qa["cont"],qb["cont"]); mn=min(qa["n50"],qb["n50"])
            try: bp=float(f[ix["breakpoints"]])
            except ValueError: continue
            add=[
                ("all",1),
                ("min_completeness_ge90",int(mc is not None and mc>=90)),
                ("min_completeness_lt70",int(mc is not None and mc<70)),
                ("min_completeness_70_90",int(mc is not None and 70<=mc<90)),
                ("max_contamination_le5",int(xc is not None and xc<=5)),
                ("max_contamination_gt5",int(xc is not None and xc>5)),
                ("HQ_comp90_contam5",int(mc is not None and mc>=90 and xc is not None and xc<=5)),
                ("min_N50_ge10kb",int(mn is not None and mn>=10000)),
                ("min_N50_ge100kb",int(mn is not None and mn>=100000)),
                ("HQ_N50_ge10kb",int(mc is not None and mc>=90 and xc is not None and xc<=5 and mn is not None and mn>=10000)),
            ]
            for s,ok in add:
                pairs[s]+=1; ge2[s]+=ok
    with open(a.output,"w") as o:
        o.write("stratum\tpairs\tpairs_ge2\tfraction_ge2\n")
        for s in strata:
            o.write(f"{s}\t{pairs[s]}\t{ge2[s]}\t{ge2[s]/pairs[s]:.8f}\n")
    print(f"rows={n} missing_quality={missing_quality}",file=__import__('sys').stderr)

if __name__=="__main__": main()
