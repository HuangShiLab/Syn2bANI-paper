#!/usr/bin/env python3
"""CheckM-stratified sensitivity analysis for the HROM hypermode result."""
from __future__ import annotations
import argparse, array, csv, gzip, json, math
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom, mannwhitneyu, spearmanr, fisher_exact


def bh_qvalues(p):
    p=np.asarray(p,dtype=float); n=len(p)
    if not n: return np.array([],dtype=float)
    order=np.argsort(p); ranked=p[order]*n/np.arange(1,n+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]
    ranked=np.minimum(ranked,1.0)
    q=np.empty(n); q[order]=ranked
    return q


def clean_species(raw,cluster):
    x=(raw or "").strip()
    if x.startswith("s__"): x=x[3:]
    if "/HROM_Genome_" in x: x=x.split("/HROM_Genome_",1)[0]
    return x or f"HROM {cluster}"


def quality_tier(mincomp,maxcont):
    if mincomp>=90 and maxcont<=5: return "high"
    if mincomp>=70 and maxcont<=10: return "medium"
    return "low"


def species_enrichment(mask,codes,names,N,n_top,min_pairs):
    total=np.bincount(codes[mask],minlength=len(names))
    struct=np.bincount(codes[mask & struct_top],minlength=len(names))
    ani=np.bincount(codes[mask & ani_top],minlength=len(names))
    both=np.bincount(codes[mask & both_top],minlength=len(names))
    rows=[]
    for code,name in enumerate(names):
        K=int(total[code])
        if K<min_pairs: continue
        ks,ka,kb=int(struct[code]),int(ani[code]),int(both[code])
        ps=float(hypergeom.sf(ks-1,N,K,n_top))
        pa=float(hypergeom.sf(ka-1,N,K,n_top))
        # logsf protects extreme enrichments from underflow.
        lps=float(hypergeom.logsf(ks-1,N,K,n_top))
        lpa=float(hypergeom.logsf(ka-1,N,K,n_top))
        l2=math.log2((ks+0.5)/(ka+0.5))
        qs=float(bh_qvalues([ps])[0]); qa=float(bh_qvalues([pa])[0])
        if qs<.05 and l2>=1: mode="structural_enriched"
        elif qa<.05 and l2<=-1: mode="ani_enriched"
        elif qs<.05 and qa<.05: mode="both_enriched"
        else: mode="not_enriched"
        rows.append({"cluster":name,"species":labels[name],"total_pairs":K,
                     "struct_top_pairs":ks,"ani_top_pairs":ka,"both_top_pairs":kb,
                     "struct_enrichment_ratio":(ks/n_top)/(K/N),
                     "ani_enrichment_ratio":(ka/n_top)/(K/N),
                     "struct_p":ps,"ani_p":pa,"struct_q":qs,"ani_q":qa,
                     "struct_neg_log10_p":-lps,"ani_neg_log10_p":-lpa,
                     "log2_struct_over_ani":l2,"neg_log10_min_p":-min(lps,lpa),
                     "mode":mode})
    return sorted(rows,key=lambda x:(-x["neg_log10_min_p"],x["cluster"]))


def mode_counts(rows):
    c=Counter(x["mode"] for x in rows)
    return {k:c.get(k,0) for k in ["structural_enriched","ani_enriched","both_enriched","not_enriched"]}


def overlap_jaccard(allrows,subrows,mode):
    a={x["cluster"] for x in allrows if x["mode"]==mode}
    b={x["cluster"] for x in subrows if x["mode"]==mode}
    return len(a&b),len(a|b),len(a&b)/len(a|b) if a|b else 0.0


# These are initialized in main after species names are known.
struct_top=ani_top=both_top=np.array([],dtype=bool)
labels={}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--census",required=True)
    ap.add_argument("--genome-quality",required=True)
    ap.add_argument("--species-metadata",required=True)
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--top-fraction",type=float,default=.05)
    ap.add_argument("--min-shared-tags",type=int,default=250)
    ap.add_argument("--min-observable-fraction",type=float,default=.5)
    ap.add_argument("--min-species-pairs",type=int,default=200)
    ap.add_argument("--correlation-sample",type=int,default=500000)
    ap.add_argument("--seed",type=int,default=42)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    # Representative cluster -> display species.
    rep={}
    with open(args.species_metadata,newline="") as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            rep[row["genome"]]=row
    labels={cl:clean_species((rep.get(cl) or {}).get("species",""),cl) for cl in rep}
    # Disambiguate repeated species names.
    c=Counter(labels.values())
    labels={k:(v if c[v]==1 else f"{v} [{k}]") for k,v in labels.items()}
    names=sorted(labels); code={x:i for i,x in enumerate(names)}

    # Genome -> CheckM/quality metadata and representative cluster.
    qmap={}
    with open(args.genome_quality,newline="") as fh:
        r=csv.reader(fh,delimiter="\t"); next(r)
        for row in r:
            if not row: continue
            genome=row[0]
            qmap[genome]={"comp":float(row[1]),"cont":float(row[2]),
                          "n50":int(row[3]),"contigs":int(row[7]),
                          "score":float(row[5]),"label":row[8],
                          "cluster":row[-1]}

    ani=array.array("d"); structural=array.array("d"); breakpoints=array.array("d")
    mincomp=array.array("d"); maxcont=array.array("d"); minN50=array.array("d")
    maxcontigs=array.array("d"); minscore=array.array("d"); codes=array.array("I")
    seen=dropped=0
    with gzip.open(args.census,"rt") as fh:
        header=fh.readline().rstrip("\n").split("\t")
        ix={x:i for i,x in enumerate(header)}
        ia,ib=ix["genome_A"],ix["genome_B"]
        iani,ist=ix["skani_ani"],ix["structural"]
        ibp,ish,iobs=ix["breakpoints"],ix["shared_tags"],ix["observable_fraction"]
        for line in fh:
            seen+=1; f=line.rstrip("\n").split("\t")
            try:
                shared=int(f[ish]); obs=float(f[iobs]); a=float(f[iani])
                st=float(f[ist]); bp=int(f[ibp])
            except ValueError:
                dropped+=1; continue
            if shared<args.min_shared_tags or obs<args.min_observable_fraction:
                dropped+=1; continue
            qa=qmap.get(f[ia]); qb=qmap.get(f[ib])
            if qa is None or qb is None or qa["cluster"]!=qb["cluster"]:
                dropped+=1; continue
            clcode=code.get(qa["cluster"])
            if clcode is None:
                dropped+=1; continue
            ani.append(a); structural.append(st); breakpoints.append(bp)
            mincomp.append(min(qa["comp"],qb["comp"]))
            maxcont.append(max(qa["cont"],qb["cont"]))
            minN50.append(min(qa["n50"],qb["n50"]))
            maxcontigs.append(max(qa["contigs"],qb["contigs"]))
            minscore.append(min(qa["score"],qb["score"]))
            codes.append(clcode)

    ani=np.frombuffer(ani,dtype=np.float64)
    structural=np.frombuffer(structural,dtype=np.float64)
    breakpoints=np.frombuffer(breakpoints,dtype=np.float64)
    mincomp=np.frombuffer(mincomp,dtype=np.float64)
    maxcont=np.frombuffer(maxcont,dtype=np.float64)
    minN50=np.frombuffer(minN50,dtype=np.float64)
    maxcontigs=np.frombuffer(maxcontigs,dtype=np.float64)
    minscore=np.frombuffer(minscore,dtype=np.float64)
    codes=np.frombuffer(codes,dtype=np.uint32)
    continuous={
        "mincomp":mincomp,"maxcont":maxcont,"minN50":minN50,
        "maxcontigs":maxcontigs,"minscore":minscore,
    }
    N=len(ani); n_top=int(math.floor(args.top_fraction*N))
    rng=np.random.default_rng(args.seed)
    idx=np.argpartition(ani,N-n_top)[N-n_top:]
    ani_top=np.zeros(N,dtype=bool); ani_top[idx]=True
    idx=np.argpartition(structural,N-n_top)[N-n_top:]
    struct_top=np.zeros(N,dtype=bool); struct_top[idx]=True
    both_top=ani_top & struct_top

    # Species enrichment in all qualifying pairs.
    total=np.bincount(codes,minlength=len(names))
    sstruct=np.bincount(codes[struct_top],minlength=len(names))
    sani=np.bincount(codes[ani_top],minlength=len(names))
    sboth=np.bincount(codes[both_top],minlength=len(names))
    allrows=[]
    for i,name in enumerate(names):
        K=int(total[i]);
        if K<args.min_species_pairs: continue
        ks=int(sstruct[i]); ka=int(sani[i]); kb=int(sboth[i])
        ps=float(hypergeom.sf(ks-1,N,K,n_top)); pa=float(hypergeom.sf(ka-1,N,K,n_top))
        lps=float(hypergeom.logsf(ks-1,N,K,n_top)); lpa=float(hypergeom.logsf(ka-1,N,K,n_top))
        l2=math.log2((ks+0.5)/(ka+0.5))
        qs=float(bh_qvalues([ps])[0]); qa=float(bh_qvalues([pa])[0])
        if qs<.05 and l2>=1: mode="structural_enriched"
        elif qa<.05 and l2<=-1: mode="ani_enriched"
        elif qs<.05 and qa<.05: mode="both_enriched"
        else: mode="not_enriched"
        allrows.append({"cluster":name,"species":labels[name],"total_pairs":K,
                        "struct_top_pairs":ks,"ani_top_pairs":ka,"both_top_pairs":kb,
                        "struct_enrichment_ratio":(ks/n_top)/(K/N),
                        "ani_enrichment_ratio":(ka/n_top)/(K/N),
                        "struct_p":ps,"ani_p":pa,"struct_q":qs,"ani_q":qa,
                        "struct_neg_log10_p":-lps,"ani_neg_log10_p":-lpa,
                        "log2_struct_over_ani":l2,"neg_log10_min_p":-min(lps,lpa),
                        "mode":mode})
    allrows.sort(key=lambda x:(-x["neg_log10_min_p"],x["cluster"]))
    with (out/"all_quality_hypermode_species.tsv").open("w") as fh:
        w=csv.DictWriter(fh,fieldnames=list(allrows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(allrows)

    # Quality bins and association tests.
    bins=[("min_completeness_lt70",mincomp<70),
          ("min_completeness_70_90",(mincomp>=70)&(mincomp<90)),
          ("min_completeness_ge90",mincomp>=90),
          ("max_contamination_le5",maxcont<=5),
          ("max_contamination_gt5",maxcont>5),
          ("min_N50_ge10kb",minN50>=10000),
          ("min_N50_lt10kb",minN50<10000)]
    bp2=breakpoints>=2
    qual_rows=[]
    for name,mask in bins:
        n=int(mask.sum()); k=int((mask&bp2).sum())
        if n==0:
            qual_rows.append({"feature":name,"n_pairs":0,
                              "n_breakpoints_ge2":0,
                              "pct_breakpoints_ge2":float("nan"),
                              "pct_all_breakpoints_ge2":100*bp2.mean(),
                              "odds_ratio":float("nan"),"fisher_p":float("nan")})
            continue
        table=np.array([[k,int(n-k)],
                        [int(bp2.sum()-k),int((~bp2).sum()-int(n-k))]])
        OR,p=fisher_exact(table)
        qual_rows.append({"feature":name,"n_pairs":n,"n_breakpoints_ge2":k,
                          "pct_breakpoints_ge2":100*k/n,
                          "pct_all_breakpoints_ge2":100*bp2.mean(),
                          "odds_ratio":OR,"fisher_p":p})
    # Continuous associations.
    sample=rng.choice(N,size=min(args.correlation_sample,N),replace=False)
    cont_rows=[]
    for trait in ["mincomp","maxcont","minN50","maxcontigs","minscore"]:
        values=continuous[trait]
        rho,p=spearmanr(values[sample],breakpoints[sample])
        u,p2=mannwhitneyu(values[bp2],values[~bp2])
        cont_rows.append({"trait":trait,"spearman_rho_vs_breakpoints":rho,
                          "spearman_p":p,
                          "median_breakpoints_ge2":float(np.median(values[bp2])),
                          "median_breakpoints_lt2":float(np.median(values[~bp2])),
                          "mannwhitney_p":p2})
    # Is either top set enriched for high-quality pairs?
    high=mincomp>=90
    for name,mask in [("structural_top",struct_top),("ani_top",ani_top),("both_top",both_top)]:
        a=int((mask&high).sum()); b=int((mask&~high).sum())
        c=int(((~mask)&high).sum()); d=int(((~mask)&~high).sum())
        OR,p=fisher_exact([[a,b],[c,d]])
        qual_rows.append({"feature":f"{name}_high_quality_enrichment",
                          "n_pairs":a+b+c+d,
                          "n_breakpoints_ge2":a,
                          "pct_breakpoints_ge2":100*a/(a+b),
                          "pct_all_breakpoints_ge2":100*high.mean(),
                          "odds_ratio":OR,"fisher_p":p})
    with (out/"quality_association.tsv").open("w") as fh:
        w=csv.DictWriter(fh,fieldnames=list(qual_rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(qual_rows)

    # Repeat enrichment in quality strata.
    strata=[("all_qualified",np.ones(N,dtype=bool)),
            ("high_quality_comp90_cont5",(mincomp>=90)&(maxcont<=5)),
            ("medium_comp70_cont10",(mincomp>=70)&(maxcont<=10)),
            ("high_completeness_only",mincomp>=90),
            ("continuous_N50_ge10kb",minN50>=10000)]
    summary=[]; stratum_rows={}
    for name,mask in strata:
        n=int(mask.sum()); nt=max(1,int(math.floor(args.top_fraction*n)))
        sani=ani[mask]; sst=structural[mask]; scodes=codes[mask]
        if n<1000:
            rows=[]; modes={}
        else:
            sidx=np.argpartition(sani,len(sani)-nt)[len(sani)-nt:]
            m_ani=np.zeros(n,dtype=bool); m_ani[sidx]=True
            sidx=np.argpartition(sst,len(sst)-nt)[len(sst)-nt:]
            m_st=np.zeros(n,dtype=bool); m_st[sidx]=True
            # Patch module globals for the helper.
            globals()["ani_top"]=m_ani; globals()["struct_top"]=m_st
            globals()["both_top"]=m_ani & m_st
            rows=species_enrichment(np.ones(n,dtype=bool),scodes,names,n,nt,args.min_species_pairs)
            modes=mode_counts(rows)
            sp=(out/f"hypermode_species_{name}.tsv")
            with sp.open("w") as fh:
                w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
                w.writeheader(); w.writerows(rows)
            stratum_rows[name]=rows
        summary.append({"stratum":name,"pairs":n,"species_tested":len(rows),
                        "top_pairs_per_axis":nt if n>=1000 else 0,
                        "structural_enriched":modes.get("structural_enriched",0),
                        "ani_enriched":modes.get("ani_enriched",0),
                        "both_enriched":modes.get("both_enriched",0),
                        "not_enriched":modes.get("not_enriched",0)})
    with (out/"quality_stratified_summary.tsv").open("w") as fh:
        w=csv.DictWriter(fh,fieldnames=list(summary[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(summary)

    # Stability of the main modes: high-quality versus all.
    stability=[]
    for mode in ["structural_enriched","ani_enriched","both_enriched"]:
        inter,jacc,frac=overlap_jaccard(allrows,stratum_rows.get("high_quality_comp90_cont5",[]),mode)
        stability.append({"mode":mode,"all_species":sum(x["mode"]==mode for x in allrows),
                          "high_quality_species":sum(x["mode"]==mode for x in stratum_rows.get("high_quality_comp90_cont5",[])),
                          "intersection":inter,"jaccard":frac})
    # Named robust high-quality candidates.
    robust=[]
    hq={x["cluster"]:x for x in stratum_rows.get("high_quality_comp90_cont5",[])}
    for x in allrows:
        y=hq.get(x["cluster"])
        if y and x["mode"]!="not_enriched" and y["mode"]==x["mode"]:
            robust.append({"cluster":x["cluster"],"species":x["species"],
                           "mode":x["mode"],"all_total_pairs":x["total_pairs"],
                           "hq_total_pairs":y["total_pairs"],
                           "all_log2":x["log2_struct_over_ani"],
                           "hq_log2":y["log2_struct_over_ani"],
                           "all_min_p":x["neg_log10_min_p"],
                           "hq_min_p":y["neg_log10_min_p"]})
    with (out/"high_quality_stable_species.tsv").open("w") as fh:
        if robust:
            w=csv.DictWriter(fh,fieldnames=list(robust[0]),delimiter="\t",lineterminator="\n")
            w.writeheader(); w.writerows(sorted(robust,key=lambda x:(x["mode"],-x["hq_min_p"])))

    report={"pairs_seen":seen,"pairs_dropped_quality_or_cluster":dropped,
            "qualifying_pairs":N,"top_pairs_per_axis":n_top,
            "quality_distribution":{"high_comp90_cont5":int(((mincomp>=90)&(maxcont<=5)).sum()),
                                    "medium_comp70_cont10":int(((mincomp>=70)&(maxcont<=10)).sum()),
                                    "low":int(((mincomp<70)|(maxcont>10)).sum())},
            "quality_association":qual_rows,
            "quality_stratified_modes":summary,
            "mode_stability_high_quality":stability,
            "robust_species_n":len(robust),
            "spearman_mincomp_vs_breakpoints":float(spearmanr(mincomp[sample],breakpoints[sample]).statistic),
            "spearman_minscore_vs_breakpoints":float(spearmanr(minscore[sample],breakpoints[sample]).statistic)}
    with (out/"quality_stratified_report.json").open("w") as fh:
        json.dump(report,fh,indent=2)

    # Compact human-readable report.
    with (out/"QUALITY_STRATIFIED_REPORT.md").open("w") as fh:
        fh.write("# HROM CheckM/quality sensitivity analysis\n\n")
        fh.write(f"- pairs seen: {seen:,}\n- qualifying pairs: {N:,}\n")
        fh.write(f"- top 5% pairs per axis: {n_top:,}\n\n")
        fh.write("## Quality-stratified mode counts\n\n")
        fh.write("| stratum | pairs | species tested | structural | ANI | both | not |\n|---|---:|---:|---:|---:|---:|---:|\n")
        for x in summary:
            fh.write(f"| {x['stratum']} | {x['pairs']:,} | {x['species_tested']:,} | "
                     f"{x['structural_enriched']:,} | {x['ani_enriched']:,} | "
                     f"{x['both_enriched']:,} | {x['not_enriched']:,} |\n")
        fh.write("\n## High-quality mode stability\n\n| mode | all | high-quality | intersection | Jaccard |\n|---|---:|---:|---:|---:|\n")
        for x in stability:
            fh.write(f"| {x['mode']} | {x['all_species']:,} | {x['high_quality_species']:,} | "
                     f"{x['intersection']:,} | {x['jaccard']:.3f} |\n")
        fh.write("\nSee TSV files for exact odds ratios, P values, per-stratum species and robust candidates.\n")

    print(json.dumps(report,indent=2,default=str))

if __name__=="__main__":
    main()
