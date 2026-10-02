#!/usr/bin/env python3
"""Joint OAPGC structural + ANI analysis by genome site.

Inputs are intentionally generic so the public OAPGC annotation package can be
mapped to a long genome metadata table before analysis.  The expected genome
metadata has genome_id, sample_id/site (or a direct site), and optional quality
columns.  The script emits pair-level axes and per-species/site summaries.
"""
import argparse, collections, gzip, json, math, random, statistics
from pathlib import Path


def read_table(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        yield fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield line.rstrip("\n").split("\t")


def is_gzip(path): return str(path).endswith(".gz")


def opener(path): return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt")


def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def site_category(x):
    x = str(x).lower()
    if any(k in x for k in ("oral", "saliva", "dental", "plaque", "tongue", "mouth")): return "oral"
    if any(k in x for k in ("airway", "airways", "lung", "sputum", "bronch", "nasopharyn", "oropharyn", "tracheal", "respiratory")): return "airway"
    return "other"


def load_genome_meta(path):
    """Return accession -> metadata dict; also infer sample_id by prefix if absent."""
    out = {}
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        lower = [x.lower() for x in header]
        required = {"genome_id", "genome", "hqmags", "accession"}
        ig = next((lower.index(x) for x in required if x in lower), None)
        if ig is None: raise SystemExit(f"no genome_id column in {path}")
        isite = next((lower.index(x) for x in ("site", "sampletype", "sample_type", "habitat", "body_site") if x in lower), None)
        isample = next((lower.index(x) for x in ("sample_id", "sample", "run", "used_name") if x in lower), None)
        icomp = next((lower.index(x) for x in ("checkm_completeness", "completeness") if x in lower), None)
        icont = next((lower.index(x) for x in ("checkm_contamination", "contamination") if x in lower), None)
        in50 = next((lower.index(x) for x in ("n50_contigs", "n50") if x in lower), None)
        icluster = next((lower.index(x) for x in ("species_cluster", "sgb", "cluster") if x in lower), None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(x for x in [ig,isite,isample,icomp,icont,in50,icluster] if x is not None): continue
            gid = f[ig]; sample = f[isample] if isample is not None else gid.split("_",1)[0]
            site = f[isite] if isite is not None else ""
            out[gid] = {"sample_id":sample, "site_raw":site, "site":site_category(site),
                        "completeness":num(f[icomp]) if icomp is not None else None,
                        "contamination":num(f[icont]) if icont is not None else None,
                        "n50":num(f[in50]) if in50 is not None else None,
                        "cluster":f[icluster] if icluster is not None else None}
    return out


def quality_pass(m, min_comp, max_cont, min_n50):
    for gid, x in m.items():
        x["quality_pass"] = ((x["completeness"] is None or x["completeness"] >= min_comp) and
                             (x["contamination"] is None or x["contamination"] <= max_cont) and
                             (x["n50"] is None or x["n50"] >= min_n50))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--structural",required=True)
    ap.add_argument("--ani-dir",required=True)
    ap.add_argument("--genome-metadata",required=True)
    ap.add_argument("--out-prefix",required=True)
    ap.add_argument("--min-comp",type=float,default=0)
    ap.add_argument("--max-cont",type=float,default=100)
    ap.add_argument("--min-n50",type=float,default=0)
    ap.add_argument("--sample-cap",type=int,default=100000)
    a=ap.parse_args()
    meta=load_genome_meta(a.genome_metadata)
    quality_pass(meta,a.min_comp,a.max_cont,a.min_n50)
    # Structural output is grouped by cluster because merge processes sorted task files.
    current=None; ani={}
    pair_path=Path(a.out_prefix+"_pair_axes.tsv.gz")
    summary_path=Path(a.out_prefix+"_species_site_summary.tsv")
    audit_path=Path(a.out_prefix+"_audit.tsv")
    reservoirs=collections.defaultdict(lambda: collections.defaultdict(list))
    exact=collections.defaultdict(lambda: collections.defaultdict(int))
    counts=collections.Counter()
    with gzip.open(pair_path,"wt") as out:
        out.write("cluster\tgenome_A\tgenome_B\tsite_pair\tANI\tSNP_distance\tbreakpoints\tge2\tbreakpoint_density\tobservable_adjacencies\tstructural\n")
        with opener(a.structural) as fh:
            header=next(fh).rstrip("\n").split("\t")
            ix={x:i for i,x in enumerate(header)}
            for line in fh:
                f=line.rstrip("\n").split("\t"); ga,gb=f[ix["genome_A"]],f[ix["genome_B"]]
                ma,mb=meta.get(ga),meta.get(gb)
                if not ma or not mb: counts["missing_genome_meta"]+=1; continue
                if ma["quality_pass"] is False or mb["quality_pass"] is False: counts["quality_filtered"]+=1; continue
                cl=ma.get("cluster") or mb.get("cluster")
                if not cl or ma.get("cluster") and mb.get("cluster") and ma["cluster"]!=mb["cluster"]: counts["cross_cluster"]+=1; continue
                if cl!=current:
                    current=cl; ani={}
                    p=Path(a.ani_dir)/f"{cl}.ani.tsv.gz"
                    if p.exists():
                        with gzip.open(p,"rt") as af:
                            next(af)
                            for al in af:
                                x=al.rstrip("\n").split("\t")
                                if len(x)>=3: ani[frozenset((x[0],x[1]))]=num(x[2])
                v=ani.get(frozenset((ga,gb)))
                if v is None: counts["missing_ani"]+=1; continue
                sa,sb=ma["site"],mb["site"]
                site_pair="-".join(sorted((sa,sb)))
                bp=num(f[ix["breakpoints"]]); density=num(f[ix["breakpoint_density"]]); obs=num(f[ix["observable_adjacencies"]]); struct=num(f[ix["structural"]])
                if bp is None: counts["missing_breakpoints"]+=1; continue
                ge2=int(bp>=2); snp=max(0.0,100.0-v)
                out.write("\t".join(map(str,[cl,ga,gb,site_pair,v,snp,bp,ge2,density if density is not None else "",obs if obs is not None else "",struct if struct is not None else ""]))+"\n")
                key=(cl,site_pair)
                exact[key]["pairs"]+=1; exact[key]["ge2"]+=ge2
                for axis,val in [("ANI",v),("SNP_distance",snp),("breakpoints",bp)]:
                    if val is None: continue
                    exact[key][f"{axis}_sum"]+=val; exact[key][f"{axis}_n"]+=1
                    arr=reservoirs[key][axis]
                    if len(arr)<a.sample_cap: arr.append(val)
                    else:
                        exact[key]["_seen"]=exact[key].get("_seen",0)+1
                        j=random.Random(exact[key]["_seen"]).randrange(exact[key]["_seen"])
                        if j<a.sample_cap: arr[j]=val
                if density is not None:
                    exact[key]["density_sum"]+=density; exact[key]["density_n"]+=1
                if obs is not None and obs>0:
                    exact[key]["normalized_bp_sum"]+=bp/obs; exact[key]["normalized_bp_n"]+=1
                counts["written"]+=1
    with open(summary_path,"w") as out:
        out.write("cluster\tsite_pair\tpairs\tpairs_ge2\tfraction_ge2\tmedian_ANI\tmean_ANI\tmedian_SNP_distance\tmean_SNP_distance\tmedian_breakpoints\tmean_breakpoints\tmedian_breakpoint_density\tmedian_normalized_breakpoints\n")
        for key in sorted(exact):
            cl,sp=key; x=exact[key]; r=reservoirs[key]
            def med(axis): return statistics.median(r[axis]) if r[axis] else ""
            def mean(axis,n): return x.get(f"{axis}_sum",0)/n if n else ""
            out.write("\t".join(map(str,[cl,sp,x["pairs"],x["ge2"],x["ge2"]/x["pairs"] if x["pairs"] else "",
                                          med("ANI"),mean("ANI",x.get("ANI_n",0)),med("SNP_distance"),mean("SNP_distance",x.get("SNP_distance_n",0)),
                                          med("breakpoints"),mean("breakpoints",x.get("breakpoints_n",0)),med("breakpoint_density"),
                                          med("normalized_breakpoints")]))+"\n")
    with open(audit_path,"w") as out:
        out.write("metric\tvalue\n")
        for k in sorted(counts): out.write(f"{k}\t{counts[k]}\n")
    print(json.dumps({"pair_rows":counts["written"],"metadata":len(meta),"output":str(pair_path)},indent=2))

if __name__=="__main__": main()
