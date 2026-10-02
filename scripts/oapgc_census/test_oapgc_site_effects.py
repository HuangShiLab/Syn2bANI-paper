#!/usr/bin/env python3
"""Species-level oral-airway effect tests for SNP and SV axes."""
import argparse,csv,gzip,math,random
from pathlib import Path
from collections import defaultdict,Counter
from scipy.stats import mannwhitneyu,fisher_exact,binomtest
import numpy as np


def bh(p):
    n=len(p); order=sorted(range(n),key=lambda i:p[i]); q=[0]*n; prev=1
    for rank,idx in reversed(list(enumerate(order,1))):
        val=min(prev,p[idx]*n/rank); q[idx]=val; prev=val
    return q


def read_eligible(summary,min_oa=30,min_within=30,min_pairs=100):
    out={}
    with open(summary) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            cl = r['cluster']
            sp = r['site_pair']
            oa = oo = aa = 0
            if sp=='airway-oral': oa=int(r['pairs'])
            elif sp=='oral-oral': oo=int(r['pairs'])
            elif sp=='airway-airway': aa=int(r['pairs'])
            old=out.setdefault(cl,{'oa':0,'within':0,'pairs':0}); old['oa']+=oa; old['within']+=oo+aa; old['pairs']+=oa+oo+aa
    return {c:x for c,x in out.items() if x['oa']>=min_oa and x['within']>=min_within and x['pairs']>=min_pairs}


def run_pass(args,label,minani,minn50):
    meta={}
    with open(args.metadata) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            try: n50=float(r['n50'])
            except: n50=float('inf')
            meta[r['genome_id']]=n50
    eligible=read_eligible(args.summary,args.min_oa,args.min_within,args.min_pairs)
    store=defaultdict(lambda: {'within_snp':[],'oa_snp':[],'within_bp':[],'oa_bp':[],
                               'within_n':0,'within_ge2':0,'oa_n':0,'oa_ge2':0})
    counts=Counter()
    with gzip.open(args.pairs,'rt') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            cl=row['cluster']; sp=row['site_pair']
            if sp not in ('oral-oral','airway-airway','airway-oral'): continue
            if sp=='airway-oral': group='oa'
            else: group='within'
            a,b=meta.get(row['genome_A']),meta.get(row['genome_B'])
            if not isinstance(a,float) or not isinstance(b,float): continue
            minn=min(a,b)
            if minn<minn50: counts[f'{label}_filtered_n50']+=1; continue
            try: ani=float(row['ANI']); snp=float(row['SNP_distance']); bp=float(row['breakpoints'])
            except: counts[f'{label}_missing_metric']+=1; continue
            if ani<minani: counts[f'{label}_filtered_ani']+=1; continue
            if cl not in eligible: counts[f'{label}_ineligible_species']+=1; continue
            ge2=int(float(row['ge2']))
            d=store[cl]; counts[f'{label}_used']+=1
            if group=='oa': d['oa_snp'].append(snp); d['oa_bp'].append(bp); d['oa_n']+=1; d['oa_ge2']+=ge2
            else: d['within_snp'].append(snp); d['within_bp'].append(bp); d['within_n']+=1; d['within_ge2']+=ge2
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    rows=[]; auc=[]
    for cl,d in store.items():
        if len(d['oa_snp'])<args.min_oa or len(d['within_snp'])<args.min_within: continue
        res={}
        for metric,a,b in [('SNP',d['oa_snp'],d['within_snp']),('BP',d['oa_bp'],d['within_bp'])]:
            u,p=mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
            aucv=u/(len(a)*len(b)); res[metric]=(p,2*aucv-1,aucv,u,len(a),len(b))
        table=[[d['oa_ge2'],d['oa_n']-d['oa_ge2']],[d['within_ge2'],d['within_n']-d['within_ge2']]]
        odds,pge=fisher_exact(table)
        rows.append({'cluster':cl,'oa_pairs':len(d['oa_snp']),'within_pairs':len(d['within_snp']),
                     'SNP_p':res['SNP'][0],'SNP_cliffs_delta':res['SNP'][1],'SNP_auc':res['SNP'][2],
                     'BP_p':res['BP'][0],'BP_cliffs_delta':res['BP'][1],'BP_auc':res['BP'][2],
                     'SV_ge2_OR':odds,'SV_ge2_p':pge})
        auc.append(res['SNP'][2]); auc.append(res['BP'][2])
    for metric in ['SNP_p','BP_p','SV_ge2_p']:
        vals=[float(r[metric]) for r in rows]; q=bh(vals)
        for r,qq in zip(rows,q): r[metric.replace('_p','_q')]=qq
    with open(out/f'species_effects_{label}.tsv','w',newline='') as fh:
        if not rows: fh.write('cluster\n'); return
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    positive_snp=sum(r['SNP_auc']>0.5 for r in rows); positive_bp=sum(r['BP_auc']>0.5 for r in rows)
    binsnp=binomtest(positive_snp,len(rows),0.5).pvalue; binbp=binomtest(positive_bp,len(rows),0.5).pvalue
    with open(out/f'overall_{label}.tsv','w') as fh:
        fh.write('metric\tvalue\n')
        fh.write(f'species_tested\t{len(rows)}\n')
        fh.write(f'median_SNP_auc\t{np.median([r["SNP_auc"] for r in rows])}\n')
        fh.write(f'median_BP_auc\t{np.median([r["BP_auc"] for r in rows])}\n')
        fh.write(f'species_OA_higher_SNP\t{positive_snp}\n')
        fh.write(f'species_OA_higher_BP\t{positive_bp}\n')
        fh.write(f'sign_test_SNP_p\t{binsnp}\n'); fh.write(f'sign_test_BP_p\t{binbp}\n')
        fh.write(f'n_OA_pairs\t{sum(r["oa_pairs"] for r in rows)}\n')
        fh.write(f'n_within_pairs\t{sum(r["within_pairs"] for r in rows)}\n')
    print(label,'species',len(rows),'medianAUC',np.median([r['SNP_auc'] for r in rows]),np.median([r['BP_auc'] for r in rows]))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--pairs',required=True); ap.add_argument('--summary',required=True); ap.add_argument('--metadata',required=True)
    ap.add_argument('--out',required=True); ap.add_argument('--min-oa',type=int,default=30); ap.add_argument('--min-within',type=int,default=30); ap.add_argument('--min-pairs',type=int,default=100)
    a=ap.parse_args()
    run_pass(a,'all_pairs',0,0)
    run_pass(a,'n50_ge10kb',0,10000)
    run_pass(a,'n50_ge100kb',0,100000)
    run_pass(a,'recent_ani_ge99p9',99.9,0)

if __name__=='__main__': main()
