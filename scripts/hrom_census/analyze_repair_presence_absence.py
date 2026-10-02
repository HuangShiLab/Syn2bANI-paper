#!/usr/bin/env python3
"""Summarize DNA-repair gene presence/absence in HROM ANI-enriched species."""
import csv
from collections import defaultdict, Counter

CAND="/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified/hrom_ani_candidate_clusters.txt"
HITS="/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified/hrom_ani_enriched_repair_genes.tsv"
META="/lustre1/g/aos_shihuang/databases/HROM/census/hrom_genome_metadata.tsv"
OUT="/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified/hrom_ani_enriched_repair_presence_absence.tsv"
GENES=["MutS","MutL","MutH","MutY","MutM","MutT","UvrD","DnaQ","Ogt_Ada","AlkB"]
KOS={"MutS":{"K03515"},"MutL":{"K03516"},"MutH":{"K03517"},"MutY":{"K03518"},"MutM":{"K03514"},"MutT":{"K03519"},"UvrD":{"K03657"},"DnaQ":{"K02320"},"Ogt_Ada":{"K03555"},"AlkB":{"K03554"}}

def classify(row):
    kos=set(x.replace("ko:","") for x in row.get("KEGG_ko","").split(",") if x and x!="-")
    text=" ".join([row.get("Description",""),row.get("Preferred_name",""),row.get("PFAMs","")]).lower()
    pref=row.get("Preferred_name","").lower(); desc=row.get("Description","").lower(); pfam=row.get("PFAMs","").lower()
    out=[]
    if "K03515" in kos or "mismatch repair protein muts" in desc or (pref.startswith("muts") and pref!="muts2") or any(x in pfam for x in ("MutS_I","MutS_II","MutS_III","MutS_IV")): out.append("MutS")
    if "K03516" in kos or "mismatch repair protein mutl" in desc or "mutl" in pref or "mutl" in pfam: out.append("MutL")
    if "K03517" in kos or "mismatch repair endonuclease muth" in desc or "muth" in pref: out.append("MutH")
    if "K03518" in kos or "a/g-specific adenine glycosylase" in desc or "muty" in pref: out.append("MutY")
    if "K03514" in kos or "formamidopyrimidine-dna glycosylase" in desc or "mutm" in pref or "fpg" in pfam: out.append("MutM")
    if "K03519" in kos or "8-oxo-dgtp" in desc or "mutt" in pref: out.append("MutT")
    if "K03657" in kos or "uvrd" in pref or "uvrd" in pfam: out.append("UvrD")
    if "K02320" in kos or "dna polymerase iii" in desc and "epsilon" in desc or "dna_pol3_epsilon" in pfam: out.append("DnaQ")
    if "K03555" in kos or "methylated-dna-protein-cysteine methyltransferase" in desc or "ogt" in pref or "ada" in pref: out.append("Ogt_Ada")
    if "K03554" in kos or "alkb" in pref or ("2-oxoglutarate" in desc and "dioxygenase" in desc): out.append("AlkB")
    return out

cands=set(x.strip() for x in open(CAND) if x.strip())
genomes_by_cluster=defaultdict(set)
with open(META) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["cluster"] in cands: genomes_by_cluster[r["cluster"]].add(r["genome"])
hits=defaultdict(lambda:defaultdict(set))
with open(HITS) as f:
    for row in csv.DictReader(f,delimiter="\t"):
        q=row["#query"]
        for cl in cands:
            if q.startswith(cl+"_"):
                genome="_".join(q.split("_")[:3])
                genomes_by_cluster[cl].add(genome)
                for gene in classify(row): hits[cl][gene].add(genome)
                break
with open(OUT,"w",newline="") as o:
    w=csv.writer(o,delimiter="\t",lineterminator="\n")
    w.writerow(["cluster","genomes"]+sum(([g,"genomes_with_"+g,"fraction_with_"+g] for g in GENES),[])+["complete_mmr","missing_core_mmr"])
    for cl in sorted(cands):
        gs=genomes_by_cluster[cl]; vals=[]
        present={}
        for g in GENES:
            n=len(hits[cl].get(g,set())); present[g]=n>0; ngen=len(gs); vals += [n,n,n/ngen if ngen else 0]
        core=["MutS","MutL","MutH"]; mmr=all(present[g] for g in core)
        w.writerow([cl,len(gs)]+vals+["yes" if mmr else "no","yes" if not mmr else "no"])
print("clusters",len(cands),"genomes",sum(len(x) for x in genomes_by_cluster.values()))
print("missing_core_mmr",sum(1 for cl in cands if not all(hits[cl].get(g,set()) for g in ["MutS","MutL","MutH"])))
