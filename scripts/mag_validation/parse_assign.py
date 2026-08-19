#!/usr/bin/env python3
"""parse_assign.py <sample_bins.fa> <paf> <genome_lengths.tsv> <species_clusters.tsv> <out.tsv>
Contig-level truth assignment from minimap2 PAF (with -c CIGAR tags).
Thresholds: query-coverage >= ASSIGN_MIN_COV, identity >= ASSIGN_MIN_ID (env or defaults).
Per bin: majority genome, alignment-based completeness of the majority genome,
contamination split (same-species vs cross-species) and contamination class.
"""
import os
import re
import sys
from collections import defaultdict

fa, paf, lengths_tsv, clusters_tsv, out_tsv = sys.argv[1:6]
MIN_COV = float(os.environ.get("ASSIGN_MIN_COV", "0.80"))
MIN_ID = float(os.environ.get("ASSIGN_MIN_ID", "0.95"))

# contig -> (binid, length)
ctg_bin, ctg_len = {}, {}
name, cur = None, 0
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                ctg_len[name] = cur
            name = line[1:].split()[0]
            ctg_bin[name] = name.split("|")[0]
            cur = 0
        else:
            cur += len(line.strip())
if name is not None:
    ctg_len[name] = cur

glen = {}
with open(lengths_tsv) as f:
    next(f)
    for l in f:
        g, L = l.rstrip("\n").split("\t")
        glen[g] = int(L)

cluster = {}
with open(clusters_tsv) as f:
    next(f)
    for l in f:
        g, c = l.rstrip("\n").split("\t")
        cluster[g] = c

CIG_RE = re.compile(r"(\d+)([=XIDM])")

# per contig: list of (genome, nmatch, qspan, tstart, tend)
hits = defaultdict(list)
with open(paf) as f:
    for l in f:
        p = l.rstrip("\n").split("\t")
        qname, qlen, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
        tname, ts, te, nmatch = p[5], int(p[7]), int(p[8]), int(p[9])
        if qname not in ctg_bin:
            continue
        genome = tname.split("|")[0]
        cov = (qe - qs) / qlen
        if cov < MIN_COV:
            continue
        cg = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
        ident = None
        if cg:
            eq = xi = ins = dele = mop = 0
            for n, op in CIG_RE.findall(cg):
                n = int(n)
                if op == "=": eq += n
                elif op == "X": xi += n
                elif op == "M": mop += n
                elif op == "I": ins += n
                else: dele += n
            if eq or xi:
                # --eqx CIGAR: exact identity
                denom = eq + xi + ins + dele
                ident = eq / denom if denom else 0.0
            else:
                # plain-M CIGAR: derive from NM (mismatches + indel bases)
                nm = next((int(t[5:]) for t in p[12:] if t.startswith("NM:i:")), None)
                if nm is not None:
                    aln = mop + ins + dele
                    ident = (aln - nm) / aln if aln else 0.0
        if ident is None:
            ident = nmatch / int(p[10])
        if ident < MIN_ID:
            continue
        hits[qname].append((genome, nmatch, qe - qs, ts, te))

bin_stats = {}
for ctg, L in ctg_len.items():
    b = bin_stats.setdefault(ctg_bin[ctg], {
        "total_bp": 0, "assigned_bp": 0, "genome_bp": defaultdict(int),
        "maj_intervals": defaultdict(list), "n_assigned": 0, "n_multihit": 0})
    b["total_bp"] += L
    h = hits.get(ctg)
    if not h:
        continue
    b["n_assigned"] += 1
    if len({x[0] for x in h}) > 1:
        b["n_multihit"] += 1
    best = max(h, key=lambda x: x[1])
    b["assigned_bp"] += L
    b["genome_bp"][best[0]] += L
    for g, _, _, ts, te in h:
        b["maj_intervals"][g].append((ts, te))

def union_len(iv):
    tot = 0
    cs = ce = None
    for s, e in sorted(iv):
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if ce is not None:
        tot += ce - cs
    return tot

with open(out_tsv, "w") as o:
    o.write("bin\ttotal_bp\tassigned_bp\tunassigned_bp\tmajority_genome\tmajority_bp\t"
            "maj_comp_est\tsame_sp_bp\tcross_sp_bp\tcontam_bp\tcontam_pct\tclass\t"
            "n_contigs_assigned\tn_contigs_multihit\n")
    for binid, b in sorted(bin_stats.items()):
        if not b["genome_bp"]:
            o.write(f"{binid}\t{b['total_bp']}\t0\t{b['total_bp']}\tNA\t0\t0\t0\t0\t0\t0\t"
                    f"unassigned\t0\t0\n")
            continue
        maj, maj_bp = max(b["genome_bp"].items(), key=lambda kv: kv[1])
        comp = union_len(b["maj_intervals"][maj]) / glen.get(maj, 1)
        same = sum(v for g, v in b["genome_bp"].items() if g != maj and cluster.get(g) == cluster.get(maj))
        cross = sum(v for g, v in b["genome_bp"].items() if g != maj and cluster.get(g) != cluster.get(maj))
        contam = same + cross
        assigned = b["assigned_bp"]
        if maj_bp / assigned >= 0.95:
            cls = "clean"
        elif same >= cross:
            cls = "strain-mixed"
        else:
            cls = "cross-species"
        o.write(f"{binid}\t{b['total_bp']}\t{assigned}\t{b['total_bp']-assigned}\t{maj}\t{maj_bp}\t"
                f"{comp:.4f}\t{same}\t{cross}\t{contam}\t{100.0*contam/assigned:.3f}\t{cls}\t"
                f"{b['n_assigned']}\t{b['n_multihit']}\n")
