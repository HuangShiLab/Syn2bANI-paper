#!/usr/bin/env python3
"""Final validation: strand structure + cs-string indel verification per pair."""
import csv
import os
import re
from collections import Counter

WD = "/Users/macstudio/Downloads/Syn2bANI-paper/case_studies/fda_argos"
os.chdir(WD)

def fasta_lens(path):
    lens, name, n = {}, None, 0
    for ln in open(path):
        if ln.startswith(">"):
            if name: lens[name] = n
            name, n = ln[1:].split()[0], 0
        else:
            n += len(ln.strip())
    if name: lens[name] = n
    return lens

def merge(ivs):
    ivs = sorted(ivs); out = []
    for s, e in ivs:
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out

def overlap_len(s1, e1, ivs):
    return sum(max(0, min(e1, e2) - max(s1, s2)) for s2, e2 in ivs)

CS_RE = re.compile(r"([:*+\-~])([^:*+\-~]+)")

def cs_indels(cs, tstart, strand):
    """Walk a cs string; return ref-anchored indels >=1000 bp: (rpos, length, 'D'|'I')."""
    out, rpos = [], tstart
    for op, arg in CS_RE.findall(cs):
        if op == ":":
            rpos += int(arg)
        elif op == "*":
            rpos += 1
        elif op == "-":          # present in ref, absent in query = DEL
            n = len(arg)
            if n >= 1000: out.append((rpos, n, "D"))
            rpos += n
        elif op == "+":          # present in query, absent in ref = INS
            n = len(arg)
            if n >= 1000: out.append((rpos, n, "I"))
        elif op == "~":          # long form ~ac10gt50ac
            m = re.match(r"([a-z]{2})(\d+)([a-z]{2})", arg)
            if m:
                n = int(m.group(2))
                rpos += n + 4
    return out

pairs = list(csv.DictReader(open("pairs99.9_ecoli.tsv"), delimiter="\t"))
meta = {r["assembly_acc"]: r for r in csv.DictReader(open("metadata_ecoli.tsv"), delimiter="\t")}
rows = []
for p in pairs:
    ref, qry = sorted([p["acc1"], p["acc2"]])
    base = f"{ref}__{qry}"
    rlens = fasta_lens(f"genomes/ecoli/{ref}.fna")
    chrom = max(rlens, key=rlens.get)
    chrom_len = rlens[chrom]

    alns, all_ref_bases = [], 0
    indels_cs = []   # (rpos, len, D/I) on chromosome
    for ln in open(f"validate/pafs_cs/{base}.paf"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12: continue
        strand, tname = f[4], f[5]
        ts, te = int(f[7]), int(f[8])
        all_ref_bases += te - ts
        cs = next((t[3:] for t in f[12:] if t.startswith("cs:Z:")), "")
        if tname == chrom:
            alns.append({"strand": strand, "ts": ts, "te": te})
            indels_cs.extend(cs_indels(cs, ts, strand))

    plus = sum(a["te"] - a["ts"] for a in alns if a["strand"] == "+")
    minus = sum(a["te"] - a["ts"] for a in alns if a["strand"] == "-")
    n_blocks = len(alns)
    minus_frac = minus / (plus + minus) if (plus + minus) else 0.0
    chrom_share = (plus + minus) / all_ref_bases if all_ref_bases else 0.0
    cls = ("global_flip_or_rotation" if minus_frac > 0.8
           else "partial_inversion" if minus_frac > 0.2 else "collinear_local")
    flags = []
    if chrom_len < 2_000_000: flags.append("draft_ref_no_dominant_chromosome")
    if chrom_share < 0.5: flags.append("chromosome_not_dominant")
    if n_blocks > 30: flags.append("many_blocks")

    minus_zone = merge([(a["ts"], a["te"]) for a in alns if a["strand"] == "-"])
    plus_zone = merge([(a["ts"], a["te"]) for a in alns if a["strand"] == "+"])
    boundaries = sorted({b for iv in (minus_zone + plus_zone) for b in iv})

    n_sv = n_off = 0
    n_rearr_artifact, n_rearr_strandok = 0, 0
    n_indel_csver, n_indel_small_or_unver, n_indel_boundary = 0, 0, 0
    type_counter = Counter()
    for ln in open(f"struct_pairs99.9/{base}.bed"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5: continue
        etype = f[3].split("_")[0]
        n_sv += 1
        type_counter[etype] += 1
        if f[0] != chrom:
            n_off += 1   # off-chromosome (plasmid/contig) event: genuine content diff
            continue
        s, e = int(f[1]), int(f[2])
        span = e - s
        if etype in ("INV", "TRA"):
            ov_m = overlap_len(s, e, minus_zone)
            ov_p = overlap_len(s, e, plus_zone)
            if minus_frac > 0.8:
                if ov_p >= 0.5 * span: n_rearr_strandok += 1
                else: n_rearr_artifact += 1
            else:
                if ov_m >= 0.5 * span: n_rearr_strandok += 1
                else:
                    near_b = any(abs(s - b) < 2000 or abs(e - b) < 2000 for b in boundaries)
                    if near_b and n_blocks > 4: n_rearr_artifact += 1
                    else: n_rearr_strandok += 1  # keep (chain-level call, no contradicting evidence)
        else:  # INS / DEL / DUP
            want = "D" if etype in ("DEL", "DUP") else "I"
            hit = any(abs(rp - s) < 3000 or abs(rp - e) < 3000 or (s <= rp <= e)
                      for rp, ln2, k in indels_cs if k == want and ln2 >= min(1000, 0.5 * span))
            near_b = any(abs(s - b) < 2000 or abs(e - b) < 2000 for b in boundaries)
            if hit: n_indel_csver += 1
            elif near_b and n_blocks <= 4: n_indel_boundary += 1
            else: n_indel_small_or_unver += 1

    genuine = n_rearr_strandok + n_indel_csver + n_off
    a1 = set(x.split("=")[0] for x in meta.get(ref, {}).get("AMR_genotypes", "").split(",") if x)
    a2 = set(x.split("=")[0] for x in meta.get(qry, {}).get("AMR_genotypes", "").split(",") if x)
    rows.append({
        "pair": base, "ref": ref, "query": qry,
        "strain_ref": p["strain1"] if p["acc1"] == ref else p["strain2"],
        "strain_query": p["strain2"] if p["acc1"] == ref else p["strain1"],
        "ani": p["ani"], "chrom_len": chrom_len, "n_blocks": n_blocks,
        "minus_fraction": round(minus_frac, 4), "chrom_share": round(chrom_share, 3),
        "class": cls, "flags": ";".join(flags) or "-",
        "n_sv_total": n_sv,
        "n_rearr_artifact": n_rearr_artifact,
        "n_rearr_strand_confirmed": n_rearr_strandok,
        "n_indel_cs_verified": n_indel_csver,
        "n_indel_small_or_unverified": n_indel_small_or_unver,
        "n_indel_boundary_artifact": n_indel_boundary,
        "n_off_chromosome": n_off,
        "n_sv_local_genuine_estimate": genuine,
        "n_cs_indels_1kb": len(indels_cs),
        "n_amr_diff_genes": len(a1 ^ a2),
        "amr_diff": ";".join(sorted(a1 ^ a2)) or "-",
        "event_types": ";".join(f"{t}:{c}" for t, c in sorted(type_counter.items())) or "none",
    })

with open("validate/validation_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

cc = Counter(r["class"] for r in rows)
print("classification:", dict(cc))
print("pairs with >=1 genuine (strand-confirmed rearr / cs-verified indel>=1kb / off-chrom):",
      sum(1 for r in rows if r["n_sv_local_genuine_estimate"] > 0), "/ 103")
for c in ["global_flip_or_rotation", "partial_inversion", "collinear_local"]:
    sub = [r for r in rows if r["class"] == c]
    g = sum(1 for r in sub if int(r["n_sv_local_genuine_estimate"]) > 0)
    tv = sum(int(r["n_indel_cs_verified"]) for r in sub)
    ts = sum(int(r["n_rearr_strand_confirmed"]) for r in sub)
    ta = sum(int(r["n_rearr_artifact"]) for r in sub)
    print(f"  {c}: n={len(sub)}, with-genuine={g}, strand-conf-rearr={ts}, cs-ver-indel={tv}, rearr-artifact={ta}")
print("\ntop flipped pairs by n_sv (deep check):")
for r in sorted([r for r in rows if r["class"] == "global_flip_or_rotation"],
                key=lambda r: -int(r["n_sv_total"]))[:10]:
    print(f"  {r['pair']}: total={r['n_sv_total']} art={r['n_rearr_artifact']} "
          f"strandOK={r['n_rearr_strand_confirmed']} csVer={r['n_indel_cs_verified']} "
          f"unv={r['n_indel_small_or_unverified']} bnd={r['n_indel_boundary_artifact']} off={r['n_off_chromosome']}")
