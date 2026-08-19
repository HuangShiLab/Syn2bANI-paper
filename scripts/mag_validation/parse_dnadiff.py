#!/usr/bin/env python3
"""parse_dnadiff.py <binid> <dd.report> -> one TSV row on stdout:
bin \t anim_ani \t anim_af_ref \t anim_af_qry(mag) \t aligned_bases_qry \t af_tier
dnadiff report layout (verified): [Bases] has a single 'AlignedBases' line with
'refbases(refpct%)   qrybases(qrypct%)'; [Alignments] first subgroup is 1-to-1
with an 'AvgIdentity' line (ANIm). af_tier: >=60 strict, 30-60 low-AF, <30 verylow.
"""
import re
import sys

binid, report = sys.argv[1], sys.argv[2]
txt = open(report).read()

m = re.search(r"^AlignedBases\s+\S+\(([\d.]+)%\)\s+(\d+)\(([\d.]+)%\)", txt, re.M)
af_ref = float(m.group(1)) if m else float("nan")
alqry = int(m.group(2)) if m else 0
af_qry = float(m.group(3)) if m else float("nan")

# 1-to-1 AvgIdentity: first AvgIdentity after the '1-to-1' line in [Alignments]
anim = float("nan")
sec = re.search(r"\[Alignments\](.*?)^1-to-1\s+\d+\s+\d+\s*$(.*?)^M-to-M", txt, re.M | re.S)
if sec:
    m2 = re.search(r"^AvgIdentity\s+([\d.]+)", sec.group(2), re.M)
    if m2:
        anim = float(m2.group(1))

tier = "strict" if af_qry >= 60 else ("low-AF" if af_qry >= 30 else "verylow-AF")
print(f"{binid}\t{anim:.4f}\t{af_ref:.2f}\t{af_qry:.2f}\t{alqry}\t{tier}")
