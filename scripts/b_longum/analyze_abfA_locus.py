#!/usr/bin/env python3
"""analyze_abfA_locus.py — validate Syn2bANI struct calls at the abfA locus
against curated metadata for 185 B. longum genomes.

Reference: FSHHK16M1_ctg (abfA+); abfA cluster at contig10:8,546-37,075
(28.5 kb), mapped by blastn of gene_list_abfa.fasta (>=95% id, >=200 bp,
1,587/1,601 genes on contig10).

Inputs (under B_LONGUM_WORK, default results/b_longum_abfA):
  struct_vs_ref/<acc>.struct.tsv   syn2bani struct <acc> vs FSHHK16M1
  metadata.tsv                     accession, abfA_status, hypba_status, ...

Output:
  abfA_locus_validation.tsv   per-genome struct call vs metadata
  ABFA_LOCUS_REPORT.md        summary stats
"""
import glob
import os
from pathlib import Path

import pandas as pd

WORK = Path(os.environ.get("B_LONGUM_WORK", "results/b_longum_abfA"))
LOCUS_CONTIG = "FSHHK16M1_ctg_contig10"
LOCUS_S, LOCUS_E = 8_546, 37_075
LOCUS_LEN = LOCUS_E - LOCUS_S
# Fraction of the locus that must be called deleted to call the cluster absent.
CALL_FRAC = 0.5


def locus_deleted_bp(tsv):
    """Total Deletion size overlapping the abfA locus on the reference."""
    deleted = 0
    if not os.path.exists(tsv):
        return None
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[2] != "Deletion" or f[6] != LOCUS_CONTIG:
                continue
            r_start, r_end, size = int(f[7]), int(f[8]), int(f[9])
            ov = max(0, min(r_end, LOCUS_E) - max(r_start, LOCUS_S))
            if ov > 0:
                deleted += min(ov, size)
    return deleted


def main():
    meta = pd.read_csv(WORK / "metadata.tsv", sep="\t").set_index("accession")
    rows = []
    for tsv in sorted(glob.glob(str(WORK / "struct_vs_ref" / "*.struct.tsv"))):
        acc = Path(tsv).name[: -len(".struct.tsv")]
        if acc not in meta.index:
            continue
        deleted = locus_deleted_bp(tsv)
        if deleted is None:
            continue
        call = "deleted" if deleted >= CALL_FRAC * LOCUS_LEN else "complete"
        rows.append({
            "accession": acc,
            "struct_deleted_bp": deleted,
            "struct_call": call,
            "abfA_status": meta.loc[acc, "abfA_status"],
            "hypba_status": meta.loc[acc, "hypba_status"],
        })
    df = pd.DataFrame(rows)
    df.to_csv(WORK / "abfA_locus_validation.tsv", sep="\t", index=False)

    ct = pd.crosstab(df["struct_call"], df["abfA_status"])
    tp = ct.get("complete", pd.Series(dtype=int)).get("deleted", 0)
    fn = ct.get("complete", pd.Series(dtype=int)).get("complete", 0)
    fp = ct.get("deleted", pd.Series(dtype=int)).get("complete", 0)
    tn = ct.get("deleted", pd.Series(dtype=int)).get("deleted", 0)
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")

    with open(WORK / "ABFA_LOCUS_REPORT.md", "w") as out:
        out.write("# abfA locus validation: Syn2bANI struct vs curated metadata\n\n")
        out.write(f"Genomes compared: {len(df)} (reference FSHHK16M1_ctg excluded)\n\n")
        out.write("## Confusion (struct call rows, metadata columns)\n\n")
        out.write(ct.to_string() + "\n\n")
        out.write(f"Sensitivity (detect curated abfA-deleted): {sens:.3f}\n\n")
        out.write(f"Specificity (curated abfA-complete not called deleted): {spec:.3f}\n")
    print(ct)
    print(f"sensitivity {sens:.3f}  specificity {spec:.3f}  (n={len(df)})")


if __name__ == "__main__":
    main()
