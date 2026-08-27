#!/usr/bin/env python3
"""Analyze ANI–synteny discordance vs abfA cluster status and phenotype in
B. longum.

Expected inputs (user to provide):
  results/b_longum_abfA/metadata.tsv
      Columns: accession, abfA_status (complete/deleted/partial),
               phenotype (effective/ineffective/unknown), [optional] notes
  results/b_longum_abfA/b_longum_s2b_matrix.tsv

Outputs:
  results/b_longum_abfA/ABfA_ANALYSIS_REPORT.md
  results/b_longum_abfA/abfA_pair_metrics.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path(os.environ.get("B_LONGUM_WORK", "results/b_longum_abfA"))


def main():
    meta_path = WORK / "metadata.tsv"
    if not meta_path.exists():
        print(f"ERROR: {meta_path} not found. Create it with columns:")
        print("  accession, abfA_status, phenotype")
        print("abfA_status: complete / deleted / partial")
        print("phenotype: effective / ineffective / unknown")
        return 1

    meta = pd.read_csv(meta_path, sep="\t")
    s2b = pd.read_csv(WORK / "b_longum_s2b_matrix.tsv", sep="\t")

    # annotate pairs
    meta = meta.set_index("accession")
    s2b["q_abfA"] = s2b["query"].map(meta["abfA_status"])
    s2b["r_abfA"] = s2b["reference"].map(meta["abfA_status"])
    s2b["q_pheno"] = s2b["query"].map(meta["phenotype"])
    s2b["r_pheno"] = s2b["reference"].map(meta["phenotype"])

    s2b["same_abfA"] = s2b["q_abfA"] == s2b["r_abfA"]
    s2b["same_pheno"] = s2b["q_pheno"] == s2b["r_pheno"]

    s2b.to_csv(WORK / "abfA_pair_metrics.tsv", sep="\t", index=False)

    L = ["# B. longum abfA ANI–synteny discordance analysis",
         f"Pairs analyzed: {len(s2b)}",
         f"Genomes with metadata: {len(meta)}",
         ""]
    for col in ["abfA_status", "phenotype"]:
        L.append(f"## {col} distribution")
        L.append(meta[col].value_counts().to_string())
        L.append("")

    # Compare ANI and synteny between same vs different abfA/phenotype groups
    for group_col, same_col in [("abfA_status", "same_abfA"), ("phenotype", "same_pheno")]:
        L.append(f"## ANI and synteny by {group_col} concordance")
        for same in [True, False]:
            sub = s2b[s2b[same_col] == same]
            if len(sub) == 0:
                continue
            L.append(f"- same {group_col}: n={len(sub)}, "
                     f"ANI mean={sub['ani'].mean():.3f}, "
                     f"synteny mean={sub['anchor_adjacency'].mean():.4f}, "
                     f"breakpoint median={sub['breakpoint_count'].median():.0f}")
        L.append("")

    # Identify high-ANI low-synteny pairs with different abfA status
    disc = s2b[(s2b["ani"] >= 98.0) & (s2b["anchor_adjacency"] < 0.98) &
               (s2b["same_abfA"] == False)].copy()
    disc = disc.sort_values("anchor_adjacency")
    L.append(f"## High-ANI (>=98%) low-synteny (<0.98) pairs with different abfA status")
    L.append(f"n = {len(disc)}")
    if len(disc):
        L.append(disc[["query", "reference", "ani", "anchor_adjacency",
                       "breakpoint_count", "q_abfA", "r_abfA", "q_pheno", "r_pheno"]]
                 .head(20).to_string(index=False))
    L.append("")

    with open(WORK / "ABfA_ANALYSIS_REPORT.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
