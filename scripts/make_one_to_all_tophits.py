#!/usr/bin/env python3
"""Extract top-1 non-self hits per query from one-to-all search outputs.

Reads $BASE/out_one_to_all/{syn2b_search,skani_search}__<QID>.tsv and writes
$BASE/gtdb_one_to_all_tophits.tsv with columns:
    query_id, query_path, syn2b_ref_path, skani_ref_path, syn2b_ani, skani_ani

Self hits (query accession == reference) are excluded so the SV stage compares
against a genuine nearest neighbour, not the query itself.
"""
import os
import sys
import csv

BASE = "/lustre1/g/aos_shihuang/Syn2bANI-paper-bench"
GENOMES = "/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all"
QUERIES = os.path.join(BASE, "gtdb_one_to_all_queries.tsv")
OUT = os.path.join(BASE, "gtdb_one_to_all_tophits.tsv")


def read_queries():
    queries = []
    with open(QUERIES) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                queries.append((parts[1], parts[5]))  # manifest_acc, path
    return queries


def top_syn2b(qid):
    """Top non-self reference from syn2bani search output; returns (acc, ani)."""
    path = os.path.join(BASE, "out_one_to_all", f"syn2b_search__{qid}.tsv")
    best = None
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ref = row["reference"]
            if ref == qid:
                continue
            if row.get("flag") == "BELOW_DETECTION":
                continue
            try:
                ani = float(row["ani_gated"])
            except (ValueError, TypeError):
                continue
            if ani != ani:  # NaN
                continue
            if best is None or ani > best[1]:
                best = (ref, ani)
    return best


def top_skani(qid):
    """Top non-self reference from skani search output; returns (ref_path, ani)."""
    path = os.path.join(BASE, "out_one_to_all", f"skani_search__{qid}.tsv")
    best = None
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ref_file = row["Ref_file"]
            if os.path.basename(ref_file).replace(".fna", "") == qid:
                continue
            try:
                ani = float(row["ANI"])
            except (ValueError, TypeError):
                continue
            if best is None or ani > best[1]:
                best = (ref_file, ani)
    return best


def main():
    queries = read_queries()
    rows = []
    problems = []
    for qid, qpath in queries:
        s2b = top_syn2b(qid)
        sk = top_skani(qid)
        if s2b is None or sk is None:
            problems.append(qid)
            continue
        s2b_ref_path = os.path.join(GENOMES, s2b[0] + ".fna")
        if not os.path.exists(s2b_ref_path):
            problems.append(f"{qid}: missing {s2b_ref_path}")
            continue
        rows.append((qid, qpath, s2b_ref_path, sk[0], f"{s2b[1]:.4f}", f"{sk[1]:.2f}"))

    with open(OUT, "w") as fh:
        fh.write("query_id\tquery_path\tsyn2b_ref\tskani_ref\tsyn2b_ani\tskani_ani\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"wrote {len(rows)} rows to {OUT}")
    for qid, qpath, s2b_ref, sk_ref, s2b_ani, sk_ani in rows:
        same = "SAME" if s2b_ref == sk_ref else "DIFF"
        print(f"  {qid}: syn2b={os.path.basename(s2b_ref)}({s2b_ani}) skani={os.path.basename(sk_ref)}({sk_ani}) [{same}]")
    if problems:
        print("PROBLEMS:", "; ".join(problems))
        sys.exit(1)


if __name__ == "__main__":
    main()
