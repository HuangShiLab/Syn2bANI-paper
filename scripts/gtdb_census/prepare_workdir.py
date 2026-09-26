#!/usr/bin/env python3
"""Prepare the census workdir: cluster accession lists, task counts, manifest.

Run once on the machine that hosts the GTDB R207 genomes:
  python3 prepare_workdir.py --taxonomy ... --genome-dir /lustre/.../gtdb_r207 \
      --tasks tasks.jsonl --workdir /lustre/.../census

Writes:
  workdir/cluster_accessions/<cluster>.txt   sorted accessions per cluster
  workdir/checkpoints/task_counts.json       cluster -> number of tasks (for the
                                              auditor's reclaim logic)
  workdir/manifest.json                      accession -> FASTA path
Missing genomes are logged to workdir/missing_genomes.txt (tasks adapt: a pair
row is simply absent if a TGT is missing; merge_census reports the shortfall).
"""
import argparse
import gzip
import json
from pathlib import Path

ACCESSION_SUFFIXES = ("_genomic.fna.gz", "_genomic.fna", ".fna.gz", ".fa.gz", ".fna", ".fa")


def strip_prefix(acc):
    return acc[3:] if acc.startswith(("GB_", "RS_")) else acc


def load_clusters(taxonomy_tsv):
    clusters = {}
    op = gzip.open if taxonomy_tsv.endswith(".gz") else open
    with op(taxonomy_tsv, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc, tax = strip_prefix(parts[0]), parts[1]
            sp = next((p[3:] for p in tax.split(";") if p.startswith("s__")), "")
            if sp:
                clusters.setdefault(sp, []).append(acc)
    for accs in clusters.values():
        accs.sort()
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--genome-dir", required=True,
                    help="root containing GTDB genome files (searched recursively)")
    ap.add_argument("--flat-dir", action="store_true",
                    help="genome files are {accession}.fna directly under "
                         "--genome-dir; construct manifest paths without "
                         "scanning (fast on large Lustre directories; missing "
                         "files are tolerated and surfaced by the worker)")
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    clusters = load_clusters(args.taxonomy)
    tasks = [json.loads(l) for l in open(args.tasks)]
    needed = {t["cluster"] for t in tasks}
    if not args.flat_dir:
        genomes = {}
        for p in Path(args.genome_dir).rglob("*"):
            if p.name.endswith(ACCESSION_SUFFIXES):
                genomes.setdefault(p.name.split("_genomic")[0].rsplit(".", 1)[0], p)
        by_stem = {}
        for stem, p in genomes.items():
            by_stem.setdefault(stem, p)
    else:
        genomes, by_stem = {}, {}

    root = Path(args.workdir)
    ca = root / "cluster_accessions"
    ck = root / "checkpoints"
    ca.mkdir(parents=True, exist_ok=True)
    ck.mkdir(parents=True, exist_ok=True)

    manifest, missing = {}, []
    counts = {}
    for cl in sorted(needed):
        accs = clusters[cl]
        (ca / (cl.replace("/", "_") + ".txt")).write_text("\n".join(accs) + "\n")
        counts[cl] = sum(1 for t in tasks if t["cluster"] == cl)
        for acc in accs:
            p = genomes.get(acc) or by_stem.get(acc)
            if p is None and args.flat_dir:
                cand = Path(args.genome_dir) / f"{acc}.fna"
                p = cand  # constructively; existence checked at digest time
            if p is None:
                missing.append(acc)
            else:
                manifest[acc] = str(p)

    (ck / "task_counts.json").write_text(
        "\n".join(f"{cl}\t{n}" for cl, n in counts.items()) + "\n")
    (root / "manifest.json").write_text(json.dumps(manifest))
    (root / "missing_genomes.txt").write_text("\n".join(missing) + "\n")
    print(f"clusters: {len(counts):,}  genomes: {len(manifest):,}  "
          f"missing: {len(missing):,}  -> {root/'manifest.json'}")


if __name__ == "__main__":
    main()
