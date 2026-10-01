#!/usr/bin/env python3
"""Restrict GTDB census scope to genomes that actually produced TGTs."""
import argparse, csv, json
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    root = Path(args.workdir)
    manifest = json.load(open(root / "manifest.json"))
    outdir = root / "effective_cluster_accessions"
    outdir.mkdir(exist_ok=True)
    audit = root / "genome_inclusion.tsv"
    included_by_cluster = {}
    counts = Counter()
    with open(audit, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "cluster", "status", "fasta_path",
                    "fasta_size", "tgt_path", "tgt_size"])
        for list_path in sorted((root / "cluster_accessions").glob("*.txt")):
            cluster = list_path.name[:-4]
            included = []
            for acc in list_path.read_text().split():
                fasta = manifest.get(acc, "")
                fasta_path = Path(fasta) if fasta else None
                fasta_ok = bool(fasta_path and fasta_path.is_file() and
                                fasta_path.stat().st_size > 0)
                tgt = root / "tgt" / f"{acc}.tgt"
                tgt_ok = tgt.is_file() and tgt.stat().st_size > 0
                if tgt_ok and fasta_ok:
                    status = "included"
                    included.append(acc)
                elif not fasta_ok:
                    status = "missing_or_empty_fasta"
                elif not tgt_ok:
                    status = "missing_or_empty_tgt"
                else:
                    status = "excluded"
                counts[status] += 1
                w.writerow([acc, cluster, status, fasta,
                            fasta_path.stat().st_size if fasta_path and fasta_path.is_file() else "",
                            tgt, tgt.stat().st_size if tgt.is_file() else ""])
            included_by_cluster[cluster] = included
            (outdir / f"{cluster}.txt").write_text(
                "".join(x + "\n" for x in included))
    expected = sum(len(x) * (len(x) - 1) // 2 for x in included_by_cluster.values())
    with open(root / "effective_scope_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in counts.items():
            fh.write(f"genome_{k}\t{v}\n")
        fh.write(f"effective_clusters\t{sum(bool(x) for x in included_by_cluster.values())}\n")
        fh.write(f"effective_genomes\t{sum(len(x) for x in included_by_cluster.values())}\n")
        fh.write(f"effective_unique_pairs\t{expected}\n")
    print(f"included={counts['included']} expected_pairs={expected}")


if __name__ == "__main__":
    main()
