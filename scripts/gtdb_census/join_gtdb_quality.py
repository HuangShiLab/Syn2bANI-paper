#!/usr/bin/env python3
"""Join GTDB census genome metadata to CheckM/assembly quality fields.

GTDB metadata accessions are prefixed (GB_GCA_..., RS_GCF_...), while the
census genome manifest uses bare accession names.  This script strips those
prefixes only for the join key and writes one quality row per census genome.
"""
import csv
import sys
from pathlib import Path

W = Path("/lustre1/g/aos_shihuang/data/syn2b_census")
M = Path("/lustre1/g/aos_shihuang/data/gtdb-r207/metadata")
OUT = W / "gtdb_genome_quality.tsv"
QUALITY_FIELDS = [
    "checkm_completeness", "checkm_contamination",
    "checkm_strain_heterogeneity", "contig_count", "n50_contigs",
    "longest_contig", "gc_percentage",
]


def bare_acc(x: str) -> str:
    for prefix in ("GB_", "RS_"):
        if x.startswith(prefix):
            return x[len(prefix):]
    return x


def main():
    meta = {}
    with open(W / "genome_metadata.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["accession"]] = row
    print(f"genome_metadata rows={len(meta)}", file=sys.stderr)

    seen = set()
    with open(OUT, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "species_cluster", "gtdb_taxonomy",
                    *QUALITY_FIELDS])
        for mf in sorted(M.glob("*_metadata_r207.tsv")):
            with open(mf, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = bare_acc(row.get("accession", ""))
                    if acc not in meta or acc in seen:
                        continue
                    seen.add(acc)
                    src = meta[acc]
                    w.writerow([acc, src["species_cluster"],
                                src["gtdb_taxonomy"],
                                *[row.get(x, "") for x in QUALITY_FIELDS]])
    print(f"joined rows={len(seen)} missing={len(meta)-len(seen)}",
          file=sys.stderr)
    if len(seen) != len(meta):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
