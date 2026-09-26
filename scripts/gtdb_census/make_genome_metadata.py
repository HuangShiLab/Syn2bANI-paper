#!/usr/bin/env python3
"""Build the census genome metadata table, retaining each genome's species
cluster AND its representative genome.

Sources (on the HPC):
  bac120/ar53_metadata_r207.tsv   GTDB metadata: gtdb_genome_representative,
                                  gtdb_taxonomy
  genomes_all/                    present representative FASTAs
  download_missing/urls_missing_genomes.tsv   NCBI URLs for member genomes

Output: <workdir>/genome_metadata.tsv
  accession  species_cluster  representative_accession  is_representative
  present_locally  ncbi_url  gtdb_taxonomy
"""
import argparse
import os
from pathlib import Path


def strip_prefix(acc):
    return acc[3:] if acc.startswith(("GB_", "RS_")) else acc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bac120", required=True)
    ap.add_argument("--ar53", required=True)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--urls", required=True)
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()

    urls = {}
    for line in open(a.urls):
        u, acc = line.rstrip("\n").split("\t")
        urls[acc] = u

    rows = []
    for meta in (a.bac120, a.ar53):
        with open(meta) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {n: i for i, n in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                acc = strip_prefix(f[col["accession"]])
                rep = strip_prefix(f[col["gtdb_genome_representative"]])
                tax = f[col["gtdb_taxonomy"]]
                sp = next((p[3:] for p in tax.split(";") if p.startswith("s__")), "")
                rows.append((acc, sp, rep, acc == rep, tax))

    out = Path(a.workdir) / "genome_metadata.tsv"
    with open(out, "w") as fh:
        fh.write("accession\tspecies_cluster\trepresentative_accession\t"
                 "is_representative\tpresent_locally\tncbi_url\tgtdb_taxonomy\n")
        n_loc = 0
        for acc, sp, rep, is_rep, tax in rows:
            present = os.path.exists(os.path.join(a.genomes_dir, f"{acc}.fna"))
            n_loc += present
            fh.write(f"{acc}\t{sp}\t{rep}\t{int(is_rep)}\t{int(present)}\t"
                     f"{urls.get(acc, '')}\t{tax}\n")
    print(f"{len(rows):,} genomes; {n_loc:,} present locally; reps: "
          f"{sum(1 for r in rows if r[3]):,}; -> {out}")


if __name__ == "__main__":
    main()
