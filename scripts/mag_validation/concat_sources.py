#!/usr/bin/env python3
"""concat_sources.py genomes.tsv source_all.fa
Concatenate source genomes into one multi-fasta with headers rewritten to
>{genome_id}|{orig} so minimap2 PAF target names carry the genome id.
"""
import sys

genomes_tsv, out_fa = sys.argv[1:3]
with open(genomes_tsv) as f, open(out_fa, "w") as out:
    for line in f:
        if line.startswith("genome_id\t"):
            continue
        gid, path = line.rstrip("\n").split("\t")
        with open(path) as fi:
            for l in fi:
                if l.startswith(">"):
                    out.write(f">{gid}|{l[1:]}")
                else:
                    out.write(l)
