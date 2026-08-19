#!/usr/bin/env python3
"""fasta_index.py genomes.tsv genome_lengths.tsv firstid_map.tsv
genomes.tsv: genome_id \t fasta_path
Outputs per-genome total length and first fasta record id (token up to whitespace).
"""
import sys

genomes_tsv, lengths_out, firstid_out = sys.argv[1:4]

def scan(path):
    total = 0
    first_id = None
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
        fh = opener(path, "rt")
    else:
        fh = opener(path)
    with fh:
        for line in fh:
            if line.startswith(">"):
                if first_id is None:
                    first_id = line[1:].split()[0]
            else:
                total += len(line.strip())
    return total, (first_id or "")

with open(genomes_tsv) as f, open(lengths_out, "w") as lo, open(firstid_out, "w") as fo:
    lo.write("genome_id\tlength_bp\n")
    fo.write("genome_id\tfirst_record_id\n")
    first = True
    for line in f:
        if first and line.startswith("genome_id\t"):
            first = False
            continue
        first = False
        gid, path = line.rstrip("\n").split("\t")
        total, fid = scan(path)
        lo.write(f"{gid}\t{total}\n")
        fo.write(f"{gid}\t{fid}\n")
