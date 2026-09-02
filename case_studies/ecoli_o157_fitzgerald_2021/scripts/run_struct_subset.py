#!/usr/bin/env python3
"""Run syn2bani struct on the 50 pairs with the highest breakpoint_count."""
import csv
import subprocess
import re
import sys
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
GENOMES_DIR = WORKDIR / "genomes"
SV_DIR = WORKDIR / "sv"
SV_DIR.mkdir(exist_ok=True)
BIN = "/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani"


def strip_acc(s):
    """NZ_CP008957.1 -> CP008957, NC_002695.2 -> NC_002695"""
    s = re.sub(r"^NZ_", "", s)
    return s.split(".")[0]


def header_to_base(fasta_path):
    with fasta_path.open() as fh:
        first = fh.readline().strip()
    if not first.startswith(">"):
        return None
    return strip_acc(first[1:].split()[0])


def main():
    # Map base nucleotide accession -> genome file
    acc_to_file = {}
    for fna in sorted(GENOMES_DIR.glob("*.fna")):
        base = header_to_base(fna)
        if base:
            acc_to_file[base] = fna

    # Load triangle.tsv
    pairs = []
    with (WORKDIR / "triangle.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            q = strip_acc(row["query"])
            r = strip_acc(row["reference"])
            if q == r:
                continue
            pairs.append({
                "q": q, "r": r,
                "breakpoint_count": int(row["breakpoint_count"]),
                "ani": float(row["ani"]),
            })

    top = sorted(pairs, key=lambda x: x["breakpoint_count"], reverse=True)[:50]

    failures = []
    for item in top:
        q_acc, r_acc = item["q"], item["r"]
        q_path = acc_to_file.get(q_acc)
        r_path = acc_to_file.get(r_acc)
        if not q_path or not r_path:
            failures.append((q_acc, r_acc, "missing genome file mapping"))
            continue
        out = SV_DIR / f"{q_acc}__vs__{r_acc}.bed"
        if out.exists():
            print(f"Skip existing {out.name}", file=sys.stderr)
            continue
        cmd = [BIN, "struct", "--bed", str(q_path), str(r_path)]
        print(f"Running {q_acc} vs {r_acc} (bp={item['breakpoint_count']}, ani={item['ani']})", file=sys.stderr)
        try:
            with out.open("w") as fh:
                subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            failures.append((q_acc, r_acc, e.stderr.decode()[:200]))
            if out.exists():
                out.unlink()

    if failures:
        with (WORKDIR / "struct_failures.log").open("w") as fh:
            for q, r, err in failures:
                fh.write(f"{q}\t{r}\t{err}\n")
        print(f"struct finished. Failures: {len(failures)}", file=sys.stderr)
    else:
        print("struct finished for all 50 pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
