#!/usr/bin/env python3
"""Run syn2bani struct on the top 50 breakpoint-count pairs."""

import csv
import subprocess
from pathlib import Path

WORK_DIR = Path("/Volumes/MoneyCat/Data/fda_argos_staphylococcus_aureus")
GENOMES_DIR = WORK_DIR / "genomes"
SV_DIR = WORK_DIR / "sv"
SV_DIR.mkdir(parents=True, exist_ok=True)

SYN2BIN = Path("/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani")


def parse_fasta_headers(fasta_path):
    seq_ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                seq_ids.append(seq_id)
    return seq_ids


def build_seq_to_fasta_map(genomes_dir):
    mapping = {}
    for fna in sorted(genomes_dir.glob("*.fna")):
        for seq_id in parse_fasta_headers(fna):
            mapping[seq_id] = fna
    return mapping


def main():
    triangle_path = WORK_DIR / "triangle.tsv"
    seq_map = build_seq_to_fasta_map(GENOMES_DIR)

    pairs = []
    with open(triangle_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            q = row["query"]
            r = row["reference"]
            bp = int(row["breakpoint_count"])
            ani = float(row["ani"])
            q_fna = seq_map.get(q)
            r_fna = seq_map.get(r)
            if q_fna and r_fna:
                pairs.append((bp, ani, q, r, q_fna, r_fna))

    pairs.sort(key=lambda x: x[0], reverse=True)
    top50 = pairs[:50]

    fail_log = WORK_DIR / "struct_failures.log"
    failures = []
    for i, (bp, ani, q, r, q_fna, r_fna) in enumerate(top50, 1):
        stem = f"{i:02d}_{q}__{r}"
        # Sanitize filename
        stem = "".join(c if c.isalnum() or c in "_-_." else "_" for c in stem)
        out_bed = SV_DIR / f"{stem}.bed"
        print(f"[{i}/50] {q} vs {r} (breakpoints={bp}, ANI={ani})")
        cmd = [
            str(SYN2BIN),
            "struct",
            "--bed",
            str(q_fna),
            str(r_fna),
        ]
        try:
            with open(out_bed, "w") as out_fh:
                subprocess.run(cmd, check=True, stdout=out_fh, stderr=subprocess.PIPE, text=True, timeout=600)
        except subprocess.CalledProcessError as e:
            failures.append((q, r, e.stderr or e.stdout or "unknown error"))
            print(f"  FAILED: {e.stderr or e.stdout}")
        except subprocess.TimeoutExpired:
            failures.append((q, r, "timeout"))
            print("  FAILED: timeout")

    with open(fail_log, "w") as fh:
        for q, r, msg in failures:
            fh.write(f"{q}\t{r}\t{msg}\n")

    print(f"\nWrote {50 - len(failures)}/50 BED files to {SV_DIR}")
    if failures:
        print(f"Failures logged to {fail_log}")


if __name__ == "__main__":
    main()
