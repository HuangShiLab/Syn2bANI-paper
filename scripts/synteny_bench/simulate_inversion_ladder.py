#!/usr/bin/env python3
"""Simulated inversion ladder with exact structural truth.

Evolves E. coli MG1655 to fixed ANI levels with counted substitutions, then
applies a counted number of non-overlapping inversions. Each inversion
contributes exactly 2 adjacency breakpoints, so the structural ground truth
is exact by construction: true_breakpoints = 2 * n_inversions.

Outputs:
  results/synteny_bench/sim/q_ani{A}_inv{N}.fasta
  results/synteny_bench/manifest.tsv   (file, ani, n_inv, true_breakpoints)

Usage: python3 simulate_inversion_ladder.py
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
REF = Path("/Users/macstudio/Downloads/Syn2bANI/prototype/mg1655.fasta")
OUT = ROOT / "results/synteny_bench/sim"
OUT.mkdir(parents=True, exist_ok=True)

ANI_LEVELS = [0.95, 0.98]
INV_COUNTS = [0, 1, 2, 4, 8, 16, 32]
INV_MIN, INV_MAX = 100_000, 400_000
SEED = 7

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_fasta_single(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    assert len(seqs) == 1, "expected single-contig reference"
    return seqs[0]


def mutate(seq_u8, ani, rng):
    n_subs = round((1 - ani) * seq_u8.size)
    pos = rng.choice(seq_u8.size, size=n_subs, replace=False)
    alt = rng.integers(0, 3, size=n_subs)
    new = seq_u8[pos].copy()
    mask = new == 65  # A
    new[mask] = BASES[1:][alt[mask]]
    mask = seq_u8[pos] == 67  # C
    b = BASES[np.array([0, 2, 3])]
    new[mask] = b[alt[mask]]
    mask = seq_u8[pos] == 71  # G
    b = BASES[np.array([0, 1, 3])]
    new[mask] = b[alt[mask]]
    mask = seq_u8[pos] == 84  # T
    new[mask] = BASES[:3][alt[mask]]
    out = seq_u8.copy()
    out[pos] = new
    return out


def apply_inversions(seq_u8, n_inv, rng):
    """n_inv non-overlapping inversions: partition the genome into n_inv
    equal slots and invert one random interval per slot (size 100-400 kb,
    capped at 90% of the slot so dense ladders stay non-overlapping)."""
    L = seq_u8.size
    if n_inv == 0:
        return seq_u8
    out = seq_u8
    slot = L // n_inv
    for i in range(n_inv):
        lo, hi = i * slot, (i + 1) * slot
        max_size = min(INV_MAX, int((hi - lo) * 0.9))
        size = int(rng.integers(min(INV_MIN, max_size), max_size + 1))
        start = int(rng.integers(lo, hi - size))
        end = start + size
        seg = out[start:end][::-1]
        # reverse complement
        comp = np.empty_like(seg)
        comp[seg == 65] = 84
        comp[seg == 67] = 71
        comp[seg == 71] = 67
        comp[seg == 84] = 65
        out = np.concatenate([out[:start], comp, out[end:]])
    return out


def write_fasta(path, name, seq_u8):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        s = seq_u8.tobytes().decode()
        for i in range(0, len(s), 80):
            fh.write(s[i:i + 80] + "\n")


def main():
    ref = read_fasta_single(REF)
    ref_u8 = np.frombuffer(ref.encode(), dtype=np.uint8)
    rows = []
    for ani in ANI_LEVELS:
        rng = np.random.default_rng(SEED + int(ani * 10000))
        evolved = mutate(ref_u8, ani, rng)
        for n_inv in INV_COUNTS:
            r2 = np.random.default_rng(SEED * 100 + n_inv)
            sim = apply_inversions(evolved, n_inv, r2)
            name = f"q_ani{ani:.2f}_inv{n_inv}"
            write_fasta(OUT / f"{name}.fasta", name, sim)
            rows.append((f"{name}.fasta", ani, n_inv, 2 * n_inv))
            print(name, "ok")
    with open(ROOT / "results/synteny_bench/manifest.tsv", "w") as fh:
        fh.write("file\tani\tn_inv\ttrue_breakpoints\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
