#!/usr/bin/env python3
"""In-silico ladder for `breakpoint_count`: exact event counts on MG1655 edits.

`breakpoint_count` is supposed to be an event count (2 per inversion, 3 per
translocation, 0 for indels and substitutions), not a chain count. Two builds
have failed that: before Syn2bANI c974f5f a fragmented reference added
`n_ref - 1`, and before v0.1.1 every chain-to-chain transition counted, so
repeat chains and collinear chain breaks inflated the number by two orders of
magnitude. This script is the regression check: it edits one genome in place
and asserts the expected count for each edit.

Usage:
    python3 analysis/breakpoint_ladder.py \
        --genome ../Syn2bANI/prototype/mg1655.fasta \
        --syn2bani ../Syn2bANI/target/release/syn2bani
Writes results/breakpoint_ladder.tsv and exits non-zero on any mismatch.
"""
import argparse
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RC = str.maketrans("ACGTacgt", "TGCAtgca")


def load(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def rc(s):
    return s.translate(RC)[::-1]


def mutate(seq, rate, seed):
    r = random.Random(seed)
    l = list(seq)
    for _ in range(int(len(l) * rate)):
        i = r.randrange(len(l))
        l[i] = r.choice([c for c in "ACGT" if c != l[i]])
    return "".join(l)


def write(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
    return path


def cases(g):
    """(name, sequence, expected breakpoint_count)."""
    out = [("identical", g, 0)]
    for kb in (500, 100, 10):
        k = kb * 1000
        out.append((f"inversion_{kb}kb", g[:1_000_000] + rc(g[1_000_000:1_000_000 + k])
                    + g[1_000_000 + k:], 2))
    for kb in (500, 100):
        k = kb * 1000
        seg, rest = g[1_000_000:1_000_000 + k], g[:1_000_000] + g[1_000_000 + k:]
        out.append((f"translocation_{kb}kb", rest[:3_000_000] + seg + rest[3_000_000:], 3))
    out.append(("deletion_10kb", g[:2_000_000] + g[2_010_000:], 0))
    r = random.Random(1)
    ins = "".join(r.choice("ACGT") for _ in range(10_000))
    out.append(("insertion_10kb", g[:2_000_000] + ins + g[2_000_000:], 0))
    s = g
    for st in (500_000, 2_000_000, 3_500_000):
        s = s[:st] + rc(s[st:st + 200_000]) + s[st + 200_000:]
    out.append(("three_inversions_200kb", s, 6))
    for pct in (1, 3):
        out.append((f"substitution_{pct}pct", mutate(g, pct / 100, 7), 0))
    inv = g[:1_000_000] + rc(g[1_000_000:1_100_000]) + g[1_100_000:]
    out.append(("inversion_100kb_on_1pct_subst", mutate(inv, 0.01, 9), 2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default=os.path.join(ROOT, "..", "Syn2bANI", "prototype", "mg1655.fasta"))
    ap.add_argument("--syn2bani", default=os.path.join(ROOT, "..", "Syn2bANI", "target", "release", "syn2bani"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "breakpoint_ladder.tsv"))
    args = ap.parse_args()

    g = load(args.genome)
    print(f"{os.path.basename(args.genome)}: {len(g):,} bp", file=sys.stderr)
    rows, bad = [], 0
    with tempfile.TemporaryDirectory() as td:
        ref = write(os.path.join(td, "ref.fa"), "ref", g)
        for name, seq, expected in cases(g):
            q = write(os.path.join(td, "q.fa"), name, seq)
            f = subprocess.run([args.syn2bani, "ani", "--verbose", q, ref],
                               capture_output=True, text=True, check=True
                               ).stdout.strip().split("\n")[-1].split("\t")
            got, blocks, ani = int(f[9]), int(f[7]), float(f[2])
            ok = got == expected
            bad += not ok
            rows.append((name, expected, got, blocks, ani, "ok" if ok else "MISMATCH"))
            print(f"{name:32s} expected {expected:2d}  got {got:4d}  blocks {blocks:3d}"
                  f"  ani {ani:.4f}  {'' if ok else 'MISMATCH'}", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("condition\texpected_breakpoints\tbreakpoint_count\tsynteny_blocks\tani\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {args.out}; {len(rows) - bad}/{len(rows)} as expected", file=sys.stderr)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
