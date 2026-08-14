#!/usr/bin/env python3
"""Measure the chain-coverage ascertainment function pi(v): the probability
that a position in a region with per-site divergence v ends up inside a chain
(span-covered) and that its tags enter the likelihood.

Also generates extra mosaic replicates (different seeds and block sizes) to
test whether pi is a stable, algorithm-determined function or case-specific.
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all
from explore_mosaic import CASES, block_rates, load_paf, coverage_mask, MOSAIC, BLOCK

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = "/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani"
EXTRA = os.path.join(MOSAIC, "extra")


def gen_extra():
    """Extra mosaic replicates: other seeds, other block sizes, other means."""
    os.makedirs(EXTRA, exist_ok=True)
    ref = read_fasta_single(os.path.join(MOSAIC, "ref.fasta"))
    n = len(ref)
    sys.path.insert(0, "/Users/macstudio/Downloads/Syn2bANI/prototype")
    from simulate_mosaic import mutate_blockwise
    from simulate_accessory import write_fasta
    from simulate import BASES
    seq_u8 = np.frombuffer(ref.encode(), dtype=np.uint8).copy()

    jobs = []
    for seed, regime, mean_ani, param, block in [
        (1001, "gamma", 0.90, 1.0, 5000),
        (1002, "gamma", 0.90, 0.5, 5000),
        (1003, "gamma", 0.85, 1.0, 5000),
        (1004, "gamma", 0.92, 1.5, 5000),
        (1005, "bimodal", 0.90, 0.70, 5000),
        (1006, "bimodal", 0.85, 0.60, 5000),
        (1007, "gamma", 0.90, 1.0, 1000),   # fine-grained heterogeneity
        (1008, "gamma", 0.90, 1.0, 20000),  # coarse-grained
        (1009, "bimodal", 0.90, 0.70, 1000),
        (1010, "gamma", 0.95, 0.7, 5000),
    ]:
        name = f"x{seed}_{regime}{param}_ani{mean_ani}_b{block//1000}k"
        jobs.append((name, seed, regime, mean_ani, param, block))
    manifest = []
    for name, seed, regime, mean_ani, param, block in jobs:
        rng = np.random.default_rng(seed)
        n_blocks = (n + block - 1) // block
        rates = block_rates(regime, mean_ani, param, n_blocks, seed)
        mut, n_sub = mutate_blockwise(seq_u8, rates, block, rng)
        true_ani = 1.0 - n_sub / n
        write_fasta(os.path.join(EXTRA, name + ".fasta"), name, mut)
        manifest.append((name, true_ani, regime, param, block))
        if not os.path.exists(os.path.join(EXTRA, name + ".ani.tsv")):
            subprocess.run([BIN, "ani", os.path.join(EXTRA, name + ".fasta"),
                            os.path.join(MOSAIC, "ref.fasta"), "--verbose",
                            "--strata-out", os.path.join(EXTRA, name + ".strata.tsv"),
                            "-o", os.path.join(EXTRA, name + ".ani.tsv")],
                           check=True, capture_output=True)
            subprocess.run([BIN, "struct", os.path.join(EXTRA, name + ".fasta"),
                            os.path.join(MOSAIC, "ref.fasta"), "--paf",
                            "-o", os.path.join(EXTRA, name + ".paf")],
                           check=True, capture_output=True)
    with open(os.path.join(EXTRA, "manifest.tsv"), "w") as fh:
        fh.write("name\ttrue_ani\tregime\tparam\tblock\n")
        for row in manifest:
            fh.write("\t".join(map(str, row)) + "\n")
    return manifest


def per_block_coverage(qname, ref_len, spans, block):
    n_blocks = (ref_len + block - 1) // block
    cov = np.zeros(ref_len, dtype=bool)
    for lo, hi in spans:
        cov[lo:hi] = True
    frac = []
    for b in range(n_blocks):
        lo, hi = b * block, min((b + 1) * block, ref_len)
        frac.append(cov[lo:hi].mean())
    return np.array(frac)


def main():
    ref = read_fasta_single(os.path.join(MOSAIC, "ref.fasta"))
    n = len(ref)

    # ---- main 9 cases ----
    print("== main 9 cases: block coverage vs true block divergence ==")
    all_v, all_cov = [], []
    for i, (label, regime, mean_ani, param) in enumerate(CASES):
        rates = block_rates(regime, mean_ani, param, (n + BLOCK - 1) // BLOCK, 90000 + i)
        qname = None
        with open(os.path.join(MOSAIC, "manifest.tsv")) as fh:
            next(fh)
            for line in fh:
                nm = line.split("\t")[0]
                if nm.endswith("__" + label):
                    qname = nm
        spans = load_paf(os.path.join(MOSAIC, qname + ".paf"))
        frac = per_block_coverage(qname, n, spans, BLOCK)
        v = -np.log(1 - rates)
        all_v.append(v)
        all_cov.append(frac)
        # binned coverage curve
        bins = np.array([0, 0.01, 0.02, 0.04, 0.07, 0.10, 0.15, 0.20, 0.30, 0.45, 1.0])
        idx = np.digitize(v, bins) - 1
        line = []
        for b in range(len(bins) - 1):
            m = idx == b
            line.append(f"{frac[m].mean():5.2f}" if m.sum() >= 3 else "    -")
        print(f"{label:<24} cov " + " ".join(line))
    print(f"{'v bin lo':<24}     " + " ".join(f"{b:5.2f}" for b in bins[:-1]))

    v = np.concatenate(all_v)
    c = np.concatenate(all_cov)
    # pooled logistic fit on block-level coverage
    from scipy.optimize import curve_fit

    def logistic(x, v50, s):
        return 1.0 / (1.0 + np.exp((x - v50) / s))

    m = c > -1
    popt, _ = curve_fit(logistic, v[m], c[m], p0=[0.12, 0.03])
    print(f"\npooled logistic pi(v) = 1/(1+exp((v-{popt[0]:.4f})/{popt[1]:.4f}))")

    # ---- extra replicates ----
    if "--gen-extra" in sys.argv:
        manifest = gen_extra()
    else:
        manifest = []
        mp = os.path.join(EXTRA, "manifest.tsv")
        if os.path.exists(mp):
            with open(mp) as fh:
                next(fh)
                for line in fh:
                    nm, tr, reg, prm, blk = line.rstrip("\n").split("\t")
                    manifest.append((nm, float(tr), reg, float(prm), int(blk)))
    if manifest:
        print("\n== extra replicates: block coverage vs divergence ==")
        print(f"{'case':<34} {'truth':>6} {'gamma':>6} {'unif':>6} {'AF':>5}")
        for name, tr, reg, prm, blk in manifest:
            rates = block_rates(reg, {"gamma": None}.get(reg, None) or
                                _mean_ani_from_name(name), prm,
                                (n + blk - 1) // blk, int(name.split("_")[0][1:]))
            spans = load_paf(os.path.join(EXTRA, name + ".paf"))
            frac = per_block_coverage(name, n, spans, blk)
            vv = -np.log(1 - rates)
            m50 = vv < popt[0]
            with open(os.path.join(EXTRA, name + ".ani.tsv")) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                row = dict(zip(header, fh.readline().rstrip("\n").split("\t")))
            print(f"{name:<34} {tr*100:6.2f} {float(row['ani']):6.2f} "
                  f"{float(row['ani_uniform']):6.2f} {float(row['af_query']):5.3f} "
                  f" covLO {frac[m50].mean():.3f} covHI {frac[~m50].mean():.3f}")


def _mean_ani_from_name(name):
    return float(name.split("_ani")[1].split("_")[0])


if __name__ == "__main__":
    main()
