#!/usr/bin/env python3
"""Cross ANI-identical FDA-ARGOS E. coli pairs with AMR genotype differences.

Reads pairs99.9_ecoli.tsv (skani 0.3.2 triangle, pairs >= 99.9 ANI) and
metadata_ecoli.tsv (AMRFinderPlus genotypes via NCBI Pathogen Detection).
Writes amr_discordance.tsv: one row per pair with both-sided AMR data,
listing genes unique to each side. Prints summary counts by ANI band.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def amr_set(acc, meta):
    r = meta.get(acc)
    if r is None:
        return None
    s = (r.get("AMR_genotypes_core") or "").strip()  # core only: acquired genes + resistance point mutations
    if s in ("", "-", "none", "None"):
        return frozenset()
    # normalise: strip completeness suffix, dedupe
    return frozenset(x.strip().split("=")[0] for x in s.replace(";", ",").split(",") if x.strip())


def main():
    meta = {}
    with open(os.path.join(HERE, "metadata_ecoli.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["assembly_acc"]] = row

    # gene prevalence across the collection: prevalent genes are intrinsic
    # (chromosomal housekeeping resistance); rare genes are acquired elements
    prevalence = {}
    n_genomes = 0
    for acc in meta:
        g = amr_set(acc, meta)
        if g is None:
            continue
        n_genomes += 1
        for gene in g:
            prevalence[gene] = prevalence.get(gene, 0) + 1
    intrinsic = {g for g, n in prevalence.items() if n >= 0.5 * n_genomes}

    pairs = []
    with open(os.path.join(HERE, "pairs99.9_ecoli.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pairs.append((row["acc1"], row["strain1"], row["acc2"], row["strain2"], float(row["ani"])))

    out_rows = []
    bands = [99.9, 99.99, 99.999]
    counts = {b: [0, 0, 0] for b in bands}  # band -> [n_both_amr, n_diff_all, n_diff_acquired]
    for a1, s1, a2, s2, ani in pairs:
        g1, g2 = amr_set(a1, meta), amr_set(a2, meta)
        if g1 is None or g2 is None:
            continue
        sym = g1 ^ g2
        diff = bool(sym)
        diff_acq = bool(sym - intrinsic)
        for b in bands:
            if ani >= b:
                counts[b][0] += 1
                if diff:
                    counts[b][1] += 1
                if diff_acq:
                    counts[b][2] += 1
        if diff:
            acq1 = sorted((g1 - g2) - intrinsic)
            acq2 = sorted((g2 - g1) - intrinsic)
            out_rows.append((f"{ani:.4f}", s1, s2, a1, a2,
                             ",".join(sorted(g1 - g2)), ",".join(sorted(g2 - g1)),
                             ",".join(acq1), ",".join(acq2)))

    out_rows.sort(key=lambda r: -float(r[0]))
    out_path = os.path.join(HERE, "amr_discordance.tsv")
    with open(out_path, "w") as fh:
        fh.write("ani\tstrain1\tstrain2\tacc1\tacc2\tonly_in_1\tonly_in_2\tacquired_only_in_1\tacquired_only_in_2\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")

    print(f"intrinsic genes (>=50% prevalence): {sorted(intrinsic)}")
    print(f"total pairs >=99.9: {len(pairs)}")
    for b in bands:
        n, d, da = counts[b]
        print(f"ANI >= {b}: {n} pairs with AMR data | any-genotype diff: {d} ({100*d/max(n,1):.1f}%) | acquired-gene diff: {da} ({100*da/max(n,1):.1f}%)")
    print(f"wrote {len(out_rows)} discordant rows -> {out_path}")


if __name__ == "__main__":
    main()
