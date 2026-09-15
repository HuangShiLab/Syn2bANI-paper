#!/usr/bin/env python3
"""Re-run a case-study collection all-vs-all with the current syn2bani binary
and summarise breakpoint_count against ANI.

The original case-study tables (results/triangle.tsv) were produced before
Syn2bANI c974f5f (2026-09-01) and ce60243 (v0.1.1), whose breakpoint_count
counted contig ends, repeat chains and collinear chain breaks as breakpoints
and reported hundreds of "breakpoints" for identical genomes. This script
regenerates the pairwise table with a current binary so the manuscript quotes
the fixed metric.

`--dedup-stem` keeps one genome per assembly accession stem: an accession list
that carries both the GenBank and the RefSeq copy of an assembly (or two
versions of one) otherwise compares a genome with itself under two names. 54
of the 122 FDA-ARGOS *S. aureus* accessions are such twins, so the honest
collection is 67 genomes and 2,211 pairs, not 122 and 7,381.

Usage:
    python3 rerun_breakpoints.py --study ecoli_o157_fitzgerald_2021 \
        --metadata results/metadata_with_lineage.tsv --group-cols assigned_lineage,host_category
    python3 rerun_breakpoints.py --study fda_argos_s_aureus \
        --metadata results/assembly_metadata.tsv --group-cols country

Expects <study>/genomes/<assembly_acc>.fna (see fetch_ncbi_fna.py). Writes
<study>/results/triangle_rerun.tsv, <study>/results/rerun_summary.md and
<study>/figures/ani_vs_breakpoints_rerun.png.
"""
import argparse
import glob
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def dedup_stem(paths):
    """One path per accession stem: prefer RefSeq (GCF), then the highest
    version. GCF_001019195.1 and GCF_001019195.2 are two releases of one
    assembly; GCA_000626615.2 and GCF_000626615.1 are two copies of one."""
    best = {}
    for p in paths:
        acc = os.path.basename(p)[:-4]
        stem = acc.split("_", 1)[1].split(".")[0]
        try:
            version = int(acc.split(".")[-1])
        except ValueError:
            version = 0
        key = (acc.startswith("GCF"), version)
        if stem not in best or key > best[stem][0]:
            best[stem] = (key, p)
    return sorted(p for _, p in best.values())


def seq_to_assembly(genomes_dir):
    m = {}
    for fna in glob.glob(os.path.join(genomes_dir, "*.fna")):
        acc = os.path.basename(fna)[:-4]
        with open(fna) as fh:
            for line in fh:
                if line.startswith(">"):
                    m[line[1:].split()[0]] = acc
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--metadata", required=True, help="TSV with an assembly_acc column, relative to the study dir")
    ap.add_argument("--group-cols", default="", help="comma-separated metadata columns to summarise by")
    ap.add_argument("--syn2bani", default=os.path.join(HERE, "..", "..", "Syn2bANI", "target", "release", "syn2bani"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dedup-stem", action="store_true",
                    help="keep one genome per assembly accession stem (GCA/GCF twins, old versions)")
    args = ap.parse_args()

    study = os.path.join(HERE, args.study)
    genomes = sorted(glob.glob(os.path.join(study, "genomes", "*.fna")))
    if not genomes:
        sys.exit(f"no genomes under {study}/genomes")
    if args.dedup_stem:
        kept = dedup_stem(genomes)
        print(f"--dedup-stem: {len(genomes)} files -> {len(kept)} distinct assemblies", file=sys.stderr)
        genomes = kept
    results = os.path.join(study, "results")
    figures = os.path.join(study, "figures")
    os.makedirs(results, exist_ok=True)
    os.makedirs(figures, exist_ok=True)
    out = os.path.join(results, "triangle_rerun.tsv")

    version = subprocess.run([args.syn2bani, "--version"], capture_output=True, text=True).stdout.strip()
    git = subprocess.run(["git", "-C", os.path.dirname(os.path.dirname(os.path.dirname(args.syn2bani))),
                          "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    print(f"{version} @ {git}; {len(genomes)} genomes", file=sys.stderr)
    cmd = [args.syn2bani, "triangle", "--edge-list", "--verbose", "-t", str(args.threads), "-o", out] + genomes
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    tri = pd.read_csv(out, sep="\t")
    m = seq_to_assembly(os.path.join(study, "genomes"))
    tri["q_acc"] = tri["query"].map(m)
    tri["r_acc"] = tri["reference"].map(m)
    tri = tri[tri.q_acc.notna() & tri.r_acc.notna() & (tri.q_acc != tri.r_acc)]
    tri["ani_pct"] = tri["ani_gated"] if "ani_gated" in tri else tri["ani"]
    meta = pd.read_csv(os.path.join(study, args.metadata), sep="\t", dtype=str)
    meta = meta.drop_duplicates("assembly_acc").set_index("assembly_acc")
    group_cols = [c for c in args.group_cols.split(",") if c]
    for c in group_cols:
        tri[f"q_{c}"] = tri.q_acc.map(meta[c])
        tri[f"r_{c}"] = tri.r_acc.map(meta[c])
    tri.to_csv(out, sep="\t", index=False)

    bp = tri.breakpoint_count.astype(int)
    rho = stats.spearmanr(tri.ani_pct, bp)
    lines = [f"# {args.study}: re-run with {version} ({git})", "",
             *([f"- accession stems deduplicated: one genome per assembly", ""] if args.dedup_stem else []),
             f"- genomes: {len(genomes)}; non-self pairs: {len(tri):,}",
             f"- ANI range: {tri.ani_pct.min():.4f} – {tri.ani_pct.max():.4f}",
             f"- breakpoint_count: median {bp.median():.0f}, IQR {bp.quantile(.25):.0f}–{bp.quantile(.75):.0f}, "
             f"max {bp.max()}, pairs at 0: {(bp == 0).mean():.1%}",
             f"- synteny_blocks: median {tri.synteny_blocks.median():.0f}",
             f"- Spearman(ANI, breakpoint_count) = {rho.statistic:.3f} (p = {rho.pvalue:.2g})",
             f"- pairs with ANI ≥ 99.9% and breakpoint_count ≥ 10: {((tri.ani_pct >= 99.9) & (bp >= 10)).sum():,}",
             ""]
    for c in group_cols:
        same = tri[tri[f"q_{c}"] == tri[f"r_{c}"]]
        diff = tri[tri[f"q_{c}"] != tri[f"r_{c}"]]
        lines += [f"## breakpoint_count by {c}", "",
                  f"- same-{c} pairs: n = {len(same):,}, median {same.breakpoint_count.median():.0f}",
                  f"- different-{c} pairs: n = {len(diff):,}, median {diff.breakpoint_count.median():.0f}",
                  "", f"| {c} (both members) | pairs | median | IQR | max |", "|---|---:|---:|---:|---:|"]
        for g, gg in same.groupby(f"q_{c}"):
            b = gg.breakpoint_count.astype(int)
            lines.append(f"| {g} | {len(gg)} | {b.median():.0f} | {b.quantile(.25):.0f}–{b.quantile(.75):.0f} | {b.max()} |")
        lines.append("")
    with open(os.path.join(results, "rerun_summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # (a) ANI vs breakpoints, coloured by the first group column when the two
    # members agree; (b) the same by the second group column, or the breakpoint
    # distribution when only one column was given.
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
    def scatter_by(ax, col):
        same = tri[f"q_{col}"] == tri[f"r_{col}"]
        ax.scatter(tri.ani_pct[~same], bp[~same], s=6, alpha=0.25, color="0.7",
                   label=f"different {col.replace('_', ' ')}")
        for g, gg in tri[same].groupby(f"q_{col}"):
            ax.scatter(gg.ani_pct, gg.breakpoint_count, s=8, alpha=0.7, label=str(g))
        ax.legend(fontsize=6, frameon=False, loc="upper left")
        ax.set_xlabel("Syn2bANI ANI (%)")
        ax.set_ylabel("breakpoint_count")
        ax.set_title(f"by {col.replace('_', ' ')}", fontsize=9)

    if group_cols:
        scatter_by(axes[0], group_cols[0])
    else:
        axes[0].scatter(tri.ani_pct, bp, s=6, alpha=0.35, color="steelblue")
        axes[0].set_xlabel("Syn2bANI ANI (%)")
        axes[0].set_ylabel("breakpoint_count")
    if len(group_cols) > 1:
        scatter_by(axes[1], group_cols[1])
    else:
        axes[1].hist(bp, bins=range(0, int(bp.max()) + 2), color="steelblue")
        axes[1].set_xlabel("breakpoint_count")
        axes[1].set_ylabel("pairs")
        axes[1].set_title(f"median {bp.median():.0f}, {(bp == 0).mean():.0%} at zero", fontsize=9)
    fig.suptitle(f"{args.study}: {len(genomes)} genomes, {len(tri):,} pairs, "
                 f"{version} ({git})", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(figures, "ani_vs_breakpoints_rerun.png"), dpi=300)
    fig.savefig(os.path.join(figures, "ani_vs_breakpoints_rerun.pdf"))


if __name__ == "__main__":
    main()
