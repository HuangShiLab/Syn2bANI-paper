#!/usr/bin/env python3
"""Find which taxa carry large inversions, from pairs already computed.

The inverted-fraction channel is only interpretable on near-closed assemblies:
each contig is deposited in an arbitrary orientation, so a fragmented pair drifts
to min(f, 1-f) ~ 0.5 with no biology in it. Measured on the >=97% ANIm pairs, the
median runs 0.4726 at >100 contigs down to 0.2880 at 1-2, and the collinear
fraction runs 0.004 up to 0.332.

So the pairs worth looking at are the near-closed ones, and the question is which
species they belong to. This annotates them with GTDB taxonomy and ranks species
by how often their near-closed pairs disagree in orientation -- i.e. where to go
looking for an inversion with a phenotype attached, chosen from data rather than
from a guess.

Contig count is not in the pair tables, so it is recovered from
`observable_fraction ~ 1 - (K-1)/S`, which MATH_REVIEW.md section 7 derives and
which held to 4 dp up to K=300. It is an estimate: the output keeps it as
`K_est` and the closed subset should be confirmed against
`ncbi_assembly_level` / `contig_count` in the metadata, which this script joins.

Usage:
    python3 scripts/annotate_closed_inversion_pairs.py \
        --results  results/gtdb50k \
        --metadata /lustre1/g/aos_shihuang/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
        --out      results/gtdb50k/closed_inversion_pairs.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def read_unique(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df.drop_duplicates(subset="pairid", keep="first")


def load_gtdb(path: Path) -> pd.DataFrame:
    """GTDB metadata, keyed by bare GCF_/GCA_ accession."""
    want = [
        "accession",
        "gtdb_taxonomy",
        "contig_count",
        "checkm_completeness",
        "checkm_contamination",
        "genome_size",
        "ncbi_assembly_level",
    ]
    df = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in want if c not in df.columns]
    if missing:
        sys.exit(f"metadata is missing columns {missing}; got {list(df.columns)[:12]} ...")
    df = df[want].copy()
    # GTDB prefixes accessions with RS_ / GB_
    df["acc"] = df.accession.str.replace(r"^(RS_|GB_)", "", regex=True)
    tax = df.gtdb_taxonomy.str.split(";", expand=True)
    for i, rank in enumerate(TAX_RANKS):
        if i < tax.shape[1]:
            df[rank] = tax[i].str.replace(r"^[a-z]__", "", regex=True)
    return df.drop(columns=["accession", "gtdb_taxonomy"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/gtdb50k")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", default="results/gtdb50k/closed_inversion_pairs.tsv")
    ap.add_argument("--max-contigs", type=float, default=3.0,
                    help="K_est below which a pair counts as near-closed (default 3)")
    ap.add_argument("--min-ani", type=float, default=97.0)
    ap.add_argument("--inverted-cut", type=float, default=0.15,
                    help="min(f,1-f) above which a pair is called inverted (default 0.15)")
    args = ap.parse_args()

    d = Path(args.results)
    m = read_unique(d / "dnadiff_inverted_fraction_high_ani_all.tsv").merge(
        read_unique(d / "syn2b_inverted_fraction_high_ani_all.tsv"), on="pairid"
    )
    m = m.merge(read_unique(d / "high_ani_truth.tsv"), on="pairid", how="left")
    m = m[m.status == "ok"].copy()
    for c in [
        "dnadiff_inverted_fraction",
        "syn2b_inverted_fraction",
        "syn2b_raw_inverted_fraction",
        "syn2b_observable_fraction",
        "syn2b_shared_tags",
        "anim_ani",
        "anim_af_ref",
    ]:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m = m.dropna(subset=["dnadiff_inverted_fraction", "syn2b_observable_fraction",
                         "syn2b_shared_tags"])

    m["K_est"] = 1 + m.syn2b_shared_tags * (1 - m.syn2b_observable_fraction)
    m["dd_major"] = np.minimum(m.dnadiff_inverted_fraction, 1 - m.dnadiff_inverted_fraction)
    m["s2b_major"] = m.syn2b_inverted_fraction

    sel = m[(m.K_est < args.max_contigs) & (m.anim_ani >= args.min_ani)].copy()
    print(f"{len(sel)} near-closed pairs at >= {args.min_ani}% ANIm "
          f"(K_est < {args.max_contigs}) out of {len(m)}", flush=True)

    gtdb = load_gtdb(Path(args.metadata))
    q = gtdb.add_prefix("q_")
    r = gtdb.add_prefix("r_")
    sel = sel.merge(q, on="q_acc", how="left").merge(r, on="r_acc", how="left")

    matched = sel.q_species.notna() & sel.r_species.notna()
    print(f"taxonomy joined for {matched.sum()}/{len(sel)} pairs", flush=True)

    sel["same_species"] = sel.q_species == sel.r_species
    sel["inverted"] = sel.dd_major >= args.inverted_cut
    cols = [
        "pairid", "q_acc", "r_acc", "anim_ani", "anim_af_ref",
        "K_est", "q_contig_count", "r_contig_count",
        "q_ncbi_assembly_level", "r_ncbi_assembly_level",
        "syn2b_shared_tags", "syn2b_observable_fraction",
        "dd_major", "s2b_major", "dnadiff_inverted_fraction",
        "syn2b_raw_inverted_fraction",
        "q_species", "r_species", "q_genus", "q_family", "q_phylum",
        "same_species", "inverted",
    ]
    out = sel[[c for c in cols if c in sel.columns]].sort_values("dd_major", ascending=False)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}", flush=True)

    print("\n=== species ranked by near-closed inverted pairs ===")
    g = (
        sel[sel.same_species]
        .groupby("q_species")
        .agg(
            n_pairs=("dd_major", "size"),
            n_inverted=("inverted", "sum"),
            median_major=("dd_major", "median"),
            max_major=("dd_major", "max"),
            median_ani=("anim_ani", "median"),
        )
    )
    g["frac_inverted"] = (g.n_inverted / g.n_pairs).round(3)
    g = g[g.n_pairs >= 3].sort_values(["n_inverted", "median_major"], ascending=False)
    print(g.head(30).round(4).to_string() if len(g) else "  (no species with >=3 pairs)")

    print("\n=== the 20 most strongly inverted near-closed pairs ===")
    top = out.head(20)
    show = [c for c in ["pairid", "anim_ani", "dd_major", "s2b_major",
                        "q_contig_count", "r_contig_count", "q_species"] if c in top]
    print(top[show].to_string(index=False))
    print("\nAgreement between Syn2b and dnadiff on this subset (majority frame):")
    ok = out.dropna(subset=["dd_major", "s2b_major"])
    if len(ok) > 2:
        print(f"  n={len(ok)}  r={np.corrcoef(ok.dd_major, ok.s2b_major)[0,1]:.4f}  "
              f"bias={float((ok.s2b_major-ok.dd_major).mean()):+.4f}  "
              f"SD={float((ok.s2b_major-ok.dd_major).std()):.4f}")


if __name__ == "__main__":
    main()
