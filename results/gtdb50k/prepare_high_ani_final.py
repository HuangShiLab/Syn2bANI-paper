#!/usr/bin/env python3
"""Prepare the final high-ANI train/test split for the unified GTDB-R207 benchmark.

Inputs:
  results/gtdb50k/high_ani_truth.tsv   dnadiff 1-to-1 ANIm truth for candidates
  results/gtdb50k/high_ani_candidates.tsv  original candidate strata

Output:
  results/gtdb50k/high_ani_pairs_final.tsv
    pairid, q_acc, r_acc, band, split, anim_ani, af_tier, stratum

Rules:
  - Keep only pairs with anim_ani >= 95.
  - Assign band: 95-97 or 97-100.
  - Genome-level bidirectional train/test split (no genome appears in both).
  - Stratify by band: each genome gets a primary band by majority vote among
    its pairs, then bands are split 60/40 by genome count.
  - A pair is train only if both genomes are in train; test only if both are in
    test; mixed pairs are dropped to prevent leakage.
"""
import os
import random
import argparse
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")

RANDOM_SEED = 20240825
TRAIN_FRAC = 0.60


def assign_band(ani):
    if pd.isna(ani):
        return None
    if 95.0 <= ani < 97.0:
        return "95-97"
    if 97.0 <= ani <= 100.0:
        return "97-100"
    return None


def genome_band_vote(df):
    """Return primary band for each genome by majority vote."""
    votes = {}
    for _, r in df.iterrows():
        for acc, band in [(r["q_acc"], r["band"]), (r["r_acc"], r["band"])]:
            votes.setdefault(acc, []).append(band)
    primary = {}
    for acc, bands in votes.items():
        primary[acc] = pd.Series(bands).mode()[0]
    return primary


def split_genomes_stratified(df, train_frac, seed):
    """Assign each genome to train/test stratified by primary band."""
    primary = genome_band_vote(df)
    genomes = pd.DataFrame({"acc": list(primary.keys()),
                            "primary_band": [primary[a] for a in primary]})
    rng = random.Random(seed)
    assignment = {}
    for band, sub in genomes.groupby("primary_band"):
        accs = sub["acc"].tolist()
        rng.shuffle(accs)
        n_train = max(1, int(round(len(accs) * train_frac)))
        for a in accs[:n_train]:
            assignment[a] = "train"
        for a in accs[n_train:]:
            assignment[a] = "test"
    return assignment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default=os.path.join(RES, "high_ani_truth.tsv"))
    parser.add_argument("--candidates", default=os.path.join(RES, "high_ani_candidates.tsv"))
    parser.add_argument("--out", default=os.path.join(RES, "high_ani_pairs_final.tsv"))
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--train-frac", type=float, default=TRAIN_FRAC)
    args = parser.parse_args()

    truth = pd.read_csv(args.truth, sep="\t")
    # parse pairid from truth into q_acc/r_acc if not present
    if "q_acc" not in truth.columns:
        parts = truth["pairid"].str.split("__", n=1, expand=True)
        truth["q_acc"] = parts[0]
        truth["r_acc"] = parts[1]

    truth["band"] = truth["anim_ani"].apply(assign_band)
    df = truth[truth["band"].notna()].copy()
    n_before = len(df)
    df = df.drop_duplicates(subset=["pairid"], keep="first")
    if len(df) < n_before:
        print(f"dropped {n_before - len(df):,} duplicate pairid rows")

    # attach original stratum label when available
    if os.path.exists(args.candidates):
        cand = pd.read_csv(args.candidates, sep="\t")
        cand["pairid"] = cand["q_acc"] + "__" + cand["r_acc"]
        df = df.merge(cand[["pairid", "stratum"]], on="pairid", how="left")
    else:
        df["stratum"] = np.nan
        print(f"warning: candidates file not found: {args.candidates}")

    # genome-level split
    assignment = split_genomes_stratified(df, args.train_frac, args.seed)
    df["split"] = df.apply(lambda r: (
        "train" if assignment.get(r["q_acc"]) == "train" and assignment.get(r["r_acc"]) == "train"
        else "test" if assignment.get(r["q_acc"]) == "test" and assignment.get(r["r_acc"]) == "test"
        else "mixed"
    ), axis=1)

    kept = df[df["split"] != "mixed"].copy()
    dropped = df[df["split"] == "mixed"]

    # deterministic ordering
    kept = kept.sort_values(["band", "pairid"]).reset_index(drop=True)

    cols = ["pairid", "q_acc", "r_acc", "band", "split", "anim_ani",
            "af_tier", "stratum"]
    kept = kept[cols]
    kept.to_csv(args.out, sep="\t", index=False, float_format="%.4f")

    print(f"input truth rows: {len(truth):,}")
    print(f"ANI>=95 pairs: {len(df):,}")
    print(f"mixed (dropped): {len(dropped):,}")
    print(f"final kept: {len(kept):,}")
    print("\nsplit by band:")
    print(kept.groupby(["band", "split"]).size().unstack(fill_value=0))
    print(f"\nunique genomes: {len(set(kept['q_acc']) | set(kept['r_acc'])):,}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
