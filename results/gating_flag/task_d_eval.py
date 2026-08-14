#!/usr/bin/env python3
"""Task D — end-to-end before/after evaluation of the gated estimator and the
recalibrated flag on the regenerated GTDB-ANIm feature matrix.

Input: ../anim_truth_2074_gated.tsv (new binary, has ani_gated/gate columns),
joined to ANIm truth exactly as in prep_data.py. Writes before_after.tsv.
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
BANDS = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]


def main():
    new = pd.read_csv(HERE.parent / "anim_truth_2074_gated.tsv", sep="\t")
    m = pd.read_csv(HERE.parent / "anim_2074_acc2seqid.tsv", sep="\t", header=None,
                    names=["accession", "seqid"])
    ev = pd.read_csv(HERE.parent / "panel_by_band/eval_pairs.tsv", sep="\t")
    f = new.merge(m, left_on="query", right_on="seqid").drop(columns=["seqid"]) \
           .rename(columns={"accession": "query_asm"})
    f = f.merge(m, left_on="reference", right_on="seqid").drop(columns=["seqid"]) \
         .rename(columns={"accession": "ref_asm"})
    df = f.merge(ev[["query_asm", "ref_asm", "band", "anim_ani"]],
                 on=["query_asm", "ref_asm"], how="inner")
    df = df[np.isfinite(df["ani"])].copy()
    truth = df["anim_ani"]
    print(f"pairs: {len(df)}; gate fallback on {(df['gate'] == 'uniform_fallback').sum()}"
          f" ({(df['gate'] == 'uniform_fallback').mean():.1%})")

    rows = {}
    for col, label in [("ani", "gamma (before)"), ("ani_uniform", "uniform"),
                       ("ani_gated", "gated (after)")]:
        err = (df[col] - truth).abs()
        rows[label] = {"all": err.mean(),
                       **{b: err[df["band"] == b].mean() for b in BANDS},
                       "bias": (df[col] - truth).mean()}
    oracle = np.minimum((df["ani"] - truth).abs(), (df["ani_uniform"] - truth).abs())
    rows["oracle"] = {"all": oracle.mean(),
                      **{b: oracle[df["band"] == b].mean() for b in BANDS},
                      "bias": np.nan}
    tab = pd.DataFrame(rows).T[["all", *BANDS, "bias"]]
    print("\nMAE vs ANIm (ANI points):")
    print(tab.round(3).to_string())
    tab.to_csv(HERE / "before_after.tsv", sep="\t")

    # new flag behaviour (gated error)
    df["err"] = (df["ani_gated"] - truth).abs()
    print("\nnew flag vs gated error:")
    for fl in ["ok", "INCONSISTENT", "BELOW_DETECTION"]:
        mm = df["flag"] == fl
        print(f"  {fl:16s} n={mm.sum():5d}  MAE {df.loc[mm, 'err'].mean():.3f}")


if __name__ == "__main__":
    main()
