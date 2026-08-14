#!/usr/bin/env python3
"""Task A (likelihood-exact gates) — evaluate LRT/BIC/boundary-mixture and
discrepancy gates against ANIm truth using the exact refits in
refit_cache.tsv (run task_a_refit.py first; needs strata_2074/).

Conclusion (see RULES.md): no likelihood-significance gate beats the simple
effect-size rule |ani_from_loss - ani_from_hist| > 5 points.
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
BANDS = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]


def main():
    rc = pd.read_csv(HERE / "refit_cache.tsv", sep="\t")
    df = pd.read_csv(HERE / "gtdb_anim_joined.tsv", sep="\t")
    df = df.merge(rc, left_on=["query_asm", "ref_asm"], right_on=["q_acc", "r_acc"],
                  how="left", suffixes=("", "_rc"))
    df = df[np.isfinite(df["ani"])].copy()
    truth = df["anim_ani"]
    df["gap_lh"] = (df["ani_from_loss"] - df["ani_from_hist"]).abs()

    def ev(use_gamma):
        est = df["ani"].where(use_gamma, df["ani_uniform"])
        err = (est - truth).abs()
        return {"all": err.mean(),
                **{b: err[df["band"] == b].mean() for b in BANDS},
                "frac_gamma": float(use_gamma.fillna(False).mean())}

    gates = {
        "always_gamma (current, LRT>3.841)": pd.Series(True, index=df.index),
        "always_uniform": pd.Series(False, index=df.index),
        "BIC: lrt > ln(n_tags)": df["lrt"] > np.log(df["n_tags_rc"]),
        "chi-bar-square boundary: lrt > 2.706": df["lrt"] > 2.706,
        "LRT > 6": df["lrt"] > 6,
        "LRT > 10": df["lrt"] > 10,
        "LRT > 20": df["lrt"] > 20,
        "|gamma-uniform| <= 2*SE else uniform":
            (df["ani"] - df["ani_uniform"]).abs() <= 2 * df["se_u"].clip(lower=0.05),
        "gap_lh <= 4 else uniform": df["gap_lh"] <= 4,
        "gap_lh <= 5 else uniform (CHOSEN)": df["gap_lh"] <= 5,
        "gap_lh <= 6 else uniform": df["gap_lh"] <= 6,
        "LRT>3.841 & gap_lh<=5": (df["lrt"] > 3.841) & (df["gap_lh"] <= 5),
        "BIC & gap_lh<=5": (df["lrt"] > np.log(df["n_tags_rc"])) & (df["gap_lh"] <= 5),
    }
    rows = {name: ev(g) for name, g in gates.items()}
    tab = pd.DataFrame(rows).T[["all", *BANDS, "frac_gamma"]].sort_values("all")
    print(tab.round(3).to_string())
    tab.to_csv(HERE / "task_a_lrt_gates.tsv", sep="\t")


if __name__ == "__main__":
    main()
