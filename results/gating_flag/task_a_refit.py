#!/usr/bin/env python3
"""Refit all 2,074 GTDB pairs from their strata and cache the exact
likelihood quantities (uniform fit, gamma fit, LRT, raw shape, NLLs)
to refit_cache.tsv. Also refits the 15 mid-ANI pairs and oral/gut if
strata are available for them.
"""
import math
import pathlib

import numpy as np
import pandas as pd

import het_fit

HERE = pathlib.Path(__file__).resolve().parent


def main():
    strata = het_fit.load_strata_dir(HERE / "strata_2074")
    rows = []
    for (q, r), st in strata.items():
        h = het_fit.estimate_het(st)
        ret = het_fit.expected_retention(h["ani_u"], st)
        rows.append(dict(
            q_acc=q, r_acc=r,
            ani_u=h["ani_u"] * 100, se_u=h["se_u"] * 100,
            ani_g_raw=het_fit.het_ani(h["d"], h["alpha_raw"]) * 100,
            alpha_raw=h["alpha_raw"], lrt=h["lrt"],
            nll_het=h["nll_het"], nll_null=h["nll_null"],
            n_tags=h["n_tags"], retention=ret,
        ))
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "refit_cache.tsv", sep="\t", index=False)
    print(f"refit {len(df)} pairs -> refit_cache.tsv")


if __name__ == "__main__":
    main()
