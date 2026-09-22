#!/usr/bin/env python3
"""Test whether island-altering cagPAI locality tracks gastric-cancer stage.

Among complete-marker genomes, locality has four levels:
island_internal, island_boundary, island_spanned and complete_collinear.
The prespecified binary contrast defines an island-altering call as
island_internal or island_boundary. The two disease contrasts are the same
ones used for cagPAI presence/rearrangement: GC versus NAG and advanced
(GC/IM) versus early (AG/NAG). Each binary test is reported crudely and
after Cochran-Mantel-Haenszel stratification by FastBAPS lineage.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
LOCALITY = RESULTS / "cagpai_states_locality.tsv"
MERGED = RESULTS / "metadata_assoc" / "merged_genome_table.tsv"
OUT = HERE.parent / "results" / "cagpai_locality_disease_association.tsv"


def or_ci_2x2(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    est = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return est, np.exp(np.log(est) - 1.96 * se), np.exp(np.log(est) + 1.96 * se)


def breslow_day_p(tables):
    """Breslow-Day homogeneity test using the Mantel-Haenszel pooled OR."""
    tables = [tuple(map(int, np.asarray(t).ravel())) for t in tables]
    num = sum(a * d / (a + b + c + d) for a, b, c, d in tables)
    den = sum(b * c / (a + b + c + d) for a, b, c, d in tables)
    if den <= 0:
        return np.nan
    psi = num / den
    stat = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        n1, n2 = a + b, c + d
        m1, m2 = a + c, b + d
        aa, bb, cc = psi - 1, -(psi * (m1 + n1) + (n2 - m1)), psi * m1 * n1
        disc = bb * bb - 4 * aa * cc
        if disc < 0 or aa == 0:
            return np.nan
        expected = (-bb - np.sqrt(disc)) / (2 * aa)
        if not (max(0, m1 - n2) - 1e-9 <= expected <= min(n1, m1) + 1e-9):
            expected = (-bb + np.sqrt(disc)) / (2 * aa)
        var = 1 / (1 / expected + 1 / (n1 - expected) +
                   1 / (m1 - expected) + 1 / (n2 - m1 + expected))
        stat += (a - expected) ** 2 / var
    from scipy.stats import chi2
    return float(chi2.sf(stat, len(tables) - 1))


def binary_tests(sub, label):
    tab = pd.crosstab(sub["altered"], sub["outcome"])
    a = int(tab.loc[1, 1])
    b = int(tab.loc[1, 0])
    c = int(tab.loc[0, 1])
    d = int(tab.loc[0, 0])
    est, lo, hi = or_ci_2x2(a, b, c, d)
    p_fisher = fisher_exact([[a, b], [c, d]])[1]

    strata = []
    for lineage, g in sub.groupby("fastbaps"):
        ct = pd.crosstab(g["altered"], g["outcome"])
        if ct.shape == (2, 2):
            strata.append((lineage, ct.values))
    st = StratifiedTable([t for _, t in strata])
    cmh_or, cmh_ci = st.oddsratio_pooled, st.oddsratio_pooled_confint()
    cmh_p = st.test_null_odds().pvalue
    bd_p = breslow_day_p([t for _, t in strata])
    return {
        "test": label,
        "n": a + b + c + d,
        "altered_cases": a,
        "nonaltered_cases": b,
        "altered_controls": c,
        "nonaltered_controls": d,
        "crude_OR": est,
        "crude_OR_lo": lo,
        "crude_OR_hi": hi,
        "crude_p_fisher": p_fisher,
        "CMH_OR": cmh_or,
        "CMH_OR_lo": cmh_ci[0],
        "CMH_OR_hi": cmh_ci[1],
        "CMH_p": cmh_p,
        "BreslowDay_p": bd_p,
        "n_strata": len(strata),
        "strata_detail": "; ".join(
            f"{lin}:a{t[0, 0]}:b{t[0, 1]}:c{t[1, 0]}:d{t[1, 1]}"
            for lin, t in strata
        ),
    }


def main():
    loc = pd.read_csv(LOCALITY, sep="\t")
    merged = pd.read_csv(MERGED, sep="\t")
    df = loc.merge(merged, left_on="genome", right_on="GCA_accession",
                   suffixes=("_locality", ""))
    df = df[df["status_extended_locality"].str.startswith("complete")].copy()
    df["locality_level"] = df["locality"].fillna("complete_collinear")
    df["altered"] = df["locality_level"].isin(
        ["island_internal", "island_boundary"]).astype(int)

    full = pd.crosstab(df["group"], df["locality_level"])
    chi2, p_full, dof, expected = chi2_contingency(full)

    rows = [{
        "test": "four-stage x four-level locality (global chi-square)",
        "n": int(full.values.sum()),
        "chi2": float(chi2),
        "df": int(dof),
        "p": float(p_full),
        "min_expected": float(expected.min()),
        "table": full.to_csv(sep=":", lineterminator="|"),
    }]

    gc = df[df["group"].isin(["GC", "NAG"])].copy()
    gc["outcome"] = (gc["group"] == "GC").astype(int)
    rows.append(binary_tests(gc, "island-altering locality: GC vs NAG"))

    adv = df[df["group"].isin(["GC", "IM", "AG", "NAG"])].copy()
    adv["outcome"] = adv["group"].isin(["GC", "IM"]).astype(int)
    rows.append(binary_tests(adv, "island-altering locality: advanced vs early"))

    locality_ps = [
        rows[0]["p"],
        rows[1]["crude_p_fisher"],
        rows[2]["crude_p_fisher"],
    ]
    from statsmodels.stats.multitest import multipletests
    locality_q = multipletests(locality_ps, method="fdr_bh")[1]
    rows[0]["q_BH"] = locality_q[0]
    rows[1]["q_BH"] = locality_q[1]
    rows[2]["q_BH"] = locality_q[2]

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"complete-marker isolates: {len(df)}")
    print(full.to_string())
    print(f"global chi-square: chi2={chi2:.3f}, df={dof}, p={p_full:.6g}, "
          f"min expected={expected.min():.3g}")
    print(out.drop(columns=["table", "strata_detail"]).to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
