#!/usr/bin/env python3
"""Re-analyze cagPAI-disease associations with crude OR + CMH side by side.

Outputs a table with, for each 2x2 contrast:
  - crude OR and 95% CI (Woolf)
  - chi-square p (no continuity correction) and p with continuity correction
  - CMH statistic, p-value, and Mantel-Haenszel OR
  - Breslow-Day test of homogeneity of ORs across FastBAPS strata
  - per-stratum cell counts and informative (discordant) counts
"""
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chi2 as chi2_dist
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STATES = ROOT / "results" / "cagpai_states_extended_filtered.tsv"
META = Path("/Volumes/MoneyCat/Data/song_2026_hpylori/metadata.csv")
OUT = ROOT / "results" / "cagpai_association_crude_and_stratified.tsv"


def breslow_day(tables_3d):
    """Return Breslow-Day chi2 and p for homogeneity of OR across strata.

    tables_3d: numpy array of shape (K, 2, 2).
    Implementation follows Breslow & Day (1980), Iterative weighted sum.
    """
    K = tables_3d.shape[0]
    if K < 2:
        return np.nan, np.nan

    A = tables_3d[:, 0, 0].astype(float)
    B = tables_3d[:, 0, 1].astype(float)
    C = tables_3d[:, 1, 0].astype(float)
    D = tables_3d[:, 1, 1].astype(float)
    N = A + B + C + D

    # MH common OR
    R = (A * D / N).sum()
    S = (B * C / N).sum()
    if S == 0:
        return np.nan, np.nan
    psi = R / S

    # Solve each stratum's expected a_i under common OR psi
    # a_i is root of quadratic: psi * a_i * (n_i - m_i - n2_i + a_i)
    #                         - (m_i - a_i) * (n2_i - a_i) = 0
    # where m_i = A+B, n2_i = A+C, n_i = N
    m = A + B
    n2 = A + C

    def solve_a(m_i, n2_i, n_i, psi):
        # Coefficients of quadratic in a:
        # a^2*(psi-1) + a*(psi*(n-n2-m) + n2 + m) - m*n2 = 0
        qa = psi - 1.0
        qb = psi * (n_i - n2_i - m_i) + n2_i + m_i
        qc = -m_i * n2_i
        if qa == 0:
            return -qc / qb
        disc = qb * qb - 4 * qa * qc
        if disc < 0:
            return m_i * n2_i / n_i  # fallback to expected
        roots = np.roots([qa, qb, qc])
        # choose root in [max(0, m_i+n2_i-n_i), min(m_i, n2_i)]
        lo = max(0.0, m_i + n2_i - n_i)
        hi = min(m_i, n2_i)
        valid = [r.real for r in roots if np.isreal(r) and lo <= r.real <= hi]
        if valid:
            return valid[0]
        return m_i * n2_i / n_i

    a_hat = np.array([solve_a(m[i], n2[i], N[i], psi) for i in range(K)])
    b_hat = m - a_hat
    c_hat = n2 - a_hat
    d_hat = N - m - n2 + a_hat

    # Variance of a_i under H0 (common OR)
    # 1 / Var(a_i) = 1/a + 1/b + 1/c + 1/d
    var_a = 1.0 / (1.0 / a_hat + 1.0 / b_hat + 1.0 / c_hat + 1.0 / d_hat)
    # Breslow-Day statistic
    bd_stat = ((A - a_hat) ** 2 / var_a).sum()
    bd_p = 1 - chi2_dist.cdf(bd_stat, K - 1)
    return float(bd_stat), float(bd_p)


def cmh_and_crude(table_dict, label):
    """table_dict: {stratum: [[a,b],[c,d]]} case+/- vs control+/-."""
    strata = list(table_dict.keys())
    tables = np.array([table_dict[s] for s in strata], dtype=float)

    # crude table = sum over strata
    crude = tables.sum(axis=0)
    a, b, c, d = crude[0, 0], crude[0, 1], crude[1, 0], crude[1, 1]
    n = crude.sum()

    crude_or = (a * d) / (b * c) if b * c > 0 else np.inf
    se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    ci_lo = np.exp(np.log(crude_or) - 1.96 * se_log_or)
    ci_hi = np.exp(np.log(crude_or) + 1.96 * se_log_or)

    chi2_val, p_crude, _, _ = chi2_contingency(crude, correction=False)
    chi2_cc, p_crude_cc, _, _ = chi2_contingency(crude, correction=True)

    # CMH manual formula
    A = tables[:, 0, 0]
    B = tables[:, 0, 1]
    C = tables[:, 1, 0]
    D = tables[:, 1, 1]
    N = A + B + C + D
    E_A = (A + B) * (A + C) / N
    V_A = (A + B) * (C + D) * (A + C) * (B + D) / (N ** 2 * (N - 1))
    cmh_stat = ((A - E_A).sum()) ** 2 / V_A.sum()
    cmh_p = 1 - chi2_dist.cdf(cmh_stat, 1)

    R = (A * D / N).sum()
    S = (B * C / N).sum()
    or_mh = R / S if S > 0 else np.inf

    # Breslow-Day
    bd_stat, bd_p = breslow_day(tables)

    # Per-stratum details: discordant counts and expected under independence
    per_stratum = []
    for i, s in enumerate(strata):
        ai, bi, ci, di = int(A[i]), int(B[i]), int(C[i]), int(D[i])
        ni = ai + bi + ci + di
        ei = (ai + bi) * (ai + ci) / ni
        # informative = discordant pairs = b + c
        informative = bi + ci
        per_stratum.append(
            f"{s}:a{ai}:b{bi}:c{ci}:d{di}:exp_a{ei:.1f}:inf{informative}"
        )

    return {
        "label": label,
        "n": int(n),
        "crude_a": int(a),
        "crude_b": int(b),
        "crude_c": int(c),
        "crude_d": int(d),
        "crude_OR": crude_or,
        "crude_OR_lo": ci_lo,
        "crude_OR_hi": ci_hi,
        "crude_chi2": float(chi2_val),
        "crude_p": float(p_crude),
        "crude_p_cc": float(p_crude_cc),
        "CMH_stat": float(cmh_stat),
        "CMH_p": float(cmh_p),
        "OR_MH": or_mh,
        "BreslowDay_stat": bd_stat,
        "BreslowDay_p": bd_p,
        "strata_detail": "; ".join(per_stratum),
    }


def build_table_dict(df, case_col, outcome_col, stratum_col="fastbaps"):
    """Build {stratum: [[a,b],[c,d]]} from a DataFrame.

    a = case & outcome
    b = case & not outcome
    c = control & outcome
    d = control & not outcome
    """
    table_dict = {}
    for st in sorted(df[stratum_col].unique()):
        s = df[df[stratum_col] == st]
        a = ((s[case_col] == True) & (s[outcome_col] == True)).sum()
        b = ((s[case_col] == True) & (s[outcome_col] == False)).sum()
        c = ((s[case_col] == False) & (s[outcome_col] == True)).sum()
        d = ((s[case_col] == False) & (s[outcome_col] == False)).sum()
        if a + b > 0 and c + d > 0:
            table_dict[st] = [[int(a), int(b)], [int(c), int(d)]]
    return table_dict


def main():
    states = pd.read_csv(STATES, sep="\t")
    meta = pd.read_csv(META)
    df = states.merge(meta, left_on="genome", right_on="assembly", how="left")
    df["group"] = df["group"].astype(str).str.strip()

    # presence: complete (collinear+rearranged) vs empty/partial
    df["present"] = df["status_extended"].isin(["complete_collinear", "complete_rearranged"])

    results = []

    # 1. presence GC vs NAG
    sub = df[df["group"].isin(["GC", "NAG"])].copy()
    sub["case"] = sub["group"] == "GC"
    td = build_table_dict(sub, "case", "present")
    results.append(cmh_and_crude(td, "presence GC vs NAG"))

    # 2. rearrangement GC vs NAG (among present)
    sub = df[df["group"].isin(["GC", "NAG"]) & df["present"]].copy()
    sub["case"] = sub["group"] == "GC"
    sub["rearranged"] = sub["status_extended"] == "complete_rearranged"
    td = build_table_dict(sub, "case", "rearranged")
    results.append(cmh_and_crude(td, "rearrangement GC vs NAG"))

    # 3. presence advanced (GC/IM) vs early (AG/NAG)
    sub = df[df["group"].isin(["GC", "IM", "AG", "NAG"])].copy()
    sub["case"] = sub["group"].isin(["GC", "IM"])
    td = build_table_dict(sub, "case", "present")
    results.append(cmh_and_crude(td, "presence advanced vs early"))

    # 4. rearrangement advanced vs early (among present)
    sub = df[df["group"].isin(["GC", "IM", "AG", "NAG"]) & df["present"]].copy()
    sub["case"] = sub["group"].isin(["GC", "IM"])
    sub["rearranged"] = sub["status_extended"] == "complete_rearranged"
    td = build_table_dict(sub, "case", "rearranged")
    results.append(cmh_and_crude(td, "rearrangement advanced vs early"))

    res_df = pd.DataFrame(results)
    print(res_df.to_string(index=False))
    res_df.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
