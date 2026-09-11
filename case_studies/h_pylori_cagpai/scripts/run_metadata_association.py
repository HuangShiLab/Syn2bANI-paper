#!/usr/bin/env python3
"""Association analysis: H. pylori structural-variation features vs new isolate
metadata (cagA, vacA, CARD antibiotic-resistance genotype).

Merges per-genome cagPAI status, SV burden (vs 26695), cagA/vacA typing and
CARD resistance calls with isolate metadata, then runs exploratory association
tests with lineage (fastbaps) stratification and BH-FDR correction across the
full test family.

Reads from /Volumes/MoneyCat/Data/song_2026_hpylori (read-only).
Writes to case_studies/h_pylori_cagpai/results/metadata_assoc/.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------- paths
DATA = Path("/Volumes/MoneyCat/Data/song_2026_hpylori")
SCRIPT_DIR = Path(__file__).resolve().parent
OUT = SCRIPT_DIR.parent / "results" / "metadata_assoc"
OUT.mkdir(parents=True, exist_ok=True)

LARGE_SV_THRESH = 50_000  # bp

# drug-class keywords in CARD Drug_class -> (column name, representative drug)
DRUG_CLASSES = {
    "macrolide antibiotic": ("res_macrolide", "clarithromycin"),
    "fluoroquinolone antibiotic": ("res_fluoroquinolone", "levofloxacin"),
    "tetracycline antibiotic": ("res_tetracycline", "tetracycline"),
    "nitroimidazole antibiotic": ("res_nitroimidazole", "metronidazole"),
    "penam": ("res_betalactam", "amoxicillin"),  # covers penam/cephalosporin/monobactam/penem
}

rng = np.random.default_rng(42)

# ================================================================= load


def load_inputs():
    caga = pd.read_excel(DATA / "cagA_final_GCA528.xlsx", sheet_name="Sheet1")
    vaca = pd.read_excel(DATA / "vacA_typing_GCA528.xlsx")
    card = pd.read_excel(DATA / "CARD_detailed_results.xlsx", sheet_name="Sheet2")
    meta = pd.read_csv(DATA / "metadata.csv")
    states = pd.read_csv(DATA / "cagpai_status" / "cagpai_states_extended.tsv",
                         sep="\t")
    return caga, vaca, card, meta, states


def sv_burden_table():
    """Per-genome SV burden from circular-origin-filtered BED files."""
    rows = []
    beddir = DATA / "struct_vs_26695_filtered"
    for bed in sorted(beddir.glob("GCA_*.vs_hp26695.bed")):
        gca = bed.name.split(".vs_")[0]
        df = pd.read_csv(bed, sep="\t", header=None, encoding="latin-1",
                         names=["chrom", "start", "end", "name", "score",
                                "strand", "ts", "te", "color"])
        df["svtype"] = df["name"].str.split("_").str[0]
        df["span"] = df["end"] - df["start"]
        rec = {
            "GCA_accession": gca,
            "sv_total_n": len(df),
            "sv_total_span_bp": int(df["span"].sum()),
            "sv_large_n": int((df["span"] > LARGE_SV_THRESH).sum()),
            "sv_large_span_bp": int(df.loc[df["span"] > LARGE_SV_THRESH, "span"].sum()),
        }
        for t in ["DEL", "INS", "INV", "TRA"]:
            rec[f"sv_n_{t}"] = int((df["svtype"] == t).sum())
            rec[f"sv_span_{t}_bp"] = int(df.loc[df["svtype"] == t, "span"].sum())
        rows.append(rec)
    return pd.DataFrame(rows)


def card_indicators(card):
    """Per-genome binary resistance indicator per drug class.

    Only Perfect/Strict hits count as detected. hp1181 is an intrinsic MFS
    efflux pump native to H. pylori (near-100% prevalence); it still maps to
    fluoroquinolone/tetracycline/nitroimidazole Drug_class strings in CARD,
    and we keep it but flag prevalence in the report.
    """
    det = card[card["Cut_Off"].isin(["Perfect", "Strict"])].copy()
    recs = {}
    for _, r in det.iterrows():
        g = r["GCA_accession"]
        rec = recs.setdefault(g, {"any_mutation_res": 0})
        if pd.notna(r["Drug_class"]):
            for cls, (col, _) in DRUG_CLASSES.items():
                if cls in r["Drug_class"]:
                    rec[col] = 1
        if pd.notna(r.get("Resistance_mutation")):
            rec["any_mutation_res"] = 1
            if "macrolide" in str(r["Drug_class"]):
                rec["res_macrolide_mutation"] = 1
    rows = []
    for g, rec in recs.items():
        rec["GCA_accession"] = g
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("GCA_accession")
    for col in [v[0] for v in DRUG_CLASSES.values()] + ["any_mutation_res",
                                                        "res_macrolide_mutation"]:
        if col not in out:
            out[col] = 0
        out[col] = out[col].fillna(0).astype(int)
    return out.reset_index()


def merge_all(caga, vaca, card, meta, states, sv):
    df = meta.rename(columns={"assembly": "GCA_accession"}).copy()
    df = df.merge(states.rename(columns={"genome": "GCA_accession"}),
                  on="GCA_accession", how="left")
    df = df.merge(caga, on="GCA_accession", how="left")
    df = df.merge(vaca, on="GCA_accession", how="left")
    df = df.merge(card, on="GCA_accession", how="left")
    df = df.merge(sv, on="GCA_accession", how="left")

    # derived cagPAI variables
    df["cagpai_complete"] = (df["status"] == "complete").astype(float)
    df.loc[df["status"].isna(), "cagpai_complete"] = np.nan
    df["cagpai_rearranged"] = np.where(
        df["status_extended"].isin(["complete_rearranged", "complete_collinear"]),
        (df["status_extended"] == "complete_rearranged").astype(float), np.nan)

    # derived cagA variables
    df["cagA_detected"] = np.where(
        df["cagA_result"].isin(["Detected", "Not_detected"]),
        (df["cagA_result"] == "Detected").astype(float), np.nan)

    def _gtype(g):
        if pd.isna(g):
            return np.nan
        return g if g in ("ABC", "ABD", "ABCC") else ("multi" if len(g) > 4 else "other")

    df["cagA_genotype_grp"] = df["cagA_genotype"].map(_gtype)
    df["east_west"] = df["east_west_type"].where(
        df["east_west_type"].isin(["Western_like", "EastAsian_like"]))
    df["east_asian"] = np.where(df["east_west"].notna(),
                                (df["east_west"] == "EastAsian_like").astype(float),
                                np.nan)

    # drug resistance indicators -> fill absent CARD hit = 0
    for col in [v[0] for v in DRUG_CLASSES.values()] + ["any_mutation_res",
                                                        "res_macrolide_mutation"]:
        df[col] = df[col].fillna(0).astype(int)

    # SV burden missing -> 0 SVs is a real value only if BED exists; all 528
    # have BEDs, so no fill needed, but guard anyway:
    for col in ["sv_total_n", "sv_large_n"]:
        df[col] = df[col].fillna(0).astype(int)

    med = df["sv_total_n"].median()
    df["sv_burden_high"] = (df["sv_total_n"] > med).astype(float)
    df["sv_large_high"] = (df["sv_large_n"] >
                           df["sv_large_n"].median()).astype(float)
    df["sv_inv_high"] = (df["sv_n_INV"] > df["sv_n_INV"].median()).astype(float)
    df["sv_span_high"] = (df["sv_total_span_bp"] >
                          df["sv_total_span_bp"].median()).astype(float)
    return df

# ======================================================= stats helpers


def or_ci_2x2(a, b, c, d):
    """OR = (a/b)/(c/d) with Woolf 95% CI; Haldane-Anscombe if zero cell."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se)


def breslow_day_p(tables):
    """Manual Breslow-Day homogeneity test (df = k-1) from 2x2 tables given
    as (a, b, c, d) tuples. Uses the Mantel-Haenszel pooled OR; returns
    np.nan if the pooled OR is not estimable."""
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
        A, B, C = psi - 1, -(psi * (m1 + n1) + (n2 - m1)), psi * m1 * n1
        disc = B * B - 4 * A * C
        if disc < 0 or A == 0:
            return np.nan
        ea = (-B - np.sqrt(disc)) / (2 * A)  # root in admissible range
        if not (max(0, m1 - n2) - 1e-9 <= ea <= min(n1, m1) + 1e-9):
            ea = (-B + np.sqrt(disc)) / (2 * A)
        var = 1 / (1 / ea + 1 / (n1 - ea) + 1 / (m1 - ea) + 1 / (n2 - m1 + ea))
        stat += (a - ea) ** 2 / var
    return float(stats.chi2.sf(stat, len(tables) - 1))


def stratified_2x2(sub, exp_col, out_col):
    """CMH pooled OR + p, Breslow-Day homogeneity, per statsmodels/manual.

    sub: complete-case df with two binary columns. Returns dict or None if
    any stratum is degenerate for CMH.
    """
    tables = []
    for _, g in sub.groupby("fastbaps"):
        ct = pd.crosstab(g[exp_col], g[out_col])
        if ct.shape != (2, 2):
            continue
        tables.append(ct.values)
    if len(tables) < 2:
        return None
    try:
        st_ = StratifiedTable(tables)
        orr, ci = st_.oddsratio_pooled, st_.oddsratio_pooled_confint()
        p_cmh = st_.test_null_odds().pvalue  # Mantel-Haenszel conditional indep.
        bd_p = breslow_day_p([(t[0, 0], t[0, 1], t[1, 0], t[1, 1])
                              for t in tables])
        return {"or": float(orr), "ci_lo": float(ci[0]), "ci_hi": float(ci[1]),
                "p": float(p_cmh), "bd_p": float(bd_p), "n_strata": len(tables)}
    except Exception:
        return None


def logit_strat(sub, out_col, exp_col):
    """Logistic regression of binary outcome on exposure + fastbaps dummies."""
    d = sub.dropna(subset=[out_col, exp_col, "fastbaps"]).copy()
    X = pd.get_dummies(d[["exp", "fastbaps"]], drop_first=True, dtype=float)
    X = sm.add_constant(X)
    try:
        m = sm.Logit(d["out"].astype(float), X).fit(disp=0)
        i = list(X.columns).index("exp")
        b = m.params.iloc[i]
        ci = m.conf_int().iloc[i]
        return {"or": float(np.exp(b)), "ci_lo": float(np.exp(ci[0])),
                "ci_hi": float(np.exp(ci[1])), "p": float(m.pvalues.iloc[i])}
    except Exception:
        return None


def run_test(rows, category, exposure, outcome, contrast, sub, exp_col,
             out_col, binary=True, effect_name="OR"):
    """2x2 (or binary-vs-multilevel: chi2 + collapsed 2x2) association."""
    d = sub.dropna(subset=[exp_col, out_col, "fastbaps"]).copy()
    n = len(d)
    rec = {"test_id": len(rows) + 1, "category": category, "exposure": exposure,
           "outcome": outcome, "contrast": contrast, "n": n,
           "effect_measure": effect_name, "est": np.nan, "ci_lo": np.nan,
           "ci_hi": np.nan, "p_value": np.nan, "method": "", "strat_or": np.nan,
           "strat_ci_lo": np.nan, "strat_ci_hi": np.nan, "strat_p": np.nan,
           "bd_p": np.nan, "min_cell": np.nan, "low_count": False, "notes": ""}

    if binary:
        d["exp"] = d[exp_col].astype(float)
        d["out"] = d[out_col].astype(float)
        ct = pd.crosstab(d["exp"], d["out"])
        if ct.shape != (2, 2):
            rec["notes"] = "degenerate table"
            rows.append(rec)
            return
        a = ct.loc[1.0, 1.0]; b = ct.loc[1.0, 0.0]
        c = ct.loc[0.0, 1.0]; dd = ct.loc[0.0, 0.0]
        orr, lo, hi = or_ci_2x2(a, b, c, dd)
        p = stats.fisher_exact([[a, b], [c, dd]])[1]
        rec.update(est=orr, ci_lo=lo, ci_hi=hi, p_value=p,
                   method="Fisher exact", min_cell=int(min(a, b, c, dd)),
                   low_count=min(a, b, c, dd) < 5)
        strat = stratified_2x2(d, "exp", "out")
        if strat:
            rec.update(strat_or=strat["or"], strat_ci_lo=strat["ci_lo"],
                       strat_ci_hi=strat["ci_hi"], strat_p=strat["p"],
                       bd_p=strat["bd_p"])
        else:
            lg = logit_strat(d, "out", "exp")
            if lg:
                rec.update(strat_or=lg["or"], strat_ci_lo=lg["ci_lo"],
                           strat_ci_hi=lg["ci_hi"], strat_p=lg["p"],
                           notes="stratified by fastbaps (logistic reg; "
                                 "CMH degenerate)")
    else:
        ct = pd.crosstab(d[exp_col], d[out_col])
        rec["min_cell"] = int(ct.values.min())
        rec["low_count"] = bool((ct.values < 5).any())
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            rec["notes"] = "degenerate table"
            rows.append(rec)
            return
        chi2, p, _, exp_ct = stats.chi2_contingency(ct)
        rec.update(p_value=p, method=f"chi2 ({ct.shape[0]-1}x{ct.shape[1]-1})")
        if (exp_ct < 5).any():
            rec["notes"] = "expected<5 in some cells"
    rows.append(rec)


def run_continuous(rows, category, exposure, outcome, contrast, sub, exp_col,
                   val_col):
    """Continuous SV-burden outcome vs binary exposure: MWU crude + logistic
    regression (high/low by median) with fastbaps covariate."""
    d = sub.dropna(subset=[exp_col, val_col, "fastbaps"]).copy()
    g1 = d.loc[d[exp_col] == 1, val_col]
    g0 = d.loc[d[exp_col] == 0, val_col]
    p = stats.mannwhitneyu(g1, g0, alternative="two-sided").pvalue
    d["out"] = (d[val_col] > d[val_col].median()).astype(float)
    d["exp"] = d[exp_col].astype(float)
    rec = {"test_id": len(rows) + 1, "category": category, "exposure": exposure,
           "outcome": outcome, "contrast": contrast, "n": len(d),
           "effect_measure": "median_diff", "est": float(g1.median() - g0.median()),
           "ci_lo": np.nan, "ci_hi": np.nan, "p_value": float(p),
           "method": "Mann-Whitney U", "strat_or": np.nan, "strat_ci_lo": np.nan,
           "strat_ci_hi": np.nan, "strat_p": np.nan, "bd_p": np.nan,
           "min_cell": int(min(len(g1), len(g0))), "low_count": min(len(g1), len(g0)) < 5,
           "notes": f"median exp={g1.median():.1f} vs unexp={g0.median():.1f}"}
    lg = logit_strat(d, "out", "exp")
    if lg:
        rec.update(strat_or=lg["or"], strat_ci_lo=lg["ci_lo"],
                   strat_ci_hi=lg["ci_hi"], strat_p=lg["p"],
                   notes=rec["notes"] + "; stratified OR = high-burden odds "
                                       "per lineage-adjusted logistic")
    rows.append(rec)

# ================================================================== main


def main():
    caga, vaca, card, meta, states = load_inputs()
    sv = sv_burden_table()
    card_ind = card_indicators(card)
    df = merge_all(caga, vaca, card_ind, meta, states, sv)

    # guard: every genome must have a real SV burden (merge key = versioned GCA)
    assert df["sv_total_span_bp"].notna().all(), "SV merge failed: check GCA keys"

    df.to_csv(OUT / "merged_genome_table.tsv", sep="\t", index=False)
    print(f"merged table: {df.shape[0]} genomes x {df.shape[1]} columns")
    print("merge losses: missing SV burden:", df["sv_total_n"].isna().sum(),
          "| missing cagPAI status:", df["status"].isna().sum(),
          "| missing metadata:", df["fastbaps"].isna().sum())

    # ------------------------------------------------ association tests
    rows = []
    d = df  # shorthand

    # --- cagPAI x cagA (cross-validation of marker panel)
    sub = d[d["cagA_result"].isin(["Detected", "Not_detected"])]
    run_test(rows, "cagPAI x cagA", "cagPAI complete (1/0)", "cagA detected",
             "complete vs empty/partial", sub, "cagpai_complete", "cagA_detected")
    run_test(rows, "cagPAI x cagA", "cagPAI status (3-level)",
             "cagA detected", "complete/partial/empty vs detected",
             sub, "status", "cagA_detected", binary=False)
    run_test(rows, "cagPAI x cagA", "cagPAI complete (1/0)", "cagA genotype",
             "ABC/ABD/ABCC/multi vs complete",
             d.dropna(subset=["cagA_genotype_grp"]),
             "cagpai_complete", "cagA_genotype_grp", binary=False)
    run_test(rows, "cagPAI x cagA", "cagPAI complete (1/0)", "east-west type",
             "EastAsian vs Western", d, "cagpai_complete", "east_asian")
    run_test(rows, "cagPAI x cagA", "cagPAI status (3-level)", "east-west type",
             "complete/partial/empty vs East/West",
             d, "status", "east_west_type", binary=False)
    sub = d[d["status_extended"].isin(["complete_rearranged", "complete_collinear"])
            & d["cagA_result"].isin(["Detected", "Not_detected"])]
    run_test(rows, "cagPAI x cagA (secondary)",
             "cagPAI rearranged (among complete)", "cagA detected",
             "rearranged vs collinear", sub, "cagpai_rearranged", "cagA_detected")
    sub = d[d["status_extended"].isin(["complete_rearranged", "complete_collinear"])]
    run_test(rows, "cagPAI x cagA (secondary)",
             "cagPAI rearranged (among complete)", "east-west type",
             "rearranged vs collinear", sub, "cagpai_rearranged", "east_asian")

    # --- cagPAI x vacA
    for vt, label in [("s_type", "vacA s1"), ("i_type", "vacA i1"),
                      ("m_type", "vacA m1")]:
        sub = d[d[vt].isin(["s1", "s2"] if vt == "s_type" else
                           (["i1", "i2"] if vt == "i_type" else ["m1", "m2"]))]
        out = (sub[vt] == ({"s_type": "s1", "i_type": "i1", "m_type": "m1"}[vt])).astype(float)
        sub = sub.assign(**{vt + "_bin": out})
        run_test(rows, "cagPAI x vacA", "cagPAI complete (1/0)", label,
                 f"{label} vs other", sub, "cagpai_complete", vt + "_bin")

    # --- SV burden x resistance
    for cls_col, (label, drug) in {
            "res_macrolide": ("macrolide resistance (CARD)", "clarithromycin"),
            "res_fluoroquinolone": ("fluoroquinolone resistance (CARD)", "levofloxacin"),
            "res_tetracycline": ("tetracycline resistance (CARD)", "tetracycline"),
            "res_nitroimidazole": ("nitroimidazole resistance (CARD)", "metronidazole"),
            "res_betalactam": ("beta-lactam resistance (CARD)", "amoxicillin"),
            "any_mutation_res": ("any mutational resistance (CARD)", "any")}.items():
        prev = d[cls_col].mean()
        if prev > 0.95 or prev < 0.05:
            rows.append({"test_id": len(rows) + 1, "category": "SV x resistance",
                         "exposure": label, "outcome": "SV burden (median split)",
                         "contrast": "skipped: prevalence outside [5%,95%]",
                         "n": int(d[cls_col].notna().sum()),
                         "effect_measure": "OR", "est": np.nan, "ci_lo": np.nan,
                         "ci_hi": np.nan, "p_value": np.nan, "method": "skipped",
                         "strat_or": np.nan, "strat_ci_lo": np.nan,
                         "strat_ci_hi": np.nan, "strat_p": np.nan, "bd_p": np.nan,
                         "min_cell": np.nan,
                         "low_count": prev < 0.05,
                         "notes": f"prevalence={prev:.3f}; non-informative"})
            continue
        run_test(rows, "SV x resistance", label, "SV burden (median split)",
                 f"high SV burden vs {drug}-resistant", d, cls_col, "sv_burden_high")
        run_test(rows, "SV x resistance", label, "large-SV count (median split)",
                 f"high large-SV count vs {drug}-resistant", d, cls_col,
                 "sv_large_high")
        run_continuous(rows, "SV x resistance", label,
                       "SV total count (continuous)",
                       f"count vs {drug}-resistant", d, cls_col, "sv_total_n")
        run_continuous(rows, "SV x resistance", label,
                       "SV total span bp (continuous)",
                       f"span vs {drug}-resistant", d, cls_col, "sv_total_span_bp")

    # --- disease group (secondary)
    sub = d[d["group"].isin(["GC", "NAG"])].copy()
    sub["gc"] = (sub["group"] == "GC").astype(float)
    run_test(rows, "disease (secondary)", "cagPAI complete (1/0)",
             "gastric cancer (GC vs NAG)", "GC vs NAG", sub, "cagpai_complete",
             "gc")
    run_test(rows, "disease (secondary)", "macrolide resistance (CARD)",
             "gastric cancer (GC vs NAG)", "GC vs NAG", sub, "res_macrolide",
             "gc")

    res = pd.DataFrame(rows)
    # BH-FDR across the full family of executed tests (skip skipped rows)
    executed = res["p_value"].notna() | res["strat_p"].notna()
    p_all = res.loc[executed, ["p_value", "strat_p"]].min(axis=1).values
    order = np.argsort(p_all)
    q = np.empty_like(p_all)
    m = len(p_all)
    prev_q = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = m - rank
        prev_q = min(prev_q, p_all[idx] * m / i)
        q[idx] = prev_q
    res.loc[executed, "q_fdr"] = q
    res.to_csv(OUT / "association_tests.tsv", sep="\t", index=False)
    print(f"\nwrote {len(res)} tests ({int(executed.sum())} executed)")
    cols = ["category", "exposure", "outcome", "n", "est", "ci_lo", "ci_hi",
            "p_value", "strat_or", "strat_ci_lo", "strat_ci_hi", "strat_p",
            "bd_p", "q_fdr", "low_count"]
    with pd.option_context("display.width", 250):
        print(res[cols].to_string(index=False))

    # ------------------------------------------------ figures
    make_figures(df, res)
    return df, res


def make_figures(df, res):
    # (a) stacked bar: cagPAI status by east-west type and cagA genotype group
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    order_status = ["complete", "partial", "empty"]
    colors = {"complete": "#2ca02c", "partial": "#ff7f0e", "empty": "#d62728"}
    for ax, col, title in [
            (axes[0], "east_west", "cagPAI status by East-West type"),
            (axes[1], "cagA_genotype_grp", "cagPAI status by cagA genotype")]:
        sub = df.dropna(subset=[col, "status"])
        if col == "east_west":
            sub = sub[sub[col].isin(["Western_like", "EastAsian_like"])]
        ct = pd.crosstab(sub[col], sub["status"])[order_status]
        ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
        frac = ct.div(ct.sum(axis=1), axis=0)
        bottom = np.zeros(len(frac))
        for st in order_status:
            ax.bar(frac.index, frac[st], bottom=bottom, color=colors[st],
                   label=st)
            bottom += frac[st].values
        ax.set_title(title)
        ax.set_ylabel("fraction")
        ax.tick_params(axis="x", rotation=45)
    axes[0].legend(title="cagPAI status", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_cagpai_by_cagA_eastwest.png", dpi=300)
    plt.close(fig)

    # (b) SV burden by macrolide / FQ resistance
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col, title in [
            (axes[0], "res_macrolide", "SV count by macrolide resistance\n(23S rRNA mutation, CARD)"),
            (axes[1], "res_fluoroquinolone", "SV count by fluoroquinolone class\n(intrinsic hp1181 efflux, CARD)")]:
        data = [df.loc[df[col] == 0, "sv_total_n"], df.loc[df[col] == 1, "sv_total_n"]]
        ticks = [f"absent\n(n={len(data[0])})", f"present\n(n={len(data[1])})"]
        bp = ax.boxplot(data, tick_labels=ticks, showfliers=False,
                        patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#1f77b4", "#d62728"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        for i, dd in enumerate(data):
            ax.scatter(np.full(len(dd), i + 1) + rng.uniform(-0.12, 0.12, len(dd)),
                       dd, s=6, alpha=0.25, color="black")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("SV count vs 26695")
    fig.tight_layout()
    fig.savefig(OUT / "fig_sv_burden_by_resistance.png", dpi=300)
    plt.close(fig)

    # (c) forest plot of headline stratified tests
    keep = (res["strat_or"].notna() & res["category"].str.contains("cagPAI|SV")
            & ~res["exposure"].str.contains("any mutational")
            & np.isfinite(res["strat_ci_lo"]) & np.isfinite(res["strat_ci_hi"])
            & (res["strat_ci_hi"] > 0))
    forest = res[keep].copy()
    fig, ax = plt.subplots(figsize=(9, 0.5 * len(forest) + 2))
    y = np.arange(len(forest))[::-1]
    ax.errorbar(forest["strat_or"], y,
                xerr=[forest["strat_or"] - forest["strat_ci_lo"],
                      forest["strat_ci_hi"] - forest["strat_or"]],
                fmt="o", color="#1f77b4", capsize=3)
    labels = [f"{r.exposure} -> {r.outcome}" for r in forest.itertuples()]
    ax.set_yticks(y, labels, fontsize=8)
    ax.axvline(1, color="grey", ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("lineage-stratified (CMH / logistic) odds ratio")
    ax.set_title("Headline stratified associations")
    fig.tight_layout()
    fig.savefig(OUT / "fig_forest_stratified.png", dpi=300)
    plt.close(fig)
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
