#!/usr/bin/env python3
"""Artifact-threshold sensitivity and reverse-complement diagnostic for H. pylori cagPAI.

Uses the *unfiltered* syn2bani struct BEDs (struct_vs_26695/) and recomputes
cagPAI extended states under a range of artifact-span thresholds.  Also checks
whether residual large structural calls are enriched for inversions, which would
suggest whole-genome reverse-complement coordinate conventions rather than
circular-origin rotation.
"""
import csv
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chi2 as chi2_dist

ROOT = Path(__file__).resolve().parent.parent
STATES = ROOT / "results" / "cagpai_states.tsv"
META = Path("/Volumes/MoneyCat/Data/song_2026_hpylori/metadata.csv")
UNFILTERED_DIR = Path("/Volumes/MoneyCat/Data/song_2026_hpylori/struct_vs_26695")
OUT_SENS = ROOT / "results" / "cagpai_artifact_threshold_sensitivity.tsv"
OUT_RC = ROOT / "results" / "cagpai_reverse_complement_diagnostic.tsv"
PLOT_SCRIPT = ROOT / "scripts" / "plot_cagpai_sensitivity.py"

REF_LENGTH = 1_667_825
CAG_START = 547_327
CAG_END = 583_481
CAG_BUFFER = 2_000
REGION = (CAG_START - CAG_BUFFER, CAG_END + CAG_BUFFER)
THRESHOLDS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]


def overlaps_cagpai(chrom, start, end):
    return chrom == "NC_000915.1" and start < REGION[1] and end > REGION[0]


def parse_bed(path):
    """Return list of all SV dicts from an unfiltered BED."""
    svs = []
    if not path.exists():
        return svs
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
            sv_type = name.split("_")[0]
            svs.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "type": sv_type,
                "span": end - start,
                "span_fraction": (end - start) / REF_LENGTH,
            })
    return svs


def classify_with_threshold(status, all_svs, threshold):
    """Return (extended_state, retained_svs, excluded_svs) for a given threshold."""
    excluded = [s for s in all_svs if s["span"] > threshold * REF_LENGTH]
    retained = [s for s in all_svs if s["span"] <= threshold * REF_LENGTH and overlaps_cagpai(s["chrom"], s["start"], s["end"])]

    if status in ("partial", "empty"):
        return status, retained, excluded

    # status == complete
    rearr = [s for s in retained if s["type"] in ("INV", "TRA")]
    if rearr:
        return "complete_rearranged", rearr, excluded
    large_del = [s for s in retained if s["type"] == "DEL" and s["span"] >= 10_000]
    if large_del:
        return "partial", large_del, excluded
    return "complete_collinear", [], excluded


def build_table_dict(df, case_col, outcome_col, stratum_col="fastbaps"):
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


def cmh(table_dict):
    """Return CMH stat, p, and MH OR."""
    tables = np.array([table_dict[s] for s in table_dict], dtype=float)
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
    return float(cmh_stat), float(cmh_p), float(or_mh)


def crude_or(table_dict):
    crude = sum(np.array(t) for t in table_dict.values())
    a, b, c, d = crude[0, 0], crude[0, 1], crude[1, 0], crude[1, 1]
    or_val = (a * d) / (b * c) if b * c > 0 else np.inf
    chi2_val, p_val, _, _ = chi2_contingency(crude, correction=False)
    return int(a + b + c + d), float(or_val), float(p_val)


def run_sensitivity(states_df, meta_df):
    rows = []
    for threshold in THRESHOLDS:
        ext_rows = []
        for _, r in states_df.iterrows():
            gid = r["genome"]
            bed = UNFILTERED_DIR / f"{gid}.vs_hp26695.bed"
            all_svs = parse_bed(bed)
            ext_state, retained, excluded = classify_with_threshold(r["status"], all_svs, threshold)
            ext_rows.append({
                "genome": gid,
                "status_extended": ext_state,
                "n_excluded": len(excluded),
                "n_retained": len(retained),
            })
        ext_df = pd.DataFrame(ext_rows)
        df = states_df.merge(meta_df, left_on="genome", right_on="assembly", how="left")
        df = df.merge(ext_df[["genome", "status_extended", "n_excluded", "n_retained"]], on="genome")
        df["present"] = df["status_extended"].isin(["complete_collinear", "complete_rearranged"])
        df["rearranged"] = df["status_extended"] == "complete_rearranged"
        df["case_gc"] = df["group"] == "GC"
        df["case_adv"] = df["group"].isin(["GC", "IM"])

        counts = df["status_extended"].value_counts().to_dict()
        n_rearr = counts.get("complete_rearranged", 0)
        n_coll = counts.get("complete_collinear", 0)
        n_excluded_total = df["n_excluded"].sum()
        n_flagged = (df["n_excluded"] > 0).sum()

        # presence GC vs NAG
        sub = df[df["group"].isin(["GC", "NAG"])]
        td = build_table_dict(sub, "case_gc", "present")
        n_pres, or_pres, p_pres = crude_or(td)
        _, p_pres_cmh, or_pres_mh = cmh(td)

        # rearrangement GC vs NAG (among present)
        sub = df[df["group"].isin(["GC", "NAG"]) & df["present"]]
        td = build_table_dict(sub, "case_gc", "rearranged")
        n_rearr_gc, or_rearr, p_rearr = crude_or(td)
        _, p_rearr_cmh, or_rearr_mh = cmh(td)

        # presence advanced vs early
        sub = df[df["group"].isin(["GC", "IM", "AG", "NAG"])]
        td = build_table_dict(sub, "case_adv", "present")
        _, or_pres_adv, p_pres_adv = crude_or(td)
        _, p_pres_adv_cmh, or_pres_adv_mh = cmh(td)

        # rearrangement advanced vs early
        sub = df[df["group"].isin(["GC", "IM", "AG", "NAG"]) & df["present"]]
        td = build_table_dict(sub, "case_adv", "rearranged")
        _, or_rearr_adv, p_rearr_adv = crude_or(td)
        _, p_rearr_adv_cmh, or_rearr_adv_mh = cmh(td)

        rows.append({
            "threshold": threshold,
            "n_complete_rearranged": n_rearr,
            "n_complete_collinear": n_coll,
            "n_excluded_calls": n_excluded_total,
            "n_flagged_genomes": n_flagged,
            "presence_gc_n": n_pres,
            "presence_gc_crude_OR": or_pres,
            "presence_gc_crude_p": p_pres,
            "presence_gc_CMHR": or_pres_mh,
            "presence_gc_CMHR_p": p_pres_cmh,
            "rearr_gc_n": n_rearr_gc,
            "rearr_gc_crude_OR": or_rearr,
            "rearr_gc_crude_p": p_rearr,
            "rearr_gc_CMHR": or_rearr_mh,
            "rearr_gc_CMHR_p": p_rearr_cmh,
            "presence_adv_crude_OR": or_pres_adv,
            "presence_adv_crude_p": p_pres_adv,
            "presence_adv_CMHR": or_pres_adv_mh,
            "presence_adv_CMHR_p": p_pres_adv_cmh,
            "rearr_adv_crude_OR": or_rearr_adv,
            "rearr_adv_crude_p": p_rearr_adv,
            "rearr_adv_CMHR": or_rearr_adv_mh,
            "rearr_adv_CMHR_p": p_rearr_adv_cmh,
        })
    return pd.DataFrame(rows)


def reverse_complement_diagnostic(states_df):
    """Check whether large unfiltered calls are inversions (RC convention) or translocations (origin rotation)."""
    rows = []
    for _, r in states_df.iterrows():
        gid = r["genome"]
        bed = UNFILTERED_DIR / f"{gid}.vs_hp26695.bed"
        svs = parse_bed(bed)
        large = [s for s in svs if s["span_fraction"] >= 0.20]
        if not large:
            continue
        largest = max(large, key=lambda s: s["span"])
        rows.append({
            "genome": gid,
            "status": r["status"],
            "n_large_calls": len(large),
            "largest_span": largest["span"],
            "largest_span_fraction": largest["span_fraction"],
            "largest_type": largest["type"],
            "largest_coords": f"{largest['chrom']}:{largest['start']}-{largest['end']}",
            "large_inv_count": sum(1 for s in large if s["type"] == "INV"),
            "large_tra_count": sum(1 for s in large if s["type"] == "TRA"),
            "large_del_count": sum(1 for s in large if s["type"] == "DEL"),
            "large_ins_count": sum(1 for s in large if s["type"] == "INS"),
        })
    return pd.DataFrame(rows)


def write_plot_script():
    """Emit a small matplotlib script to visualise the sensitivity curve."""
    if PLOT_SCRIPT.exists():
        return
    PLOT_SCRIPT.write_text("""#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sens = pd.read_csv(ROOT / "results" / "cagpai_artifact_threshold_sensitivity.tsv", sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)

ax = axes[0]
ax.plot(sens["threshold"], sens["n_complete_rearranged"], "o-", label="complete_rearranged")
ax.plot(sens["threshold"], sens["n_excluded_calls"], "s--", label="excluded calls")
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("Count")
ax.set_title("cagPAI rearrangement count vs. artifact threshold")
ax.legend()

ax = axes[1]
import math
pvals = sens["rearr_gc_crude_p"].replace(0, 1e-300)
neglogp = [-math.log10(p) for p in pvals]
ax.plot(sens["threshold"], neglogp, "o-", label="rearrangement GC vs NAG")
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("-log10(p) crude")
ax.set_title("Association p-value sensitivity")
ax.axhline(-math.log10(0.05), color="red", linestyle="--", linewidth=0.8)
ax.legend()

fig.savefig(ROOT / "results" / "cagpai_artifact_threshold_sensitivity.png", dpi=300)
print("Wrote", ROOT / "results" / "cagpai_artifact_threshold_sensitivity.png")
""")
    PLOT_SCRIPT.chmod(0o755)


def main():
    states = pd.read_csv(STATES, sep="\t")
    meta = pd.read_csv(META)
    meta["group"] = meta["group"].astype(str).str.strip()

    print("Running artifact-threshold sensitivity...", file=sys.stderr)
    sens_df = run_sensitivity(states, meta)
    sens_df.to_csv(OUT_SENS, sep="\t", index=False)
    print(sens_df.to_string(index=False))
    print(f"\nWrote {OUT_SENS}", file=sys.stderr)

    print("\nRunning reverse-complement diagnostic...", file=sys.stderr)
    rc_df = reverse_complement_diagnostic(states)
    rc_df.to_csv(OUT_RC, sep="\t", index=False)
    print(rc_df.to_string(index=False))
    print(f"\nWrote {OUT_RC}", file=sys.stderr)

    write_plot_script()
    print(f"\nPlot script: {PLOT_SCRIPT}", file=sys.stderr)


if __name__ == "__main__":
    main()
