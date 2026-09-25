#!/usr/bin/env python3
"""ANI-identical but biologically discordant pairs in the 528-genome
H. pylori cohort (song 2026).

Crosses the skani all-vs-all ANI matrix with the merged genome table
(cagPAI structural state, cagA EPIYA genotype, vacA genotype, CARD
resistance genotype, disease group, FastBAPS lineage, SV burden) to find
pairs that ANI would rank as identical but that differ in structure-linked
biology.

Outputs (results/pairwise_ani/):
  high_ani_pairs_discordance.tsv   one row per ANI >= 97 pair
  discordance_rates_by_band.tsv    discordance rates per ANI band
  showcase_pairs.tsv               highest-ANI biology-discordant pairs
  ANI_DISCORDANT_BIOLOGY.md        summary for the manuscript
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "case_studies" / "h_pylori_cagpai" / "results"
OUT = RES / "pairwise_ani"

MERGED = RES / "metadata_assoc" / "merged_genome_table.tsv"
SKANI = OUT / "skani_triangle.tsv"

MIN_ANI = 97.0
BANDS = [(97, 99), (99, 99.5), (99.5, 99.9), (99.9, 100.01)]


def load_metadata():
    m = pd.read_csv(MERGED, sep="\t")
    m = m.rename(columns={"GCA_accession": "acc"})
    return m.set_index("acc")


def load_skani():
    s = pd.read_csv(SKANI, sep="\t")
    # skani triangle output columns: Ref_file Query_file ANI Align_fraction_* ...
    cols = {c.lower(): c for c in s.columns}
    s = s.rename(columns={cols["ref_file"]: "ref", cols["query_file"]: "qry",
                          cols["ani"]: "ani"})
    for side in ("ref", "qry"):
        s[side] = s[side].str.replace(r"\.fna$", "", regex=True)
        s[side] = s[side].str.split("/").str[-1]
    return s[["ref", "qry", "ani"]]


def main():
    meta = load_metadata()
    s = load_skani()
    s = s[s.ani >= MIN_ANI].copy()
    n_gen = len(meta)
    print(f"{len(s):,} pairs at ANI >= {MIN_ANI} across {n_gen} genomes")

    def flag_pair(row, col):
        a = meta.at[row.ref, col] if row.ref in meta.index else np.nan
        b = meta.at[row.qry, col] if row.qry in meta.index else np.nan
        if pd.isna(a) or pd.isna(b):
            return np.nan
        return 1 if str(a) != str(b) else 0

    for col, out in [
        ("status_extended", "d_cagpai_state"),
        ("cagA_genotype_grp", "d_cagA"),
        ("vacA_genotype", "d_vacA"),
        ("res_macrolide", "d_res_clarithromycin"),
        ("res_fluoroquinolone", "d_res_levofloxacin"),
        ("group", "d_disease_group"),
        ("fastbaps", "d_lineage"),
    ]:
        s[out] = s.apply(flag_pair, axis=1, args=(col,))

    # SV burden difference
    s["sv_burden_abs_diff"] = s.apply(
        lambda r: abs(meta.at[r.ref, "sv_total_n"] - meta.at[r.qry, "sv_total_n"])
        if r.ref in meta.index and r.qry in meta.index else np.nan, axis=1)

    s.to_csv(OUT / "high_ani_pairs_discordance.tsv", sep="\t", index=False)

    # ---- rates by band ------------------------------------------------------
    rows = []
    dcols = [c for c in s.columns if c.startswith("d_")]
    for lo, hi in BANDS:
        sub = s[(s.ani >= lo) & (s.ani < hi)]
        if len(sub) < 20:
            continue
        row = {"band": f"{lo:g}-{hi:g}", "n_pairs": len(sub)}
        for c in dcols:
            v = sub[c].dropna()
            row[c.replace("d_", "pct_diff_")] = 100 * v.mean() if len(v) else np.nan
        row["median_sv_burden_abs_diff"] = sub.sv_burden_abs_diff.median()
        rows.append(row)
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / "discordance_rates_by_band.tsv", sep="\t", index=False)

    # ---- showcase pairs ------------------------------------------------------
    # biology-relevant discordance at the highest ANI
    s["n_biology_discordant"] = s[["d_cagpai_state", "d_cagA", "d_vacA",
                                   "d_res_clarithromycin", "d_res_levofloxacin"]].sum(axis=1)
    show = (s[(s.n_biology_discordant > 0)]
            .sort_values(["ani", "n_biology_discordant"], ascending=[False, False])
            .head(40))
    show = show.assign(
        ref_cagpai=show.ref.map(meta.status_extended),
        qry_cagpai=show.qry.map(meta.status_extended),
        ref_cagA=show.ref.map(meta.cagA_genotype_grp),
        qry_cagA=show.qry.map(meta.cagA_genotype_grp),
        ref_vacA=show.ref.map(meta.vacA_genotype),
        qry_vacA=show.qry.map(meta.vacA_genotype),
        ref_group=show.ref.map(meta.group),
        qry_group=show.qry.map(meta.group),
        ref_lineage=show.ref.map(meta.fastbaps),
        qry_lineage=show.qry.map(meta.fastbaps),
    )
    show.to_csv(OUT / "showcase_pairs.tsv", sep="\t", index=False)

    # ---- summary -------------------------------------------------------------
    ge99 = s[s.ani >= 99]
    ge999 = s[s.ani >= 99.9]
    def pct(sub, c):
        v = sub[c].dropna()
        return f"{100 * v.mean():.1f}% (n={len(v):,})" if len(v) else "n/a"

    lines = [
        "# ANI-identical but biologically discordant pairs (H. pylori, n=528)",
        "",
        f"Pairs at ANI >= 97: {len(s):,}; >= 99: {len(ge99):,}; >= 99.9: {len(ge999):,}",
        "",
        "## Discordance rates",
        "",
        "| band | n | cagPAI state | cagA genotype | vacA | clarithro-R | levo-R | disease group | same lineage |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in rates.iterrows():
        lines.append(
            f"| {r.band} | {int(r.n_pairs):,} | {r.pct_diff_cagpai_state:.1f}% | "
            f"{r.pct_diff_cagA:.1f}% | {r.pct_diff_vacA:.1f}% | "
            f"{r.pct_diff_res_clarithromycin:.1f}% | {r.pct_diff_res_levofloxacin:.1f}% | "
            f"{r.pct_diff_disease_group:.1f}% | {100 - r.pct_diff_lineage:.1f}% |")
    lines += [
        "",
        f"At ANI >= 99.9 (near-clonal): cagPAI state differs in {pct(ge999, 'd_cagpai_state')}",
        f"pairs; cagA genotype in {pct(ge999, 'd_cagA')}; vacA in {pct(ge999, 'd_vacA')};",
        f"clarithromycin resistance genotype in {pct(ge999, 'd_res_clarithromycin')}.",
        "",
        "Showcase pairs: `showcase_pairs.tsv` (sorted by ANI, then number of",
        "discordant biology columns).",
        "",
        "Files: high_ani_pairs_discordance.tsv, discordance_rates_by_band.tsv,",
        "showcase_pairs.tsv.",
    ]
    (OUT / "ANI_DISCORDANT_BIOLOGY.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
