#!/usr/bin/env python3
"""SV validation on real genome pairs: syn2bani `struct` calls vs dnadiff truth.

Reads, per pair directory under results/sv_validation/<pair>/:
  out.1coords   dnadiff 1-to-1 alignment coords (nucmer show-coords)
  out.qdiff     dnadiff query-side feature diff (GAP/DUP/INV/JMP/BRK)
  out.rdiff     dnadiff reference-side feature diff
  struct_sv.tsv syn2bani struct SV calls
  time.log      /usr/bin/time -v for dnadiff (wall clock)

Truth model:
  - Inversions: reverse-orientation blocks in out.1coords, clustered
    (merge when both ref and qry gaps to the previous block are <= MERGE_GAP);
    a cluster is a rearranged region, its span the inversion span.
  - Indels >= 1 kb: qdiff/rdiff GAP lines with |diff| >= 1000 (query-side
    positive diff = insertion in query) and DUP lines >= 1000 (qdiff DUP =
    extra query copy = insertion; rdiff DUP = extra reference copy = deletion
    in query). GAP events appear in both files with identical |diff|; the
    rdiff copy is dropped (multiset dedup on size).

Matching is span-based and one-to-many, because the two tools report at
different granularity: dnadiff fragments large accessory regions into many
GAP/BRK/DUP pieces, syn2bani calls the whole junction once. A call matches a
truth event of the same type when the truth position falls inside the call's
span on the appropriate genome (+/- MATCH_TOL). A truth event is covered when
any same-type call's span contains it.

Writes sv_inversion_compare.tsv, sv_indel_compare.tsv, sv_summary.tsv into
results/sv_validation/ and prints a human-readable summary.
"""

from pathlib import Path

import pandas as pd

SVDIR = Path(__file__).resolve().parent.parent / "results" / "sv_validation"
PAIRS = ["MG1655_vs_W3110", "MG1655_vs_Sakai", "CT18_vs_LT2"]
MERGE_GAP = 50_000
INDEL_MIN = 1_000
MATCH_TOL = 10_000

# Wall times (s): syn2bani struct = Mac Studio, min of 5 runs, 2026-08-14;
# dnadiff = HPC login node (MUMMER/3.23), /usr/bin/time -v elapsed.
RUNTIME = {
    # pair: (syn2bani_struct_s, dnadiff_s)
    "MG1655_vs_W3110": (0.048, 8.24),
    "MG1655_vs_Sakai": (0.054, 9.59),
    "CT18_vs_LT2": (0.050, 8.52),
}


def read_coords(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 13 or not f[0].strip().isdigit():
            continue
        rows.append(
            dict(
                r_start=int(f[0]),
                r_end=int(f[1]),
                q_start=int(f[2]),
                q_end=int(f[3]),
                ident=float(f[6]),
            )
        )
    df = pd.DataFrame(rows)
    df["reverse"] = df.q_start > df.q_end
    return df


def inversion_clusters(coords: pd.DataFrame) -> pd.DataFrame:
    rev = coords[coords.reverse].sort_values("r_start")
    clusters = []
    for _, b in rev.iterrows():
        qs, qe = min(b.q_start, b.q_end), max(b.q_start, b.q_end)
        if (
            clusters
            and b.r_start - clusters[-1]["r_end"] <= MERGE_GAP
            and qs - clusters[-1]["q_end"] <= MERGE_GAP
        ):
            c = clusters[-1]
            c["r_end"] = max(c["r_end"], b.r_end)
            c["q_start"] = min(c["q_start"], qs)
            c["q_end"] = max(c["q_end"], qe)
            c["n_blocks"] += 1
        else:
            clusters.append(
                dict(r_start=b.r_start, r_end=b.r_end, q_start=qs, q_end=qe, n_blocks=1)
            )
    df = pd.DataFrame(clusters)
    if len(df):
        df["q_size"] = df.q_end - df.q_start
    return df


def read_diff(path: Path, side: str) -> pd.DataFrame:
    """GAP: name GAP start end gap othergap diff -> indel, size |diff|.
    DUP: name DUP start end len -> copy-number event, size len.
    BRK: name BRK start end len -> unaligned stretch, size len (this is how
    dnadiff surfaces prophage/mobile regions it cannot align; a junction-level
    caller sees the same segment as one indel).
    `side` records which genome `pos` refers to."""
    rows = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        if f[1] == "GAP" and len(f) >= 7:
            diff = int(f[6])
            if abs(diff) < INDEL_MIN:
                continue
            ins_in_query = diff > 0 if side == "q" else diff < 0
            rows.append(
                dict(pos=(int(f[2]) + int(f[3])) // 2, size=abs(diff),
                     truth_type="Insertion" if ins_in_query else "Deletion",
                     side=side, src=f"{side}diff:GAP")
            )
        elif f[1] in ("DUP", "BRK"):
            size = int(f[4])
            if size < INDEL_MIN:
                continue
            rows.append(
                dict(pos=(int(f[2]) + int(f[3])) // 2, size=size,
                     truth_type="Insertion" if side == "q" else "Deletion",
                     side=side, src=f"{side}diff:{f[1]}")
            )
    return pd.DataFrame(rows)


def truth_events(pair_dir: Path) -> pd.DataFrame:
    q = read_diff(pair_dir / "out.qdiff", "q")
    r = read_diff(pair_dir / "out.rdiff", "r")
    # GAPs are reported from both sides with identical |diff|; drop the rdiff
    # copies one-for-one by size.
    q_gap_sizes = sorted(q[q.src.str.endswith("GAP")]["size"].tolist())
    keep = []
    for _, t in r.iterrows():
        if t.src.endswith("GAP") and t["size"] in q_gap_sizes:
            q_gap_sizes.remove(t["size"])  # consumed one qdiff copy
            continue
        keep.append(t)
    ev = pd.concat([q, pd.DataFrame(keep)], ignore_index=True)
    return ev.sort_values("pos", ignore_index=True)


def span_hit(lo: int, hi: int, pos: float, tol: int = MATCH_TOL) -> bool:
    return lo - tol <= pos <= hi + tol


def main() -> None:
    inv_rows, indel_rows, summ_rows = [], [], []
    for pair in PAIRS:
        d = SVDIR / pair
        coords = read_coords(d / "out.1coords")
        truth_inv = inversion_clusters(coords)
        tev = truth_events(d)

        sv = pd.read_csv(d / "struct_sv.tsv", sep="\t")
        sv_indel = sv[sv.sv_type.isin(["Insertion", "Deletion"])].copy()
        sv_inv = sv[sv.sv_type == "Inversion"].copy()
        sv_tra = sv[sv.sv_type == "Translocation"].copy()

        # --- inversions: overlap of query spans with reverse-block clusters ---
        for _, c in sv_inv.iterrows():
            best, best_ov = None, 0
            for ti, t in truth_inv.iterrows():
                ov = min(c.q_end, t.q_end) - max(c.q_start, t.q_start)
                if ov > best_ov:
                    best, best_ov = ti, ov
            row = dict(pair=pair, syn_q_start=c.q_start, syn_q_end=c.q_end,
                       syn_size=c["size"], matched=best is not None)
            if best is not None:
                t = truth_inv.loc[best]
                row.update(truth_q_start=t.q_start, truth_q_end=t.q_end,
                           truth_q_size=t.q_size, overlap_bp=best_ov,
                           start_err=c.q_start - t.q_start, end_err=c.q_end - t.q_end)
            inv_rows.append(row)

        # --- indels: span-based, one-to-many both ways ---
        covered = set()
        for _, c in sv_indel.iterrows():
            hits = []
            for ti, t in tev.iterrows():
                if t.truth_type != c.sv_type:
                    continue
                lo, hi = (c.q_start, c.q_end) if t.side == "q" else (c.r_start, c.r_end)
                if span_hit(lo, hi, t.pos):
                    hits.append(ti)
            row = dict(pair=pair, sv_type=c.sv_type,
                       syn_q_start=c.q_start, syn_q_end=c.q_end,
                       syn_r_start=c.r_start, syn_r_end=c.r_end,
                       syn_size=c["size"], n_truth_covered=len(hits),
                       matched=bool(hits))
            if hits:
                covered.update(hits)
                # size agreement against the single best (nearest) event
                best = min(
                    hits,
                    key=lambda ti: abs(
                        tev.loc[ti].pos
                        - ((c.q_start + c.q_end) / 2 if tev.loc[ti].side == "q"
                           else (c.r_start + c.r_end) / 2)
                    ),
                )
                t = tev.loc[best]
                row.update(truth_pos=t.pos, truth_size=t["size"],
                           size_ratio=c["size"] / t["size"], truth_src=t.src,
                           truth_size_sum=int(tev.loc[hits]["size"].sum()))
            indel_rows.append(row)
        for ti, t in tev.iterrows():
            if ti not in covered:
                indel_rows.append(dict(pair=pair, sv_type=t.truth_type, matched=False,
                                       truth_pos=t.pos, truth_size=t["size"],
                                       missed=True, truth_src=t.src))

        missed = [t for ti, t in tev.iterrows() if ti not in covered]
        syn_t, dna_t = RUNTIME[pair]
        summ_rows.append(dict(
            pair=pair,
            ani_note="",
            syn_inversions=len(sv_inv),
            syn_translocations=len(sv_tra),
            truth_inv_clusters=len(truth_inv),
            syn_indels=len(sv_indel),
            indel_tp=int(sum(1 for r in indel_rows
                             if r["pair"] == pair and r.get("matched") and not r.get("missed"))),
            indel_fp=int(sum(1 for r in indel_rows
                             if r["pair"] == pair and not r.get("matched") and not r.get("missed"))),
            truth_indels=len(tev),
            truth_covered=len(covered),
            truth_missed=len(missed),
            truth_missed_ge5kb=int(sum(1 for t in missed if t["size"] >= 5000)),
            syn_struct_s=syn_t,
            dnadiff_s=dna_t,
            speedup=round(dna_t / syn_t),
        ))

    inv_df = pd.DataFrame(inv_rows)
    indel_df = pd.DataFrame(indel_rows)
    summ = pd.DataFrame(summ_rows).drop(columns=["ani_note"])
    inv_df.to_csv(SVDIR / "sv_inversion_compare.tsv", sep="\t", index=False)
    indel_df.to_csv(SVDIR / "sv_indel_compare.tsv", sep="\t", index=False)
    summ.to_csv(SVDIR / "sv_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("=== summary ===")
    print(summ.to_string(index=False))
    print("\n=== inversions ===")
    print(inv_df.to_string(index=False))
    one2one = indel_df[(indel_df.matched == True) & (indel_df.n_truth_covered == 1)]  # noqa: E712
    print(f"\n=== 1:1 matched indels (n={len(one2one)}) ===")
    print(f"median size ratio {one2one.size_ratio.median():.3f}; "
          f"within [0.8, 1.25]: {(one2one.size_ratio.between(0.8, 1.25)).sum()}/{len(one2one)}")
    print("\n=== false-positive calls ===")
    print(indel_df[(indel_df.matched == False) & (indel_df.missed.isna())].to_string(index=False))  # noqa: E712
    print("\n=== missed truth events, by pair and size ===")
    mm = indel_df[indel_df.missed == True]  # noqa: E712
    print(mm.groupby(["pair", "sv_type"]).truth_size.describe()[["count", "50%", "max"]].to_string())
    print("\nlarge missed (>=10 kb):")
    print(mm[mm.truth_size >= 10_000].to_string(index=False))


if __name__ == "__main__":
    main()
