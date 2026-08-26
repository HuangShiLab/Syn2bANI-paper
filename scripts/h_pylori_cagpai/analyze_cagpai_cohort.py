#!/usr/bin/env python3
"""analyze_cagpai_cohort.py — cagPAI-state classification + ANI/af discordance
analysis for the H. pylori cohort (Song et al. 2026, PNAS).

Inputs
------
--struct-dir   per-genome struct TSVs from s2_struct_vs_26695.slurm
               (<id>.vs_hp26695.struct.tsv, <id>.vs_hp26695_delpai.struct.tsv)
--triangle     all-vs-all edge list from s3_triangle_dist.slurm
--metadata     (optional) TSV with columns: accession, stage (NAG/AG/IM/GC)
--outdir       output directory

cagPAI locus in 26695 (NC_000915.1): 547,328-583,481 (36,154 bp).

Classification (from <id>.vs_hp26695.struct.tsv, query = strain):
  overlapping Deletion size at the cagPAI locus
    >= 30 kb  -> empty
    5-30 kb   -> partial
    < 5 kb    -> complete
"""
import argparse
import glob
import os
import sys
from collections import Counter

CAG_S, CAG_E = 547_328, 583_481
EMPTY_MIN = 30_000
PARTIAL_MIN = 5_000


def classify_strain(struct_tsv):
    """Sum Deletion sizes overlapping the cagPAI locus on the 26695 reference."""
    del_bp = 0
    events = []
    if not os.path.exists(struct_tsv):
        return None, 0, events
    with open(struct_tsv) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            sv_type, r_start, r_end, size = f[2], int(f[7]), int(f[8]), int(f[9])
            ov = max(0, min(r_end, CAG_E) - max(r_start, CAG_S))
            if sv_type == "Deletion" and ov > 0:
                del_bp += ov
                events.append((sv_type, r_start, r_end, size))
    if del_bp >= EMPTY_MIN:
        state = "empty"
    elif del_bp >= PARTIAL_MIN:
        state = "partial"
    else:
        state = "complete"
    return state, del_bp, events


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", required=True)
    ap.add_argument("--triangle", required=True)
    ap.add_argument("--metadata")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ani-bins", default="95,97,99")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # --- classify strains -------------------------------------------------
    states = {}
    for path in glob.glob(os.path.join(args.struct_dir, "*.vs_hp26695.struct.tsv")):
        sid = os.path.basename(path).split(".vs_hp26695")[0]
        state, del_bp, _ = classify_strain(path)
        if state:
            states[sid] = (state, del_bp)
    with open(os.path.join(args.outdir, "cagpai_states.tsv"), "w") as out:
        out.write("accession\tcagpai_state\tcagpai_deleted_bp\n")
        for sid, (st, bp) in sorted(states.items()):
            out.write(f"{sid}\t{st}\t{bp}\n")
    print("cagPAI states:", Counter(st for st, _ in states.values()))

    # --- metadata ---------------------------------------------------------
    stage = {}
    if args.metadata:
        with open(args.metadata) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ai, si = hdr.index("accession"), hdr.index("stage")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                stage[f[ai]] = f[si]

    # --- pairwise discordance --------------------------------------------
    bins = [float(x) for x in args.ani_bins.split(",")]
    bin_stat = {b: [0, 0] for b in bins}  # b -> [n_pairs, n_discordant]
    discordant_rows = []
    with open(args.triangle) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        qi, ri = header.index("query"), header.index("reference")
        ai_ = header.index("ani")
        aqi, ari = header.index("af_query"), header.index("af_reference")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, r = f[qi], f[ri]
            if q not in states or r not in states:
                continue
            ani = float(f[ai_])
            af_asym = abs(float(f[aqi]) - float(f[ari]))
            sq, sr = states[q][0], states[r][0]
            disc = sq != sr and "partial" not in (sq, sr)
            for b in bins:
                if ani >= b:
                    bin_stat[b][0] += 1
                    bin_stat[b][1] += disc
            if disc and ani >= bins[0]:
                discordant_rows.append(
                    (q, r, ani, af_asym, sq, sr,
                     stage.get(q, "NA"), stage.get(r, "NA")))
    with open(os.path.join(args.outdir, "discordant_pairs.tsv"), "w") as out:
        out.write("query\treference\tani\taf_asym\tcagpai_q\tcagpai_r\tstage_q\tstage_r\n")
        for row in sorted(discordant_rows, key=lambda x: -x[2]):
            out.write("\t".join(str(x) for x in row) + "\n")

    print("\nANI bin    n_pairs   n_discordant(complete vs empty)   rate")
    for b in bins:
        n, d = bin_stat[b]
        print(f">= {b:<5} {n:>9,} {d:>12,} {d / n * 100 if n else 0:>8.2f}%")

    if stage and discordant_rows:
        cross = sum(
            1 for *_, sq_st, sr_st in discordant_rows
            if sq_st != "NA" and sr_st != "NA"
            and (sq_st in ("IM", "GC")) != (sr_st in ("IM", "GC")))
        labelled = sum(1 for *_, a, b in discordant_rows if a != "NA" and b != "NA")
        print(f"\ndiscordant pairs crossing early/late stage boundary: "
              f"{cross}/{labelled} ({cross / labelled * 100 if labelled else 0:.1f}%)")

    # --- figures ----------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figures", file=sys.stderr)
        return

    pairs, conc, disc_pairs = [], [], []
    with open(args.triangle) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        qi, ri = header.index("query"), header.index("reference")
        ai_ = header.index("ani")
        aqi, ari = header.index("af_query"), header.index("af_reference")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, r = f[qi], f[ri]
            if q not in states or r not in states:
                continue
            ani = float(f[ai_])
            if ani < 90:
                continue
            af_asym = abs(float(f[aqi]) - float(f[ari]))
            row = (ani, af_asym)
            if states[q][0] != states[r][0] and "partial" not in (
                    states[q][0], states[r][0]):
                disc_pairs.append(row)
            else:
                conc.append(row)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    if conc:
        ax.scatter(*zip(*conc), s=3, alpha=0.3, c="grey", label="cagPAI-concordant")
    if disc_pairs:
        ax.scatter(*zip(*disc_pairs), s=8, alpha=0.7, c="crimson",
                   label="cagPAI-discordant")
    ax.set_xlabel("Syn2bANI ANI (%)")
    ax.set_ylabel("|ΔAF| (query vs reference)")
    ax.legend(markerscale=3, fontsize=8)
    ax.set_title("A  High-ANI pairs hide cagPAI presence/absence")

    ax = axes[1]
    xs = [f"≥{b:g}" for b in bins]
    rates = [bin_stat[b][1] / bin_stat[b][0] * 100 if bin_stat[b][0] else 0
             for b in bins]
    ax.bar(xs, rates, color="steelblue")
    ax.set_xlabel("ANI threshold (%)")
    ax.set_ylabel("cagPAI-discordant pairs (%)")
    ax.set_title("B  Discordance rate among ANI-search hits")

    ax = axes[2]
    if stage and discordant_rows:
        stage_pairs = Counter(
            tuple(sorted((r[6], r[7]))) for r in discordant_rows
            if r[6] != "NA" and r[7] != "NA")
        labels = [f"{a}–{b}" for a, b in stage_pairs]
        ax.barh(labels, list(stage_pairs.values()), color="darkorange")
        ax.set_xlabel("discordant pairs")
        ax.set_title("C  Stage composition of discordant pairs")
    else:
        ax.text(0.5, 0.5, "metadata not provided", ha="center", va="center")
        ax.set_title("C  Stage composition of discordant pairs")
    fig.tight_layout()
    figpath = os.path.join(args.outdir, "cagpai_discordance.pdf")
    fig.savefig(figpath)
    fig.savefig(figpath.replace(".pdf", ".png"), dpi=300)
    print(f"\nfigure: {figpath}")


if __name__ == "__main__":
    main()
