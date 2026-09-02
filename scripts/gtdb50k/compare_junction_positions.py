#!/usr/bin/env python3
"""Do Syn2b and dnadiff put rearrangement breakpoints in the same places?

`MATH_REVIEW.md` claims they do, at a ~2790x cost difference, on the strength of two
hand-checked coordinates. This turns that into a statistic over every gtdb50k pair,
from the table `collect_junction_coordinates.py` writes.

MATCHING
--------
Nearest-neighbour in both directions double-counts: one Syn2b junction can be the
nearest to three dnadiff boundaries and look like three successes. So the two sets
are matched **one-to-one**, greedily by increasing distance, and everything left
over is reported as unmatched on the side it came from:

    matched          both methods placed a boundary, distance = |syn2b - dnadiff|
    syn2b only       Syn2b placed one where dnadiff did not
    dnadiff only     dnadiff placed one where Syn2b did not

Neither leftover is automatically an error. dnadiff aligns; Syn2b reads landmark
order. An event smaller than landmark spacing is invisible to Syn2b by construction
(the detection floor, not a bug), and dnadiff's 1-to-1 filter drops alignments Syn2b
still has landmarks in. The point of the table is the *distance distribution of the
matched set* -- that is what "the same coordinates" has to mean.

WHAT LIMITS THE DISTANCE
------------------------
Syn2b reports the position of the left landmark of the broken adjacency, so its
error is bounded by landmark spacing, not by anything about the algorithm. Measured
on a closed E. coli K-12 control with three known 200 kb inversions:

    BcgI            2,872 landmarks   spacing 1,582 bp   median err   614 bp
    4-enzyme        6,079 landmarks   spacing   747 bp   median err   273 bp
    FracMinHash 750 6,034 landmarks   spacing   753 bp   median err   456 bp
    FracMinHash 200 22,708 landmarks  spacing   200 bp   median err   248 bp
    FracMinHash  50 90,394 landmarks  spacing    50 bp   median err    44 bp

So a median matched distance near the panel's landmark spacing is the *expected*
result and the correct claim to make. A distance much larger than spacing is the
finding that would need explaining.

USAGE
-----
    python3 compare_junction_positions.py \\
        --coords   $WORK/junction_coords/junction_coordinates.tsv \\
        --truth    $WORK/high_ani_truth.tsv \\
        --syn2b    $WORK/syn2b_inverted_fraction_50k.tsv \\
        --out      $WORK/junction_coords/position_agreement.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

THRESHOLDS = [500, 1_000, 2_000, 5_000, 10_000, 50_000]


def parse_positions(s):
    if not s:
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if tok:
            try:
                out.append(int(tok))
            except ValueError:
                pass
    return sorted(out)


def match_one_to_one(a, b):
    """Greedy one-to-one match by increasing distance.

    Returns (pairs, a_only, b_only) where pairs is a list of (a_pos, b_pos, dist).
    Quadratic in the number of junctions per pair, which is fine: the gtdb50k
    maximum is in the hundreds, not the thousands.

    Greedy is not the minimum-total-distance assignment -- on a crossing pair like
    a = [1000, 2000], b = [1900, 2100] it takes (2000, 1900) first and is then forced
    into (1000, 2100), for 1200 against the optimal 1000. It is used anyway because
    the bias runs the safe way: greedy can only *inflate* matched distances relative
    to the optimal assignment, never shrink them, so the agreement reported here is a
    lower bound. Crossings also require two events closer together than the spread of
    the two methods' estimates, where the pairing is ambiguous in any case.
    """
    if not a or not b:
        return [], list(a), list(b)
    cand = sorted(
        ((abs(x - y), i, j) for i, x in enumerate(a) for j, y in enumerate(b))
    )
    used_a, used_b, pairs = set(), set(), []
    for d, i, j in cand:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        pairs.append((a[i], b[j], d))
    return (pairs,
            [x for i, x in enumerate(a) if i not in used_a],
            [y for j, y in enumerate(b) if j not in used_b])


def median(v):
    v = sorted(v)
    if not v:
        return float("nan")
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def pct(v, q):
    v = sorted(v)
    if not v:
        return float("nan")
    k = min(len(v) - 1, int(round(q * (len(v) - 1))))
    return v[k]


def load_lookup(path, key, cols):
    """Read `cols` from a TSV keyed on `key`. Missing file is not an error."""
    if not path or not Path(path).exists():
        return {}
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = r.get(key)
            if k and k not in out:
                out[k] = {c: r.get(c) for c in cols if c in r}
    return out


def fnum(d, k):
    try:
        return float(d.get(k, ""))
    except (TypeError, ValueError):
        return float("nan")


def band_of(ani):
    if ani != ani:
        return "unknown"
    for lo, hi in [(0, 90), (90, 95), (95, 97), (97, 99), (99, 100.01)]:
        if lo <= ani < hi:
            return f"{lo}-{hi if hi <= 100 else 100}"
    return "unknown"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--coords", required=True,
                   help="junction_coordinates.tsv from collect_junction_coordinates.py")
    p.add_argument("--truth", default=None,
                   help="optional TSV with pairid + anim_ani, for banding")
    p.add_argument("--syn2b", default=None,
                   help="optional run_syn2b_inverted_fraction.py output, for "
                        "syn2b_shared_tags and syn2b_observable_fraction -- landmark "
                        "spacing and contig count, which bound the distance")
    p.add_argument("--genome-length", type=float, default=4.5e6,
                   help="assumed genome length for the spacing estimate, bp")
    p.add_argument("--out", required=True)
    p.add_argument("--require-count-match", action="store_true", default=True,
                   help="use only pairs where the derived dnadiff count matches "
                        "dd.report exactly (default on; the derivation is not "
                        "dnadiff's own output)")
    p.add_argument("--all-pairs", dest="require_count_match", action="store_false",
                   help="include pairs that failed the derivation check")
    args = p.parse_args()

    truth = load_lookup(args.truth, "pairid", ["anim_ani"])
    synt = load_lookup(args.syn2b, "pairid",
                       ["syn2b_shared_tags", "syn2b_observable_fraction"])

    rows, skipped_check, no_data, bad_frame = [], 0, 0, 0
    with open(args.coords) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["syn2b_n"] in ("-1", "") or r["dnadiff_n"] in ("-1", ""):
                no_data += 1
                continue
            # Only frames the collector could put in Syn2b's genome-wide coordinates.
            # A multi-contig reference whose TGT was unavailable yields dnadiff
            # positions counted from the start of each contig, which would differ
            # from Syn2b's by the cumulative length of every preceding one.
            if r.get("frame") not in (None, "", "single_contig", "genome_wide"):
                bad_frame += 1
                continue
            if args.require_count_match and r.get("count_diff") != "0":
                skipped_check += 1
                continue
            a = parse_positions(r["syn2b_pos"])
            b = parse_positions(r["dnadiff_pos"])
            if not a and not b:
                continue
            pairs, a_only, b_only = match_one_to_one(a, b)
            t = truth.get(r["pairid"], {})
            s = synt.get(r["pairid"], {})
            m = fnum(s, "syn2b_shared_tags")
            rows.append({
                "pairid": r["pairid"],
                "anim_ani": t.get("anim_ani", ""),
                "band": band_of(fnum(t, "anim_ani")),
                "n_syn2b": len(a),
                "n_dnadiff": len(b),
                "n_matched": len(pairs),
                "n_syn2b_only": len(a_only),
                "n_dnadiff_only": len(b_only),
                "median_dist": f"{median([d for _, _, d in pairs]):.0f}" if pairs else "",
                "max_dist": max((d for _, _, d in pairs), default=""),
                "shared_tags": s.get("syn2b_shared_tags", ""),
                "spacing_bp": f"{args.genome_length / m:.0f}" if m == m and m > 0 else "",
                "observable_fraction": s.get("syn2b_observable_fraction", ""),
                "dists": ",".join(str(d) for _, _, d in pairs),
            })

    if not rows:
        print("No usable pairs. Check --coords, and whether the derivation check "
              "passed for any pair (rerun with --all-pairs to see).", file=sys.stderr)
        return 2

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [k for k in rows[0] if k != "dists"] + ["dists"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}  ({len(rows)} pairs)")
    if skipped_check:
        print(f"  excluded {skipped_check} pairs whose derived dnadiff count "
              f"disagreed with dd.report (use --all-pairs to include them)")
    if no_data:
        print(f"  excluded {no_data} pairs missing one side entirely")
    if bad_frame:
        print(f"  excluded {bad_frame} pairs whose dnadiff coordinates could not be "
              f"placed in Syn2b's genome-wide frame (multi-contig reference, no TGT)")

    # ---- the headline: distance distribution of the matched set --------------
    all_d = [d for r in rows for d in parse_positions(r["dists"])]
    tot_syn = sum(r["n_syn2b"] for r in rows)
    tot_dd = sum(r["n_dnadiff"] for r in rows)
    tot_m = sum(r["n_matched"] for r in rows)
    print(f"\nboundaries: Syn2b {tot_syn:,}, dnadiff {tot_dd:,}, matched one-to-one {tot_m:,}")
    if tot_syn:
        print(f"  Syn2b boundaries with a dnadiff partner   {100*tot_m/tot_syn:5.1f}%")
    if tot_dd:
        print(f"  dnadiff boundaries with a Syn2b partner   {100*tot_m/tot_dd:5.1f}%")

    if all_d:
        print(f"\nmatched-pair distance, n = {len(all_d):,}")
        # median() and pct() disagree at q=0.5 on even n; use median() there so the
        # headline and the per-band table below report the same number.
        print(f"  {'median':<8} {median(all_d):>10,.0f} bp")
        for label, q in [("p75", 0.75), ("p90", 0.90), ("p99", 0.99)]:
            print(f"  {label:<8} {pct(all_d, q):>10,.0f} bp")
        print()
        for t in THRESHOLDS:
            n = sum(1 for d in all_d if d <= t)
            print(f"  within {t:>7,} bp   {n:>8,}  ({100*n/len(all_d):5.1f}%)")

        spacings = [float(r["spacing_bp"]) for r in rows if r["spacing_bp"]]
        if spacings:
            print(f"\nmedian landmark spacing across these pairs: {median(spacings):,.0f} bp")
            print("Syn2b's coordinate error is bounded by that spacing, so a matched "
                  "distance of the same order is the expected result, not a weak one.")

    # ---- banded, because agreement should not be uniform in divergence -------
    print(f"\n{'band':<10} {'pairs':>7} {'syn2b':>8} {'dnadiff':>8} {'matched':>8} "
          f"{'med dist':>9} {'<=5kb':>7}")
    bands = {}
    for r in rows:
        bands.setdefault(r["band"], []).append(r)
    for band in sorted(bands):
        rs = bands[band]
        d = [x for r in rs for x in parse_positions(r["dists"])]
        s_ = sum(r["n_syn2b"] for r in rs)
        dd_ = sum(r["n_dnadiff"] for r in rs)
        m_ = sum(r["n_matched"] for r in rs)
        w5 = 100 * sum(1 for x in d if x <= 5000) / len(d) if d else float("nan")
        print(f"{band:<10} {len(rs):>7,} {s_:>8,} {dd_:>8,} {m_:>8,} "
              f"{median(d):>9,.0f} {w5:>6.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
