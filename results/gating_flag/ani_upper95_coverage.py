#!/usr/bin/env python3
"""Coverage of the ani_upper95 one-sided bound on the 105 BELOW_DETECTION
GTDB-ANIm pairs.

Replicates the Rust formula in Syn2bANI src/core/mle.rs exactly:

  p_upper   = Clopper-Pearson one-sided 95% upper limit on the found fraction
              (rule of three when found == 0)
  a_hi      = largest a with expected_retention(a) == p_upper
              (bisection on [0.5, 0.999999]; p_e(a) = sum_{m<=2} C(b,m)
              (1-a)^m a^(k-m) per enzyme stratum; this panel has no IUPAC-
              degenerate site positions)

Inputs:
  - per-pair chain-restricted strata: strata_2074/<QASM>__<RASM>.tsv
    (84 pairs that formed chains)
  - prerun TSV for the BELOW_DETECTION list and n_anchors (21 chains-empty
    pairs: found = n_anchors, an upper approximation of distinct anchored
    query tags; denominator = whole-genome in-panel query tag count, which is
    NOT available off-HPC, so those bounds are evaluated at N in {3k,5k,8k})
  - ANIm truth: gtdb_anim_joined.tsv column anim_ani

Output: ani_upper95_coverage.tsv + a summary on stdout.
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PRERUN = os.path.join(HERE, "anim_2074_prerun.tsv")
JOINED = os.path.join(HERE, "gtdb_anim_joined.tsv")
STRATA = os.path.join(HERE, "strata_2074")
ACC2SEQ = os.path.join(HERE, "..", "anim_2074_acc2seqid.tsv")

# Default panel (src/enzyme/registry.rs): name -> (tag_len, exact_site).
# All four anchors are fully specific A/C/G/T, so d2 = d3 = 0.
GEOM = {
    "BcgI": (32, 6),
    "AlfI": (32, 6),
    "AloI": (27, 7),
    "FalI": (27, 6),
}
TOL = 2
A_LO, A_HI = 0.50, 0.999_999
ALPHA = 0.05


def ln_binom(n, m):
    return math.lgamma(n + 1) - math.lgamma(m + 1) - math.lgamma(n - m + 1)


def p_found(a, tag_len, site_len, tol=TOL):
    b = tag_len - site_len
    return sum(
        math.exp(ln_binom(b, m) + m * math.log(1 - a) + (tag_len - m) * math.log(a))
        for m in range(tol + 1)
    )


def expected_retention(a, strata):
    # strata: list of (tag_len, site_len, total)
    tot = sum(n for _, _, n in strata)
    return sum(n * p_found(a, k, s) for k, s, n in strata) / tot


def binom_cdf_le(n, k, p):
    return sum(
        math.exp(
            ln_binom(n, i) + (i * math.log(p) if i else 0.0)
            + ((n - i) * math.log(1 - p) if i < n else 0.0)
        )
        for i in range(k + 1)
    )


def cp_upper(k, n, alpha=ALPHA):
    if n == 0:
        return float("nan")
    if k >= n:
        return 1.0
    if k == 0:
        return 1.0 - alpha ** (1.0 / n)
    lo, hi = k / n, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if binom_cdf_le(n, k, mid) > alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ani_upper_bound(strata, confidence=0.95):
    # strata: list of (tag_len, site_len, found, total)
    found = sum(f for _, _, f, _ in strata)
    total = sum(t for _, _, _, t in strata)
    if total == 0:
        return float("nan")
    p_u = cp_upper(found, total, 1.0 - confidence)
    if p_u >= 1.0:
        return 1.0
    curve = [(k, s, t) for k, s, _, t in strata]
    if p_u <= expected_retention(A_LO, curve):
        return A_LO
    if p_u >= expected_retention(A_HI, curve):
        return 1.0
    lo, hi = A_LO, A_HI
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if expected_retention(mid, curve) < p_u:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    seq2acc = {}
    with open(ACC2SEQ) as fh:
        for line in fh:
            acc, seq = line.strip().split("\t")[:2]
            seq2acc[seq] = acc

    anim = {}
    with open(JOINED) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            anim[(r["query"], r["reference"])] = float(r["anim_ani"])

    below = []
    with open(PRERUN) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["flag"] == "BELOW_DETECTION":
                below.append(r)
    print(f"BELOW_DETECTION pairs: {len(below)}")

    rows = []
    for r in below:
        q, ref = r["query"], r["reference"]
        truth = anim.get((q, ref))
        fn = os.path.join(STRATA, f"{seq2acc[q]}__{seq2acc[ref]}.tsv")
        strata = []
        with open(fn) as fh:
            for s in csv.DictReader(fh, delimiter="\t"):
                tag_len = int(s["tag_len"])
                body = int(s["body_len"])
                hist = [int(x) for x in s["hist"].split(",")]
                found = sum(hist)
                total = found + int(s["n_miss"])
                strata.append((tag_len, tag_len - body, found, total))
        if strata:
            bound = ani_upper_bound(strata)
            kind = "chained"
            found = sum(f for _, _, f, _ in strata)
            total = sum(t for _, _, _, t in strata)
        else:
            # Chains-empty: found = n_anchors (upper approx of distinct
            # anchored query tags), N = whole-genome query tag count —
            # unavailable off-HPC; evaluate a sensitivity range instead.
            bound = float("nan")
            kind = "empty"
            found = int(r["n_anchors"])
            total = 0
        rows.append(
            dict(query=q, reference=ref, kind=kind, found=found, total=total,
                 bound=100.0 * bound, truth=truth)
        )

    out = os.path.join(HERE, "ani_upper95_coverage.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["query", "reference", "kind", "found", "total",
                    "ani_upper95", "anim_truth", "covered"])
        for x in rows:
            cov = "" if x["truth"] is None or math.isnan(x["bound"]) else int(
                x["truth"] <= x["bound"] + 1e-9)
            w.writerow([x["query"], x["reference"], x["kind"], x["found"],
                        x["total"] or "",
                        "" if math.isnan(x["bound"]) else f'{x["bound"]:.4f}',
                        x["truth"], cov])

    chained = [x for x in rows if x["kind"] == "chained" and x["truth"] is not None]
    cov = sum(x["truth"] <= x["bound"] + 1e-9 for x in chained)
    print(f"\nchained BELOW_DETECTION pairs: {len(chained)}")
    print(f"  coverage (truth <= bound): {cov}/{len(chained)} = {cov/len(chained):.3f}")
    misses = [x for x in chained if x["truth"] > x["bound"] + 1e-9]
    for x in sorted(misses, key=lambda x: x["truth"] - x["bound"], reverse=True):
        print(f'  MISS {x["query"]} {x["reference"]}: truth {x["truth"]:.2f} '
              f'> bound {x["bound"]:.2f} (by {x["truth"]-x["bound"]:.2f})')
    bs = sorted(x["bound"] for x in chained)
    print(f"  bound range: {bs[0]:.2f}..{bs[-1]:.2f}, median {bs[len(bs)//2]:.2f}")

    empty = [x for x in rows if x["kind"] == "empty" and x["truth"] is not None]
    print(f"\nchains-empty pairs: {len(empty)} (found = n_anchors; N unknown, "
          f"sensitivity over panel-typical tag counts)")
    panel = [(32, 6), (32, 6), (27, 7), (27, 6)]  # BcgI, AlfI, AloI, FalI
    for n_tot in (3000, 5000, 8000):
        ok = 0
        for x in empty:
            # Spread N uniformly over the 4 panel enzymes.
            strata = [(k, s, x["found"] // 4, n_tot // 4) for k, s in panel]
            strata[0] = (32, 6, x["found"] - 3 * (x["found"] // 4), n_tot // 4)
            b = 100.0 * ani_upper_bound(strata)
            x.setdefault("bounds", {})[n_tot] = b
            ok += x["truth"] <= b + 1e-9
        print(f"  N={n_tot}: coverage {ok}/{len(empty)} = {ok/len(empty):.3f}")
    for x in sorted(empty, key=lambda x: -x["truth"]):
        b = x["bounds"]
        flag = "" if all(x["truth"] <= v for v in b.values()) else "  <- truth above bound"
        print(f'  {x["query"]} {x["reference"]}: truth {x["truth"]:.2f}, '
              f'anchors {x["found"]}, bound(N=3k/5k/8k) '
              f'{b[3000]:.2f}/{b[5000]:.2f}/{b[8000]:.2f}{flag}')


if __name__ == "__main__":
    main()
