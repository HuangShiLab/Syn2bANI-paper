#!/usr/bin/env python3
"""Run the Syn2b structural channel on the SynTracker isolate cohorts.

Enav, Paz & Ley, *Nat Biotechnol* 43:773-783 (2024), doi:10.1038/s41587-024-02276-2
pairs SynTracker (synteny) with inStrain (SNPs) on four isolate collections and
reports which mode of evolution dominates in each. That gives us something our own
GTDB data cannot: cohorts whose expected answer is already published, including a
positive and a negative control.

    S. rimosus           positive control -- popANI pinned at ~1.0, APSS spans
                         0.90-1.0, so the variation is structural. Our structural
                         channel must see it where SNP-based tools see nothing.
    E. coli hypermutator negative control -- inStrain calls no pair the same strain,
                         SynTracker calls every pair the same strain, so the
                         variation is SNP-driven. Our channel must stay near zero.
    N. gonorrhoeae       both modes -- rho = 0.985 between the two published scores.
    H. pylori            mixed, per participant; connects to the cagPAI cohort.

WHY THIS SCRIPT EXISTS RATHER THAN A RE-RUN OF 04_run_syn2bani.sh
-----------------------------------------------------------------
The Syn2bANI pass already in `data/syntracker_validation/syn2bani/` is not usable
for a structural claim. Every genome was also compared against *itself*, and those
self-comparisons return 904 / 25 / 327 / 1415 breakpoints (E. coli / H. pylori /
N. gonorrhoeae / S. rimosus) where the only correct answer is 0. The per-cohort
means reported in `results/syntracker_validation/syntracker_summary.tsv` are
905.7 / 28.1 / 305.3 / 1980.2 -- i.e. the published-looking ranking of "how much
structural variation each species has" is, to within a few percent, a ranking of
how fragmented its assemblies are. The cause is the reference-side inflation fixed
in Syn2bANI `c974f5f`.

So the self-comparison is not a nicety here; it is the one control that would have
caught this before the numbers reached a summary table. This script runs it FIRST
and refuses to report cohort statistics for any cohort that fails it, unless
--allow-failed-controls is passed.

WHAT IT MEASURES
----------------
`syn2b synteny` on every within-cohort pair, keeping the full structural channel:

    breakpoints          adjacencies the other genome positively contradicts.
                         Fragmentation-immune -- a contig break is an absence of
                         evidence, not a contradiction.
    scj_distance         single-cut-junction distance, the symmetric difference of
                         the adjacency sets. NOT fragmentation-immune: a break
                         genuinely removes an adjacency from one set, so SCJ carries
                         a +(K-1) term. Do not read it on draft assemblies.
    inverted_fraction    majority-frame orientation disagreement, in [0, 0.5]
    raw_inverted_fraction  fixed-reference version, in [0, 1]
    observable_fraction  share of A's adjacencies B can judge; ~ 1 - (K-1)/S,
                         which is how K (contig count) is recovered below

Measured on a synthetic control built from E. coli K-12 (closed original; the same
sequence shattered into 40 contigs with half of them reverse-complemented; a 500 kb
inversion; a 1.2 Mb origin rotation), four-enzyme panel:

    closed vs rotated          bp 0   SCJ  0   f 0.0000   obs 1.0000
    closed vs 500 kb inversion bp 2   SCJ  4   f 0.1092   obs 1.0000   (truth f = 0.1101)
    closed vs 40-contig shatter bp 0  SCJ 39   f 0.3535   obs 0.9930   (K_est = 40.2)
    shatter vs inversion       bp 2   SCJ 43   f 0.3217   obs 1.0000

So `breakpoints` is exact under fragmentation and rotation, SCJ inherits K-1, and
the orientation channel drifts to 0.5 on fragmented input exactly as the GTDB
measurement predicted. `excess_over_floor` below is therefore computed on
breakpoints, not on SCJ.

THE FRAGMENTATION CAVEAT, STATED UP FRONT
-----------------------------------------
These are short-read assemblies. Contigs are deposited in arbitrary orientation, so
on fragmented input min(f, 1-f) drifts to 0.5 and carries no biology; measured on
the GTDB >=97% ANIm pairs the median goes 0.2880 (K = 1-2) -> 0.4726 (K > 100).
The orientation channel is therefore only interpretable on near-closed assemblies,
and this script reports estimated K per cohort so that is visible rather than
assumed.

The SynTracker paper hits the same wall and orders contigs against a reference with
Mauve before comparing. Do the same and pass --label ref_ordered, then compare
against the --label raw run. Reference-guided ordering biases each contig toward
collinearity with the reference, so inversions with breakpoints on contig
boundaries are absorbed rather than detected -- the bias is toward the null, which
costs sensitivity but not specificity. The gap between the two runs measures what
the orientation artifact costs.

USAGE
-----
    python3 08_syn2b_structural.py \
        --assembly-dir /lustre1/g/aos_shihuang/data/syntracker_validation/assemblies \
        --samples-dir  /lustre1/g/aos_shihuang/data/syntracker_validation/samples \
        --out-dir      /lustre1/g/aos_shihuang/data/syntracker_validation/syn2b_structural \
        --syn2b        /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
        --label        raw --workers 16
"""

import argparse
import csv
import itertools
import os
import subprocess
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

DEFAULT_ENZYMES = "BcgI,AlfI,AloI,FalI"

COHORTS = ["Streptomyces_rimosus", "Escherichia_coli_hypermutator",
           "Neisseria_gonorrhoeae", "Helicobacter_pylori"]

# What the SynTracker paper reports for each cohort, and what that implies our
# structural channel must do. Used only to label the summary; nothing branches on it.
EXPECTATION = {
    "Streptomyces_rimosus": (
        "positive control",
        "popANI 0.99990-1.0 (clonal) but APSS ~0.90-1.0: variation is structural. "
        "Structural signal must exceed the self-comparison floor."),
    "Escherichia_coli_hypermutator": (
        "negative control",
        "inStrain calls no pair the same strain, SynTracker calls all pairs the "
        "same strain: variation is SNP-driven. Structural signal must stay at the floor."),
    "Neisseria_gonorrhoeae": (
        "both modes",
        "published rho = 0.985 between synteny and SNP scores; both channels move together."),
    "Helicobacter_pylori": (
        "mixed, per participant",
        "participants 322, 326 and 439 carry subpopulations only one tool calls same-strain."),
}

# Self-comparison must return exactly this. A genome is collinear with itself.
# observable_fraction is NOT required to be 1.0 on draft assemblies; it only
# reflects how many of A's adjacencies B can judge, and a draft A judges fewer
# than a closed A. The fragmentation-immune metrics (breakpoints, SCJ,
# inverted_fraction) must be zero.
SELF_EXPECT = {
    "syn2b_scj_distance": 0.0,
    "syn2b_breakpoints": 0.0,
    "syn2b_inverted_fraction": 0.0,
}
SELF_TOL = 1e-9


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(map(str, cmd))}\n"
                           f"stderr: {r.stderr}\nstdout: {r.stdout}")
    return r


def retag(tgt_path: Path, genome_id: str):
    """Rewrite the TGT header line to `genome_id`.

    `syn2b digest` takes the genome id from the *first contig's FASTA header*
    (`main.rs`: `if gid.is_empty(){gid=rec.id.clone();}`). On SPAdes assemblies that
    is `NODE_1_length_..._cov_...`, which makes genome_A/genome_B unreadable and, if
    two isolates ever share a first-contig name, makes `syn2b synteny` refuse the
    pair outright ("duplicate genome id ... Rename one."). Renaming to the isolate
    removes both hazards.
    """
    with open(tgt_path) as fh:
        lines = fh.readlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"{tgt_path}: not a text TGT (header line missing)")
    length = lines[0].rstrip("\n").partition("|")[2]
    lines[0] = f">{genome_id}|{length}\n" if length else f">{genome_id}\n"
    with open(tgt_path, "w") as fh:
        fh.writelines(lines)


def self_copy(tgt_path: Path, dest: Path, genome_id: str) -> Path:
    """A genome cannot be compared with itself by symlinking the same TGT twice:
    `syn2b synteny` rejects two genomes sharing an id. Write a renamed copy."""
    dest.write_text(tgt_path.read_text())
    retag(dest, genome_id)
    return dest


def digest_one(task):
    isolate, asm_path, tgt_dir, syn2b, enzymes = task
    out_tgt = Path(tgt_dir) / f"{isolate}.tgt"
    if out_tgt.exists() and out_tgt.stat().st_size > 0:
        return isolate, None
    try:
        # Text, not binary: `syn2b synteny` reads only the text TGT path and does
        # not sniff the format, so a binary TGT fails with the unhelpful
        # "stream did not contain valid UTF-8".
        run([syn2b, "digest", "-i", str(asm_path), "-o", str(out_tgt),
             "-e", enzymes, "-f", "text"], timeout=1800)
        retag(out_tgt, isolate)
        return isolate, None
    except Exception as e:                       # noqa: BLE001 - reported, not raised
        return isolate, str(e)


def synteny_pair(task):
    """Compare two TGTs. `a` becomes genome_A, the fixed reference for raw_inverted_fraction."""
    cohort, a, b, tgt_dir, work_dir, syn2b = task
    tgt_dir, work_dir = Path(tgt_dir), Path(work_dir)
    a_tgt, b_tgt = tgt_dir / f"{a}.tgt", tgt_dir / f"{b}.tgt"
    if not a_tgt.exists() or not b_tgt.exists():
        return {"cohort": cohort, "genome_A": a, "genome_B": b, "status": "missing_tgt"}

    pairdir = work_dir / "pairs" / cohort / f"{a}__{b}"
    pairdir.mkdir(parents=True, exist_ok=True)
    for stale in pairdir.glob("*.tgt"):
        stale.unlink()
    # Name so that `a` sorts first and becomes genome_A, the fixed reference for
    # raw_inverted_fraction.
    (pairdir / f"a_{a}.tgt").symlink_to(a_tgt.resolve())
    if a == b:
        # syn2b refuses two genomes with the same id, so the self-comparison needs
        # a renamed copy rather than a second symlink.
        self_copy(b_tgt, pairdir / f"b_{b}__self.tgt", f"{b}__self")
    else:
        (pairdir / f"b_{b}.tgt").symlink_to(b_tgt.resolve())

    out_csv = pairdir / "synteny.csv"
    try:
        run([syn2b, "synteny", "-i", str(pairdir), "-o", str(out_csv)], timeout=600)
    except Exception as e:                       # noqa: BLE001
        return {"cohort": cohort, "genome_A": a, "genome_B": b,
                "status": f"synteny_error: {e}"}

    try:
        with open(out_csv) as fh:
            rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
        if not rows:
            return {"cohort": cohort, "genome_A": a, "genome_B": b, "status": "no_data"}
        r = rows[0]
    except Exception as e:                       # noqa: BLE001
        return {"cohort": cohort, "genome_A": a, "genome_B": b,
                "status": f"parse_error: {e}"}

    out = {"cohort": cohort, "genome_A": a, "genome_B": b, "status": "ok",
           "is_self": int(a == b)}
    for key in ("breakpoints", "scj_distance", "breakpoint_density",
                "inverted_fraction", "raw_inverted_fraction",
                "orientation_mismatches", "orientation_mismatches_raw",
                "orientation_uninformative", "observable_fraction",
                "observable_adjacencies", "structural", "shared_tags",
                "repeats_dropped", "landmarks_collapsed", "circular"):
        out[f"syn2b_{key}"] = r.get(key, "NA")
    return out


def read_samples(samples_dir: Path, cohort: str):
    """Return (isolate, group) pairs; group is the host/participant if the file has one."""
    path = samples_dir / f"samples_{cohort}.tsv"
    if not path.exists():
        return []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        group_col = next((c for c in ("host", "participant", "subject")
                          if c in (reader.fieldnames or [])), None)
        return [(r["isolate"], r[group_col] if group_col else "") for r in reader]


def fnum(row, key):
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return float("nan")


def spearman(xs, ys):
    """Rank correlation with average ranks for ties. No scipy dependency."""
    pts = [(x, y) for x, y in zip(xs, ys) if x == x and y == y]
    if len(pts) < 3:
        return float("nan")

    def rank(vals):
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[order[k]] = avg
            i = j + 1
        return ranks

    rx, ry = rank([p[0] for p in pts]), rank([p[1] for p in pts])
    n = len(pts)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def median(vals):
    vals = sorted(v for v in vals if v == v)
    if not vals:
        return float("nan")
    n = len(vals)
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2.0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--assembly-dir", required=True,
                   help="directory of {isolate}.fna assemblies")
    p.add_argument("--samples-dir", required=True,
                   help="directory of samples_{cohort}.tsv metadata")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--syn2b", default="/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b")
    p.add_argument("--enzymes", default=DEFAULT_ENZYMES,
                   help="enzyme panel for `syn2b digest`. The TGT cache keys on "
                        "isolate alone, so use a separate --out-dir per panel.")
    p.add_argument("--label", default="raw",
                   help="tag for this assembly set, e.g. raw or ref_ordered. Goes "
                        "into the output filenames and a column.")
    p.add_argument("--cohorts", nargs="*", default=COHORTS)
    p.add_argument("--workers", type=int, default=min(16, cpu_count() or 1))
    p.add_argument("--allow-failed-controls", action="store_true",
                   help="report cohort statistics even when self-comparisons are "
                        "not perfectly collinear. Off by default: a nonzero self "
                        "floor is exactly what invalidated the earlier pass.")
    args = p.parse_args()

    asm_dir = Path(args.assembly_dir)
    samples_dir = Path(args.samples_dir)
    out_dir = Path(args.out_dir)
    tgt_dir = out_dir / "tgt"
    tgt_dir.mkdir(parents=True, exist_ok=True)

    # ---- enumerate ----------------------------------------------------------
    cohort_samples, groups = {}, {}
    for cohort in args.cohorts:
        found = []
        for isolate, group in read_samples(samples_dir, cohort):
            asm = asm_dir / f"{isolate}.fna"
            if asm.exists() and asm.stat().st_size > 0:
                found.append(isolate)
                groups[(cohort, isolate)] = group
            else:
                print(f"WARN {cohort}: missing assembly {asm}", file=sys.stderr)
        cohort_samples[cohort] = found
        print(f"{cohort}: {len(found)} assemblies", flush=True)

    # ---- digest -------------------------------------------------------------
    tasks = [(iso, asm_dir / f"{iso}.fna", tgt_dir, args.syn2b, args.enzymes)
             for cohort in args.cohorts for iso in cohort_samples[cohort]]
    seen, dedup = set(), []
    for t in tasks:
        if t[0] not in seen:
            seen.add(t[0])
            dedup.append(t)
    print(f"digesting {len(dedup)} assemblies ({args.enzymes}) "
          f"with {args.workers} workers ...", flush=True)
    with Pool(args.workers) as pool:
        for i, (iso, err) in enumerate(pool.imap_unordered(digest_one, dedup, chunksize=4), 1):
            if err:
                print(f"  digest FAILED {iso}: {err}", file=sys.stderr)
            if i % 20 == 0:
                print(f"  digested {i}/{len(dedup)}", flush=True)

    # ---- self-comparison control, first and blocking -------------------------
    print("\n=== self-comparison control ===", flush=True)
    self_tasks = [(cohort, iso, iso, tgt_dir, out_dir, args.syn2b)
                  for cohort in args.cohorts for iso in cohort_samples[cohort]]
    with Pool(args.workers) as pool:
        self_rows = list(pool.imap_unordered(synteny_pair, self_tasks, chunksize=2))

    control = {}
    for cohort in args.cohorts:
        rows = [r for r in self_rows if r["cohort"] == cohort and r.get("status") == "ok"]
        if not rows:
            control[cohort] = {"n": 0, "passed": False, "worst": {}}
            continue
        worst = {}
        for key, want in SELF_EXPECT.items():
            worst[key] = max(abs(fnum(r, key) - want) for r in rows)
        passed = all(v <= SELF_TOL for v in worst.values())
        control[cohort] = {
            "n": len(rows), "passed": passed, "worst": worst,
            "median_bp": median([fnum(r, "syn2b_breakpoints") for r in rows]),
            "median_scj": median([fnum(r, "syn2b_scj_distance") for r in rows]),
            "median_obsfrac": median([fnum(r, "syn2b_observable_fraction") for r in rows]),
        }
        obs = control[cohort]["median_obsfrac"]
        shared = median([fnum(r, "syn2b_shared_tags") for r in rows])
        # observable_fraction ~ 1 - (K-1)/S  =>  K ~ 1 + (1-obs)*S
        control[cohort]["K_est"] = (1.0 + (1.0 - obs) * shared) if obs == obs else float("nan")

    print(f"{'cohort':<32} {'n':>4} {'self bp':>9} {'self SCJ':>9} "
          f"{'obs_frac':>9} {'K_est':>8}  control")
    for cohort in args.cohorts:
        c = control[cohort]
        if not c["n"]:
            print(f"{cohort:<32} {0:>4} {'-':>9} {'-':>9} {'-':>9} {'-':>8}  NO DATA")
            continue
        print(f"{cohort:<32} {c['n']:>4} {c['median_bp']:>9.1f} {c['median_scj']:>9.1f} "
              f"{c['median_obsfrac']:>9.4f} {c['K_est']:>8.0f}  "
              f"{'PASS' if c['passed'] else 'FAIL'}")

    failed = [c for c in args.cohorts if not control[c]["passed"]]
    if failed:
        print(f"\nSelf-comparison control FAILED for: {', '.join(failed)}", file=sys.stderr)
        print("A genome compared against itself must give 0 breakpoints, SCJ 0 and "
              "inverted_fraction 0. (observable_fraction is allowed to be <1 on draft "
              "assemblies; it estimates contig count.) A nonzero floor on the three "
              "fragmentation-immune metrics means the cohort's structural numbers "
              "measure assembly fragmentation, not biology -- which is precisely how "
              "the earlier Syn2bANI pass produced a species ranking that matched its "
              "own self-comparison floor to within a few percent.", file=sys.stderr)

    out_dir.mkdir(parents=True, exist_ok=True)
    ctrl_path = out_dir / f"self_control_{args.label}.tsv"
    with open(ctrl_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "label", "n", "passed", "median_breakpoints",
                    "median_scj", "median_observable_fraction", "K_est",
                    *(f"max_dev_{k}" for k in SELF_EXPECT)])
        for cohort in args.cohorts:
            c = control[cohort]
            w.writerow([cohort, args.label, c["n"], int(c["passed"]),
                        c.get("median_bp", ""), c.get("median_scj", ""),
                        c.get("median_obsfrac", ""), c.get("K_est", ""),
                        *(c["worst"].get(k, "") for k in SELF_EXPECT)])
    print(f"\nwrote {ctrl_path}")

    if failed and not args.allow_failed_controls:
        print("\nRefusing to report cohort statistics. Fix the floor, or re-run with "
              "--allow-failed-controls to inspect the numbers anyway.", file=sys.stderr)
        return 2

    # ---- all within-cohort pairs --------------------------------------------
    pair_tasks = []
    for cohort in args.cohorts:
        for a, b in itertools.combinations(sorted(cohort_samples[cohort]), 2):
            pair_tasks.append((cohort, a, b, tgt_dir, out_dir, args.syn2b))
    print(f"\ncomparing {len(pair_tasks)} within-cohort pairs ...", flush=True)
    with Pool(args.workers) as pool:
        rows = []
        for i, r in enumerate(pool.imap_unordered(synteny_pair, pair_tasks, chunksize=4), 1):
            rows.append(r)
            if i % 200 == 0:
                print(f"  {i}/{len(pair_tasks)}", flush=True)

    for r in rows:
        r["label"] = args.label
        r["same_group"] = int(groups.get((r["cohort"], r["genome_A"]), "") ==
                              groups.get((r["cohort"], r["genome_B"]), "") != "")
        r["group_A"] = groups.get((r["cohort"], r["genome_A"]), "")
        r["group_B"] = groups.get((r["cohort"], r["genome_B"]), "")

    # Take the union across rows, not rows[0]: an error row carries no syn2b_*
    # keys, and DictWriter(extrasaction="ignore") would drop the whole structural
    # channel without saying so.
    syn2b_cols = []
    for r in rows:
        for k in r:
            if k.startswith("syn2b_") and k not in syn2b_cols:
                syn2b_cols.append(k)
    fields = ["cohort", "label", "genome_A", "genome_B", "group_A", "group_B",
              "same_group", "is_self", "status"] + syn2b_cols
    pairs_path = out_dir / f"syn2b_structural_pairs_{args.label}.tsv"
    with open(pairs_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {pairs_path}  ({len(rows)} pairs)")

    # ---- per-cohort summary, floor-subtracted -------------------------------
    summary_path = out_dir / f"syn2b_structural_summary_{args.label}.tsv"
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "label", "role", "n_pairs", "self_floor_bp",
                    "median_bp", "excess_over_floor", "median_scj",
                    "median_inverted_fraction", "median_min_f_1mf",
                    "median_observable_fraction", "K_est",
                    "rho_scj_vs_invfrac", "expectation"])
        print(f"\n{'cohort':<32} {'role':<18} {'n':>5} {'floor':>7} {'med bp':>8} "
              f"{'excess':>8} {'med f':>7} {'K_est':>7}")
        for cohort in args.cohorts:
            ok = [r for r in rows if r["cohort"] == cohort and r.get("status") == "ok"]
            if not ok:
                continue
            bp = [fnum(r, "syn2b_breakpoints") for r in ok]
            scj = [fnum(r, "syn2b_scj_distance") for r in ok]
            invf = [fnum(r, "syn2b_inverted_fraction") for r in ok]
            minf = [min(v, 1.0 - v) for v in
                    (fnum(r, "syn2b_raw_inverted_fraction") for r in ok) if v == v]
            obs = [fnum(r, "syn2b_observable_fraction") for r in ok]
            floor = control[cohort].get("median_bp", float("nan"))
            med_bp = median(bp)
            role, note = EXPECTATION.get(cohort, ("", ""))
            w.writerow([cohort, args.label, role, len(ok), f"{floor:.1f}",
                        f"{med_bp:.1f}", f"{med_bp - floor:.1f}",
                        f"{median(scj):.1f}", f"{median(invf):.4f}",
                        f"{median(minf):.4f}", f"{median(obs):.4f}",
                        f"{control[cohort].get('K_est', float('nan')):.0f}",
                        f"{spearman(scj, invf):.4f}", note])
            print(f"{cohort:<32} {role:<18} {len(ok):>5} {floor:>7.1f} {med_bp:>8.1f} "
                  f"{med_bp - floor:>8.1f} {median(invf):>7.4f} "
                  f"{control[cohort].get('K_est', float('nan')):>7.0f}")
    print(f"\nwrote {summary_path}")
    print("\n`excess_over_floor` is the only column that can carry biology. If it is "
          "~0 for S. rimosus the positive control has failed; if it is large for the "
          "E. coli hypermutator set the negative control has failed.")
    print("`median_scj` carries a +(K-1) fragmentation term that `median_bp` does "
          "not -- on draft assemblies read the breakpoint columns, not SCJ.")
    print("`median_min_f_1mf` is only interpretable where K_est is small; on "
          "fragmented input it drifts to 0.5 and carries no biology.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
