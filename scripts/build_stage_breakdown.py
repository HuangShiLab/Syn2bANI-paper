#!/usr/bin/env python3
"""Build stage-level runtime decomposition + Stage-1 screen recall tables.

Outputs:
  results/efficiency_v8/stage_breakdown.tsv  (Task A)
  results/efficiency_v8/screen_recall.tsv    (Task B)

Data sources (all local, read-only):
  results/efficiency_v8/runtime_scaling.tsv      syn2bani ani wall times (FASTA & sketch modes)
  results/efficiency_v8/sketch_benchmark.tsv     syn2bani sketch command wall times
  results/efficiency_v8/syn2b_struct_benchmark.tsv  `syn2bani struct` (SV stage) wall times
  results/db_scale/DB_REWRITE_VALIDATION.md      screen pass rates / FRR numbers (regex-parsed)
  results/db_scale/rewrite/rewrite_scaling*.tsv  rewritten-pipeline db-scale wall times
  results/gtdb50k/                               checked: no 50k-scale screen numbers exist

Important honesty note: syn2bani emits no per-stage timing (its stdout/stderr logs are
empty; `--verbose` only adds per-pair diagnostic columns). Stage 1 (digestion) and
stage 2 (sketch build/write) are only measured *combined* via `syn2bani sketch`;
stages 3-5 (screen / chaining / MLE+calibration) are only measured *combined* via
`syn2bani ani`. Derived rows and bounds are labelled as such in `timing_status`.
"""

import csv
import re
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent.parent
EFF = ROOT / "results" / "efficiency_v8"
DB_SCALE_MD = ROOT / "results" / "db_scale" / "DB_REWRITE_VALIDATION.md"

# ---------------------------------------------------------------------------
# Task A: stage breakdown
# ---------------------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def med(rows, key="wall_s"):
    vals = sorted(float(r[key]) for r in rows)
    return median(vals), vals[0], vals[-1]


def stage_breakdown():
    runtime = read_tsv(EFF / "runtime_scaling.tsv")
    sketch = read_tsv(EFF / "sketch_benchmark.tsv")
    struct = read_tsv(EFF / "syn2b_struct_benchmark.tsv")

    out = []  # dict rows

    def add(stage, scale, n_genomes, n_pairs, rows, status, note, source, key="wall_s"):
        if rows:
            m, lo, hi = med(rows, key)
            out.append({
                "stage": stage, "scale": scale,
                "n_genomes": n_genomes, "n_pairs": n_pairs,
                "wall_s_median": f"{m:.3f}", "wall_s_min": f"{lo:.3f}",
                "wall_s_max": f"{hi:.3f}",
                "timing_status": status, "derivation_note": note, "source": source,
            })
        else:
            out.append({
                "stage": stage, "scale": scale,
                "n_genomes": n_genomes, "n_pairs": n_pairs,
                "wall_s_median": "NA", "wall_s_min": "NA", "wall_s_max": "NA",
                "timing_status": status, "derivation_note": note, "source": source,
            })

    syn = [r for r in runtime if r["tool"] == "syn2bani"]
    by_mode_n = {}
    for r in syn:
        by_mode_n.setdefault((r["mode"], r["n_genomes"]), []).append(r)

    # Measured whole-pipeline anchors for every panel size (n=22 panel = 484 pairs).
    for (mode, n), rows in sorted(by_mode_n.items(), key=lambda kv: (kv[0][0], int(kv[0][1]))):
        if mode == "ani_fasta":
            stage, note = (
                "total_end_to_end_fasta_mode",
                "stages 1+3+4+5 combined: in-memory digestion of both lists + screen + chaining + MLE (no sketch I/O)",
            )
        else:
            stage, note = (
                "total_end_to_end_sketch_mode",
                "stages 3+4+5 combined + .s2ba load: screen + chaining + MLE (digestion done offline by `sketch`)",
            )
        add(stage, f"n{n}", n, rows[0]["n_pairs"], rows, "measured",
            note + "; 3 reps, /usr/bin/time -l",
            f"results/efficiency_v8/runtime_scaling.tsv (tool=syn2bani, mode={mode}, n={n})")

    # Sketch command: stages 1+2 combined (digest + sketch build + disk write).
    s2 = [r for r in sketch if r["tool"] == "syn2bani"]
    for n in sorted({r["n_genomes"] for r in s2}, key=int):
        rows = [r for r in s2 if r["n_genomes"] == n]
        add("sketch_cmd_digest_build_write", f"n{n}", n, "NA", rows, "measured",
            "stages 1+2 combined (in-silico digestion + sketch build + disk write); "
            "NOT separable into digestion vs build/write from these data",
            f"results/efficiency_v8/sketch_benchmark.tsv (tool=syn2bani, n={n})")

    # Derived: in-memory digestion ~= ani_fasta - ani_sketches (medians per n).
    for n in sorted({k[1] for k in by_mode_n}, key=int):
        fa = by_mode_n.get(("ani_fasta", n))
        sk = by_mode_n.get(("ani_sketches", n))
        if fa and sk:
            est = med(fa)[0] - med(sk)[0]
            out.append({
                "stage": "in_silico_digestion", "scale": f"n{n}",
                "n_genomes": n, "n_pairs": "NA",
                "wall_s_median": f"{est:.3f}", "wall_s_min": "NA", "wall_s_max": "NA",
                "timing_status": "derived_upper_bound",
                "derivation_note": (
                    "median(ani_fasta) - median(ani_sketches); upper bound because the "
                    "sketch-mode term also pays .s2ba load I/O and FASTA parsing is bundled "
                    "in here. Cross-check is inconsistent at sub-second scale: the whole "
                    "`sketch` command (digest+build+write of the same genomes) measured "
                    "0.50-0.51 s at n=22, so treat ~0.9 s as I/O-noisy."
                ),
                "source": f"results/efficiency_v8/runtime_scaling.tsv (derived, n={n})",
            })

    # n=22 panel: stages with no independent timing -> honest bounds.
    n22_sk = med(by_mode_n[("ani_sketches", "22")])[0]
    n2_sk = med(by_mode_n[("ani_sketches", "2")])[0]
    marg = (n22_sk - n2_sk) / (484 - 4)
    out.append({
        "stage": "screen_pass1_index_lookup", "scale": "n22", "n_genomes": "22",
        "n_pairs": "484", "wall_s_median": "NA", "wall_s_min": "NA", "wall_s_max": "NA",
        "timing_status": "not_separately_timed",
        "derivation_note": (
            f"merged with chaining+MLE: screen+chain+MLE together = {n22_sk:.2f} s (sketch mode). "
            f"Bound: 0 < t_screen < {n22_sk:.2f} s. Cross-scale anchor: rewritten triangle at "
            "n=500 screens 124,750 pairs AND refines 21,595 in 8.5 s total, so screen is not "
            "the dominant cost even at 250x the pair count; no exact split available locally."
        ),
        "source": "results/efficiency_v8/runtime_scaling.tsv (no per-stage log output exists)",
    })
    out.append({
        "stage": "chaining_anchor_dp", "scale": "n22", "n_genomes": "22",
        "n_pairs": "484", "wall_s_median": "NA", "wall_s_min": "NA", "wall_s_max": "NA",
        "timing_status": "not_separately_timed",
        "derivation_note": (
            "single chain_ani::compute call per pair, merged with screen and MLE in the "
            f"{n22_sk:.2f} s sketch-mode total; no independent timer in any local log"
        ),
        "source": "results/efficiency_v8/runtime_scaling.tsv",
    })
    out.append({
        "stage": "mle_fit_plus_calibration", "scale": "n22", "n_genomes": "22",
        "n_pairs": "484", "wall_s_median": "NA", "wall_s_min": "NA", "wall_s_max": "NA",
        "timing_status": "not_separately_timed",
        "derivation_note": (
            "MLE fit shares the per-pair chain_ani::compute call; --calibrate applies a linear "
            "model (sub-ms/pair, not measurable at this scale). Merged total for stages 3+4+5 "
            f"= {n22_sk:.2f} s. Marginal per-pair cost of stages 3+4+5 derived as "
            f"{marg*1000:.2f} ms/pair = (median(ani_sk,n22) - median(ani_sk,n2)) / (484-4)."
        ),
        "source": "results/efficiency_v8/runtime_scaling.tsv (derived)",
    })

    # Extra SV stage (`syn2bani struct`), n=22 panel.
    st22 = [r for r in struct if r["n_genomes"] == "22"]
    add("structural_sv_stage_syn2b_struct", "n22", "22", "484", st22, "measured",
        "extra structural/SV stage beyond the 5 ANI stages; run via python driver spawning "
        "one `syn2bani struct` process per pair (16-way pool) -> includes process-spawn overhead",
        "results/efficiency_v8/syn2b_struct_benchmark.tsv (n=22 rows)", key="struct_wall_s")

    # DB-scale anchors: rewritten screen+refine pipeline, HPC amd 32 threads.
    md = DB_SCALE_MD.read_text()
    anchor_note = "screen+refine combined (stages 3+4+5), no finer split available; syn2bani built from rewrite working tree on top of 15da386"
    m = re.search(r"triangle n=500 \| \*\*([\d.]+) s", md)
    if m:
        out.append({"stage": "screen_plus_refine_db_scale", "scale": "n500", "n_genomes": "500",
                    "n_pairs": "124750", "wall_s_median": m.group(1), "wall_s_min": "NA",
                    "wall_s_max": "NA", "timing_status": "measured_total_only",
                    "derivation_note": anchor_note + "; 21,595 pairs pass screen (17.3%)",
                    "source": "results/db_scale/DB_REWRITE_VALIDATION.md (HPC scale check table)"})
    m = re.search(r"triangle n=2000 \| \*\*([\d.]+) s", md)
    if m:
        out.append({"stage": "screen_plus_refine_db_scale", "scale": "n2000", "n_genomes": "2000",
                    "n_pairs": "1999000", "wall_s_median": m.group(1), "wall_s_min": "NA",
                    "wall_s_max": "NA", "timing_status": "measured_total_only",
                    "derivation_note": anchor_note + "; 305,996 pairs pass screen (15.3%)",
                    "source": "results/db_scale/DB_REWRITE_VALIDATION.md (HPC scale check table)"})
    m = re.search(r"triangle n=5000 \| \*\*([\d.]+) s", md)
    if m:
        out.append({"stage": "screen_plus_refine_db_scale", "scale": "n5000", "n_genomes": "5000",
                    "n_pairs": "12497500", "wall_s_median": m.group(1), "wall_s_min": "NA",
                    "wall_s_max": "NA", "timing_status": "measured_total_only",
                    "derivation_note": anchor_note + "; 1,570,949 pairs pass screen (12.6%)",
                    "source": "results/db_scale/DB_REWRITE_VALIDATION.md (HPC scale check table)"})
    m = re.search(r"search 100×5000 \| \*\*([\d.]+) s", md)
    if m:
        out.append({"stage": "screen_plus_refine_search", "scale": "100x5000", "n_genomes": "5100",
                    "n_pairs": "500000", "wall_s_median": m.group(1), "wall_s_min": "NA",
                    "wall_s_max": "NA", "timing_status": "measured_total_only",
                    "derivation_note": anchor_note + "; 184 hits at --min-ani 0.8",
                    "source": "results/db_scale/DB_REWRITE_VALIDATION.md (HPC scale check table)"})
    m = re.search(r"5000-genome DB was re-sketched[^;]*`s2ba_n5000_new`, 479 MB store, (\d+) s", md)
    if m:
        out.append({"stage": "sketch_cmd_digest_build_write", "scale": "n5000", "n_genomes": "5000",
                    "n_pairs": "NA", "wall_s_median": m.group(1), "wall_s_min": "NA",
                    "wall_s_max": "NA", "timing_status": "measured",
                    "derivation_note": "stages 1+2 combined; new default 4-enzyme panel; 479 MB store, 1.2 GB RSS",
                    "source": "results/db_scale/DB_REWRITE_VALIDATION.md (screen calibration section)"})

    path = EFF / "stage_breakdown.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {path} ({len(out)} rows)")


# ---------------------------------------------------------------------------
# Task B: Stage-1 screen recall / pass rate
# ---------------------------------------------------------------------------

def screen_recall():
    rows = []

    def add(metric, scale, num, den, value, source, notes):
        rows.append({"metric": metric, "scale": scale, "numerator": num,
                     "denominator": den, "value": value, "source": source,
                     "notes": notes})

    md_path = DB_SCALE_MD
    md = md_path.read_text()

    # Pass rates (regex over the validation report; numbers cross-checked against
    # results/db_scale/rewrite/triangle_n{500,2000,5000}_v3.tsv row counts where applicable).
    pr = re.search(
        r"Screen pass rates: ([\d,]+)/([\d,]+) \(([\d.]+)%\) at n=500;\s+"
        r"([\d,]+)/([\d,]+) \(([\d.]+)%\)\s+at n=2000;\s+"
        r"([\d,]+)/([\d,]+) \(([\d.]+)%\) at n=5000", md)
    if pr:
        g = [pr.group(i).replace(",", "") for i in range(1, 10)]
        for scale, num, den, pct in (("triangle n=500", g[0], g[1], g[2]),
                                     ("triangle n=2000", g[3], g[4], g[5]),
                                     ("triangle n=5000", g[6], g[7], g[8])):
            add("stage1_screen_pass_rate", scale, num, den, f"{float(pct)/100:.4f}",
                f"{md_path} (Screen pass rates paragraph)",
                "all-vs-all on random GTDB-R207 subsets; survivors are mostly same-family "
                "background that refine then measures or marks BELOW_DETECTION; raw per-pair "
                "screen dumps exist on HPC only (results/db_scale/rewrite/screen_{acc,n500}_w*.tsv)")
    else:
        add("stage1_screen_pass_rate", "triangle n=500/2000/5000", "NA", "NA", "NA",
            str(md_path), "PATTERN NOT FOUND in report -- verify manually")

    # FRR on 500 validated true pairs at the shipped gate.
    frr = re.search(r"\| shared≥3 AND cont≥0\.001 \| \*\*0/500\*\* \| ([\d.]+)% \|", md)
    if frr:
        add("stage1_screen_false_reject_rate_validated_true_pairs",
            "500 true pairs ANI 80-100 (accuracy_pairs.tsv)", "0", "500", "0.0000",
            f"{md_path} (screen calibration table)",
            "shipped gate W=18, shared>=3 AND containment>=0.001; weakest 80-85 band true pair "
            "has shared=29, cont=0.0061 (10x/6x above floors)")

    # Pass-all reconciliation at n=2000.
    rec = re.search(r"retains \*\*([\d,]+)/([\d,]+)\*\* pairs.*?false-reject rate on\s+estimator-reportable pairs: \*\*1/844 = ([\d.]+)%", md, re.S)
    if rec:
        add("stage1_screen_false_reject_rate_pass_all_reconciliation",
            "triangle n=2000 (forced pass-all screen run)", "1", "844",
            f"{float(rec.group(3))/100:.4f}",
            f"{md_path} (triangle n=2000 output reconciliation)",
            "lost pair JADJDU010000001.1 x DRLG01000001.1 (gated 98.44) had "
            "af_query=af_reference=0.0000 -- tiny shared-island call, INCONSISTENT-flagged; "
            "max |dANI| = 0.0 on retained pairs")

    # search recall vs skani.
    if re.search(r"184 hits at `--min-ani 0.8`", md):
        add("search_recall_vs_skani_hits", "search 100 queries x 5000-genome GTDB-R207 DB",
            "12", "12", "1.0000",
            f"{md_path} (search 100 queries x 5000 DB)",
            "all 12 skani-reported pairs found by rewritten search; on the 7 intersection pairs "
            "with independent MLE truth the new search matches exactly (91.09-95.11)")
        add("search_hits_min_ani_0.8", "search 100 queries x 5000-genome GTDB-R207 DB",
            "184", "500000", "0.000368",
            f"{md_path} (search 100 queries x 5000 DB)",
            "184 reported hits at --min-ani 0.8; screen-pass count for search not reported "
            "(hits = screen-passed AND refined AND ani_gated>=0.8); 13/184 have max(AF)>=0.15 "
            "(skani-comparable), 171 are low-AF calls")

    # Legacy screen contrast (pre-rewrite).
    if re.search(r"legacy exact-tag `min_af=0.1` screen rejected \*\*94.4%\*\* of true ≥80%", md):
        add("legacy_exact_tag_screen_true_pair_rejection_PRE_REWRITE",
            "500 true pairs ANI 80-100", "472", "500", "0.9440",
            f"{md_path} + results/db_scale/DB_SCALE_BENCHMARK.md (dist accuracy section)",
            "CONTRAST ROW: the retired v7 exact-tag min_af=0.1 screen; the new W=18 screen "
            "rejects 0/500 of the same pairs")

    # Post-refine reported rows (context for what survives the whole pipeline).
    for scale, f, reported in (("triangle n=500", "rewrite/triangle_n500_v3.tsv", "90"),
                               ("triangle n=2000", "rewrite/triangle_n2000_v3.tsv", "843"),
                               ("triangle n=5000", "rewrite/triangle_n5000_v3.tsv", "3410")):
        p = ROOT / "results" / "db_scale" / f
        if p.exists():
            n_rows = sum(1 for _ in open(p)) - 1  # minus header
            add("estimator_backed_rows_edge_list", scale, str(n_rows),
                {"triangle n=500": "124750", "triangle n=2000": "1999000",
                 "triangle n=5000": "12497500"}[scale],
                "see notes",
                f"results/db_scale/{f} (row count)",
                "pairs with a finite gated estimate after screen+refine (edge-list rows); "
                f"report text quotes {reported} for this scale -- matches" if str(n_rows) == reported
                else f"ROW COUNT {n_rows} != report's {reported} -- verify")

    # The 50k question, answered honestly.
    add("gtdb50k_scale_screen_pass_rate", "GTDB-R207 50k (all-vs-all or one-to-all)",
        "NA", "NA", "NOT_FOUND",
        "results/gtdb50k/ (whole directory), results/efficiency_v8/GTDB_SCALE_BENCHMARK_PLAN.md",
        "no 50k-scale screen numbers exist locally: the GTDB50k held-out benchmark ran `ani` on "
        "45,967 PRE-SELECTED pairs (results/gtdb50k/GTDB50K_HELDOUT_REPORT.md), i.e. the screen "
        "was never exercised at that scale; the full 65,703-genome all-vs-all (~2.16B pairs) was "
        "never run; the 7-query x 65,703 one-to-all benchmark is planned but its results file "
        "results/gtdb_one_to_all_scaling.tsv is absent from this repo (plan only, see GTDB_SCALE_"
        "BENCHMARK_PLAN.md 'Next steps'). Largest scale WITH screen data: n=5000 all-vs-all above.")

    path = EFF / "screen_recall.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    stage_breakdown()
    screen_recall()
