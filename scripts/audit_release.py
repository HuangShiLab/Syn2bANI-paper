#!/usr/bin/env python3
"""Binary release audit for the Syn2bANI paper package.

The audit is deliberately independent of pandas and reads the released tables
with the standard library so that a clean checkout can run it without an
environment installation. Exit status zero means all mandatory gates pass.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "paper/manuscript/manuscript_nature_methods.md"
SUPPLEMENT = ROOT / "paper/manuscript/SUPPLEMENTARY.md"
REPRO = ROOT / "REPRODUCIBILITY_MANIFEST.md"
FAILURES: list[str] = []
CHECKS = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        FAILURES.append(f"{name}: {detail}" if detail else name)


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def read_rows(p: Path, delimiter="\t"):
    with p.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def finite_mae(rows, observed_col: str, truth_col: str, factor: float = 1.0):
    errors = []
    for r in rows:
        try:
            observed = float(r[observed_col])
            truth = float(r[truth_col]) * factor
        except (ValueError, TypeError, KeyError):
            continue
        errors.append(observed - truth)
    return len(errors), sum(abs(x) for x in errors) / len(errors), sum(errors) / len(errors)


def main() -> int:
    check("root-files", all((ROOT / x).exists() for x in [
        "PROJECT_FRAME.md", "PAPER_ARCHITECTURE.md", "CLAIM_EVIDENCE_MATRIX.md",
        "REPRODUCIBILITY_MANIFEST.md", "paper/manuscript/manuscript_nature_methods.md",
    ]), "core release packet missing")

    m = read(MANUSCRIPT)
    supp = read(SUPPLEMENT)
    check("tool-name", "Syn2b-ANI" not in m and "Syn2b-ANI" not in supp, "old tool name")
    check("em-dashes", "—" not in m and "—" not in supp, "prose em dash")

    # Exact truth and comparator version values.
    for phrase in [
        "0.073 ANI points versus 0.477 for skani 0.3.2",
        "20–200 contigs",
        "Benchmark-only v6 calibration is not embedded",
        "MAE: Syn2bANI 0.073%, skani 0.3.2 0.477%, FastANI 0.742%",
        "0.619 points against ANIm versus 0.957 for skani and 0.977 for FastANI",
        "183 non-reference genomes (127 curated complete and 56 curated deleted/peripheral-truncated)",
    ]:
        check(f"manuscript-phrase:{phrase[:38]}", phrase in m, "missing or stale claim")

    for stale in ["0.377 for skani", "20–201 contigs", "abfA loss", "Manuscript draft —", "/Users/", "/Volumes/"]:
        check(f"stale:{stale}", stale not in m, "stale host/process/value string")

    # Display and citation architecture.
    body = m[:m.index("## Online Methods")]
    ref_text = m[m.index("## References"):m.index("## Figure Legends")]
    refs = [int(x) for x in re.findall(r"^(\d+)\.\s", ref_text, re.M)]
    check("references-46", refs == list(range(1, 47)), f"references {len(refs)} or nonconsecutive")
    check("main-figures", sorted(set(map(int, re.findall(r"\bFig\. (\d+)", body)))) == list(range(1, 7)), "main figure refs")
    check("ed-figures", sorted(set(map(int, re.findall(r"Extended Data Fig\. (\d+)", body)))) == list(range(1, 6)), "ED figure refs")
    check("main-tables", sorted(set(map(int, re.findall(r"\bTable (\d+)", body)))) == [1, 2], "main table refs")

    # Repository paths cited from either manuscript document must resolve.
    for label, text in [("main-paths", m), ("supplement-paths", supp)]:
        missing = []
        for token in re.findall(r"`([^`\n]+)`", text):
            if re.match(r"^(?:results|paper|case_studies|scripts|analysis)/", token) and " " not in token and "--" not in token:
                if not (ROOT / token).exists():
                    missing.append(token)
        check(label, not missing, f"missing {sorted(set(missing))[:10]}")

    # Verify all checksums declared in the reproducibility manifest.
    checksum_text = "\n".join(re.findall(r"```text\n(.*?)```", read(REPRO), re.S))
    checksum_entries = []
    for line in checksum_text.splitlines():
        if "  " in line and "/" in line:
            sha, rel = line.split("  ", 1)
            checksum_entries.append((sha, rel.strip()))
    check("manifest-has-checksums", len(checksum_entries) >= 20, f"only {len(checksum_entries)} entries")
    bad_sha = []
    for sha, rel in checksum_entries:
        p = ROOT / rel
        if not p.exists():
            bad_sha.append(f"{rel}:missing")
            continue
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        if actual != sha:
            bad_sha.append(f"{rel}:checksum")
    check("manifest-checksums", not bad_sha, f"bad {bad_sha[:10]}")

    # Recalculate Figure 2 source values.
    sim = ROOT / "results/simulations"
    man = read_rows(sim / "simindel/manifest.tsv")
    syn = read_rows(sim / "simindel_results_4e.tsv")
    cross = read_rows(sim / "simindel_cross_tool_4e.tsv")
    truth = {r["name"]: float(r["true_ani"]) * 100 for r in man}
    syn_by_name = {r["query"]: r for r in syn}
    common = sorted(set(truth) & set(syn_by_name) & {r["name"] for r in cross})
    check("simulation-common-n", len(common) == 12, f"n={len(common)}")
    for method, rows, col in [
        ("Syn2bANI", syn_by_name, "ani"),
        ("skani", {r["name"]: r for r in cross}, "skani_ani"),
        ("FastANI", {r["name"]: r for r in cross}, "fastani_ani"),
    ]:
        errs = [abs(float(rows[n][col]) - truth[n]) for n in common]
        expected = {"Syn2bANI": 0.0732, "skani": 0.4775, "FastANI": 0.7423}[method]
        check(f"simulation-mae:{method}", abs(sum(errs) / len(errs) - expected) < 0.00055,
              f"got {sum(errs) / len(errs):.6f}, expected {expected}")

    # S10 current rerun values.
    saureus = read_rows(ROOT / "case_studies/fda_argos_s_aureus/results/triangle_rerun.tsv")
    bp = [float(r["breakpoint_count"]) for r in saureus]
    check("saureus-n", len(saureus) == 2211, f"n={len(saureus)}")
    check("saureus-zero", sum(x == 0 for x in bp) / len(bp) > 0.604 and sum(x == 0 for x in bp) / len(bp) < 0.606, "zero fraction")
    check("saureus-max", max(bp) == 7, f"max={max(bp)}")

    # OAPGC overall effect values.
    oapgc = {r["metric"]: float(r["value"]) for r in read_rows(ROOT / "results/census/oapgc/site_effects/overall_all_pairs.tsv")}
    for k, expected in {
        "species_tested": 200,
        "median_SNP_auc": 0.5737193778971896,
        "median_BP_auc": 0.49415237243802246,
        "n_OA_pairs": 1264327,
        "n_within_pairs": 16067030,
    }.items():
        check(f"oapgc:{k}", abs(oapgc.get(k, float("nan")) - expected) < 1e-9, f"got {oapgc.get(k)}")

    # GTDB joint table row count includes header.
    with (ROOT / "results/census/gtdb/gtdb_species_ani_sv_joint.tsv").open() as fh:
        joint_rows = sum(1 for _ in fh) - 1
    check("gtdb-joint-clusters", joint_rows == 21973, f"rows={joint_rows}")

    # Inspect textual content of the packaged main Figure 2 and S10 title through pdftotext.
    if shutil := __import__("shutil").which("pdftotext"):
        fig2 = subprocess.check_output(["pdftotext", str(ROOT / "paper/figures/main/final/figure2.pdf"), "-"], text=True)
        for expected in ["MAE Syn2bANI 0.073", "MAE skani 0.477", "MAE FastANI 0.742"]:
            check(f"fig2-render:{expected}", expected in fig2, "rendered label missing")
        s10 = subprocess.check_output(["pdftotext", str(ROOT / "paper/figures/supplementary/fig_s10_saureus_breakpoints.pdf"), "-"], text=True)
        check("fig-s10-render", "60.5% of pairs have zero breakpoints" in s10 and "0.1.0" not in s10, "stale S10 rendering")
    else:
        check("pdftotext-available", False, "pdftotext is required for release rendering gate")

    # Code/data release plan is explicit and not represented only by a commit hash.
    check("software-archive-plan", "Zenodo DOI for the code and embedded model will be added before acceptance" in m, "missing code archive release condition")
    check("census-archive-plan", "complete pair-level census tables and their checksums/accession manifests will be deposited as an immutable Zenodo archive before acceptance" in m, "missing census archive release condition")

    result = {
        "verdict": "RELEASE READY" if not FAILURES else "RELEASE BLOCKED",
        "checks": CHECKS,
        "failures": FAILURES,
    }
    out = ROOT / "reviews/RELEASE_CHECK_LATEST.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    raise SystemExit(main())
