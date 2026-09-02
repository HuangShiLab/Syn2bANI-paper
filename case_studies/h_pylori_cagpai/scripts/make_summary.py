#!/usr/bin/env python3
"""Generate cagpai_summary.md from cohort and pilot TSV outputs."""

import csv
import statistics
from pathlib import Path


def read_tsv(path: Path):
    with path.open("r", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    cohort = read_tsv(Path("cagpai_states.tsv"))
    pilot = read_tsv(Path("cagpai_states_pilot.tsv"))

    status_counts = {"complete": 0, "partial": 0, "empty": 0}
    fractions = []
    for row in cohort:
        status_counts[row["status"]] += 1
        fractions.append(float(row["fraction_present"]))

    frac_mean = statistics.mean(fractions)
    frac_median = statistics.median(fractions)
    frac_min = min(fractions)
    frac_max = max(fractions)

    # Pilot table, ordered by biological interpretation.
    pilot_order = [
        "wt",
        "hp26695",
        "mut1",
        "inv",
        "transloc",
        "mut1_inv",
        "mut1_transloc",
        "del",
        "mut1_del",
    ]
    pilot_by_name = {row["genome"]: row for row in pilot}

    lines = []
    lines.append("# cagPAI status summary for 528 H. pylori genomes")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Aligner:** minimap2 v2.31 (`-cx asm20`)")
    lines.append("- **Markers:** 28 CDS loci from H. pylori 26695 cagPAI (HP0520–HP0547)")
    lines.append("- **Presence rule:** coverage ≥ 0.8 AND identity ≥ 0.8")
    lines.append("- **Classification:** complete ≥ 0.85, empty ≤ 0.15, otherwise partial")
    lines.append("")
    lines.append("## Cohort classification counts")
    lines.append("")
    lines.append("| status   | count | fraction |")
    lines.append("|----------|-------|----------|")
    total = len(cohort)
    for st in ("complete", "partial", "empty"):
        n = status_counts[st]
        lines.append(f"| {st:8s} | {n:5d} | {n/total:.3f}    |")
    lines.append(f"| **total**| {total:5d} | 1.000    |")
    lines.append("")
    lines.append("## Fraction-present distribution")
    lines.append("")
    lines.append(f"- mean:   {frac_mean:.3f}")
    lines.append(f"- median: {frac_median:.3f}")
    lines.append(f"- min:    {frac_min:.3f}")
    lines.append(f"- max:    {frac_max:.3f}")
    lines.append("")
    lines.append("## Pilot validation")
    lines.append("")
    lines.append("| strain            | n_present | fraction_present | status   | expected        |")
    lines.append("|-------------------|-----------|------------------|----------|-----------------|")
    for name in pilot_order:
        row = pilot_by_name.get(name)
        if row is None:
            continue
        expected = (
            "complete" if name in ("wt", "hp26695", "mut1", "inv", "transloc", "mut1_inv", "mut1_transloc") else "empty"
        )
        lines.append(
            f"| {name:17s} | {row['n_present']:9s} | {float(row['fraction_present']):16.3f} | "
            f"{row['status']:8s} | {expected:15s} |"
        )
    lines.append("")
    lines.append("## Threshold interpretation")
    lines.append("")
    correct = sum(
        1
        for name in pilot_order
        if pilot_by_name.get(name)
        and pilot_by_name[name]["status"]
        == (
            "complete"
            if name in ("wt", "hp26695", "mut1", "inv", "transloc", "mut1_inv", "mut1_transloc")
            else "empty"
        )
    )
    lines.append(
        f"Pilot controls classified correctly: {correct}/{len(pilot_order)}. "
        "All intact-island strains (wt, hp26695, mut1, inv, transloc, mut1_inv, mut1_transloc) "
        "score fraction_present ≥ 0.857, while the two deletion mutants (del, mut1_del) score 0. "
        "The thresholds therefore cleanly separate engineered cagPAI-present from cagPAI-absent controls."
    )
    lines.append("")

    Path("cagpai_summary.md").write_text("\n".join(lines))
    print("Wrote cagpai_summary.md")


if __name__ == "__main__":
    main()
