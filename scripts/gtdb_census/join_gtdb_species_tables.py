#!/usr/bin/env python3
"""Join GTDB structural and ANI species summaries into one QC table."""

import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        key = next((x for x in rows.fieldnames if x in ("species_cluster", "cluster")), None)
        if key is None:
            raise ValueError(f"no species key in {path}")
        return {row[key]: row for row in rows}


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structural", required=True)
    ap.add_argument("--ani", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    structural = read_tsv(args.structural)
    ani = read_tsv(args.ani)
    fields = [
        "species_cluster", "taxonomy", "n_genomes", "expected_pairs",
        "structural_rows", "structural_coverage", "sv_ge2_pairs",
        "sv_ge2_fraction", "median_breakpoints", "ani_pairs", "ani_coverage",
        "mean_skani_ani", "min_skani_ani", "max_skani_ani",
        "ani_sv_layer", "joint_qc_note",
    ]
    tmp = Path(str(args.output) + ".tmp")
    with tmp.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for species in sorted(set(structural) | set(ani)):
            srow = structural.get(species, {})
            arow = ani.get(species, {})
            expected = int(fnum(srow.get("effective_pairs", arow.get("expected_pairs", ""))) or 0)
            if expected == 0 and not srow.get("rows") and not arow.get("ani_pairs"):
                # QC writers emit zero-genome placeholder rows for taxonomy
                # clusters outside the effective scope; exclude these ghosts.
                continue
            struct_rows = int(fnum(srow.get("rows", "")) or 0)
            ani_rows = int(fnum(arow.get("ani_pairs", "")) or 0)
            struct_cov = struct_rows / expected if expected else None
            ani_cov = ani_rows / expected if expected else None
            pct_ge2 = fnum(srow.get("pct_ge2", ""))
            sv_ge2_pairs = round(struct_rows * pct_ge2 / 100.0) if pct_ge2 is not None else ""
            if expected and struct_cov == 1 and ani_cov == 1:
                layer, note = "complete", ""
            elif expected and ani_cov == 1:
                layer = "structural_partial"
                note = "structural layer partial"
            elif expected and struct_cov == 1:
                layer = "ani_partial"
                note = "ANI layer partial"
            elif expected:
                layer, note = "partial", "both layers incomplete"
            else:
                layer, note = "unplanned", "not represented in effective scope"
            writer.writerow({
                "species_cluster": species,
                "taxonomy": arow.get("taxonomy", ""),
                "n_genomes": fnum(srow.get("effective_genomes", arow.get("n_genomes", ""))) or "",
                "expected_pairs": expected or "",
                "structural_rows": struct_rows or "",
                "structural_coverage": struct_cov if struct_cov is not None else "",
                "sv_ge2_pairs": sv_ge2_pairs,
                "sv_ge2_fraction": pct_ge2 / 100 if pct_ge2 is not None else "",
                "median_breakpoints": srow.get("median_breakpoints", ""),
                "ani_pairs": ani_rows or "",
                "ani_coverage": ani_cov if ani_cov is not None else "",
                "mean_skani_ani": arow.get("mean_skani_ani", ""),
                "min_skani_ani": arow.get("min_skani_ani", ""),
                "max_skani_ani": arow.get("max_skani_ani", ""),
                "ani_sv_layer": layer,
                "joint_qc_note": note,
            })
    tmp.replace(args.output)


if __name__ == "__main__":
    main()
