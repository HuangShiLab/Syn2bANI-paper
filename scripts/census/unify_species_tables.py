#!/usr/bin/env python3
"""Unify GTDB/HROM/OAPGC species-level structural tables into one schema."""
import csv
from collections import Counter
from pathlib import Path

ROOT = Path("results/census")
OUTDIR = ROOT / "unified"
OUTDIR.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fnum(x, default=""):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    # ---- structural census table ----
    rows = []
    for r in read_tsv(ROOT / "gtdb/effective/effective_species_summary.tsv"):
        reported = int(r["rows"])
        ge2 = round(reported * float(r["pct_ge2"]) / 100.0)
        rows.append({
            "database": "GTDB_R207",
            "species_id": r["cluster"],
            "n_genomes": r["effective_genomes"],
            "planned_pairs": r["effective_pairs"],
            "reported_pairs": reported,
            "sv_ge2_pairs": ge2,
            "sv_ge2_fraction": ge2 / reported if reported else "",
            "median_breakpoints": r["median_breakpoints"],
            "ani_layer": "pending",
            "quality_note": "effective scope; one E. coli row-count exception",
        })

    for r in read_tsv(ROOT / "hrom/species_summary.tsv"):
        reported = int(r["rows"])
        pct = fnum(r["pct_pairs_ge2_junctions"])
        ge2 = round(reported * pct / 100.0)
        rows.append({
            "database": "HROM",
            "species_id": r["cluster"],
            "n_genomes": r["genomes"],
            "planned_pairs": r["expected_pairs"],
            "reported_pairs": reported,
            "sv_ge2_pairs": ge2,
            "sv_ge2_fraction": ge2 / reported if reported else "",
            "median_breakpoints": r["median_junctions"],
            "ani_layer": "complete",
            "quality_note": "",
        })

    oapgc_counts = {r["cluster"]: r for r in read_tsv(ROOT / "oapgc/oapgc_species_structural_counts.tsv")}
    oapgc_meta = {r["cluster"]: r for r in read_tsv("data/oapgc/oapgc_genome_metadata.tsv")}
    clusters = set(oapgc_meta)
    clusters.update(oapgc_counts)
    for cl in sorted(clusters):
        c = oapgc_counts.get(cl, {"pairs": 0, "pairs_ge2": 0, "fraction_ge2": 0})
        m = oapgc_meta.get(cl, {})
        rows.append({
            "database": "OAPGC",
            "species_id": cl,
            "n_genomes": sum(1 for r in read_tsv("data/oapgc/oapgc_genome_metadata.tsv") if False),  # filled below
            "planned_pairs": c.get("pairs", 0),
            "reported_pairs": c.get("pairs", 0),
            "sv_ge2_pairs": c.get("pairs_ge2", 0),
            "sv_ge2_fraction": c.get("fraction_ge2", ""),
            "median_breakpoints": "",
            "ani_layer": "complete",
            "quality_note": m.get("site_raw", ""),
        })

    # avoid re-reading metadata inside loop
    cluster_counts = Counter()
    for r in read_tsv("data/oapgc/oapgc_genome_metadata.tsv"):
        cluster_counts[r["cluster"]] += 1
    for row in rows:
        if row["database"] == "OAPGC":
            row["n_genomes"] = cluster_counts[row["species_id"]]

    fields = ["database", "species_id", "n_genomes", "planned_pairs", "reported_pairs",
              "sv_ge2_pairs", "sv_ge2_fraction", "median_breakpoints", "ani_layer", "quality_note"]
    write_tsv(OUTDIR / "species_structural_census_unified.tsv", rows, fields)

    # ---- effect/mode long table ----
    effects = []
    for r in read_tsv(ROOT / "hrom/quality_stratified/hypermode_species_high_quality_comp90_cont5.tsv"):
        effects.append({
            "database": "HROM",
            "species_id": r["cluster"],
            "species_name": r.get("species", ""),
            "analysis": "within_species_top5_high_quality",
            "axis": "SV",
            "effect_size": r.get("struct_enrichment_ratio", ""),
            "effect_direction": "higher_in_SV_top" if fnum(r.get("struct_enrichment_ratio")) > 1 else "lower_in_SV_top",
            "p_value": r.get("struct_p", ""),
            "q_value": r.get("struct_q", ""),
        })
        effects.append({
            "database": "HROM",
            "species_id": r["cluster"],
            "species_name": r.get("species", ""),
            "analysis": "within_species_top5_high_quality",
            "axis": "ANI",
            "effect_size": r.get("ani_enrichment_ratio", ""),
            "effect_direction": "higher_in_ANI_top" if fnum(r.get("ani_enrichment_ratio")) > 1 else "lower_in_ANI_top",
            "p_value": r.get("ani_p", ""),
            "q_value": r.get("ani_q", ""),
        })

    for tag, path in [
        ("all_pairs", ROOT / "oapgc/site_effects/species_effects_all_pairs.tsv"),
        ("n50_ge10kb", ROOT / "oapgc/site_effects/species_effects_n50_ge10kb.tsv"),
        ("n50_ge100kb", ROOT / "oapgc/site_effects/species_effects_n50_ge100kb.tsv"),
    ]:
        for r in read_tsv(path):
            for axis, auc, q in [
                ("SNP_distance", "SNP_auc", "SNP_q"),
                ("breakpoints", "BP_auc", "BP_q"),
            ]:
                effects.append({
                    "database": "OAPGC",
                    "species_id": r["cluster"],
                    "species_name": "",
                    "analysis": f"oral_airway_vs_within_{tag}",
                    "axis": axis,
                    "effect_size": r[auc],
                    "effect_direction": "higher_oral_airway" if fnum(r[auc]) > 0.5 else "lower_oral_airway",
                    "p_value": "",
                    "q_value": r[q],
                })

    efields = ["database", "species_id", "species_name", "analysis", "axis",
               "effect_size", "effect_direction", "p_value", "q_value"]
    write_tsv(OUTDIR / "species_effects_unified_long.tsv", effects, efields)


from collections import Counter

if __name__ == "__main__":
    main()
