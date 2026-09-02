#!/usr/bin/env python3
"""Assign E. coli O157:H7 lineage based on closest reference by anchor adjacency."""
import csv
import re
from pathlib import Path
from collections import defaultdict

WORKDIR = Path(__file__).resolve().parent

# Reference strain accessions (base nucleotide_acc) and their Fitzgerald lineage labels
REFS = {
    "NC_002695": "Ia",   # Sakai
    "CP018252": "Ic",    # strain 9000
    "NC_013008": "I/II", # TW14359
    "CP015832": "II",    # strain 180
}

REF_NAMES = {
    "NC_002695": "sakai",
    "CP018252": "9000",
    "NC_013008": "tw14359",
    "CP015832": "180",
}


def normalize_acc(s):
    """Map versioned RefSeq/GenBank accessions to the base nucleotide accession."""
    s = re.sub(r"^NZ_", "", s)
    return s.split(".")[0]


def load_triangle():
    pairs = defaultdict(dict)
    with (WORKDIR / "triangle.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            q = normalize_acc(row["query"])
            r = normalize_acc(row["reference"])
            key = tuple(sorted([q, r]))
            # Store both orientations so lookup is easy
            pairs[q][r] = {
                "ani": float(row["ani"]),
                "breakpoint_count": int(row["breakpoint_count"]),
                "anchor_adjacency": float(row["anchor_adjacency"]),
            }
            pairs[r][q] = pairs[q][r]
    return pairs


def load_metadata():
    with (WORKDIR / "metadata.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {row["nucleotide_acc"]: row for row in reader}
    return rows


def assign_lineage(acc, pairs):
    """Return dict of ref metrics and the chosen lineage for a single genome."""
    # Reference strains are anchored to their own known lineage
    if acc in REFS:
        assigned = REFS[acc]
    else:
        assigned = None

    metrics = {}
    candidates = []
    for ref_acc, ref_label in REFS.items():
        pair = pairs.get(acc, {}).get(ref_acc)
        if pair is None:
            metrics[f"ref_{REF_NAMES[ref_acc]}_ani"] = None
            metrics[f"ref_{REF_NAMES[ref_acc]}_breakpoints"] = None
            metrics[f"ref_{REF_NAMES[ref_acc]}_anchor_adjacency"] = None
            continue
        metrics[f"ref_{REF_NAMES[ref_acc]}_ani"] = pair["ani"]
        metrics[f"ref_{REF_NAMES[ref_acc]}_breakpoints"] = pair["breakpoint_count"]
        metrics[f"ref_{REF_NAMES[ref_acc]}_anchor_adjacency"] = pair["anchor_adjacency"]
        candidates.append({
            "ref_acc": ref_acc,
            "lineage": ref_label,
            "anchor_adjacency": pair["anchor_adjacency"],
            "ani": pair["ani"],
            "breakpoints": pair["breakpoint_count"],
        })

    if assigned is None:
        if not candidates:
            assigned = "unassigned"
        else:
            # Highest anchor_adjacency, then highest ANI, then fewest breakpoints
            best = max(
                candidates,
                key=lambda c: (c["anchor_adjacency"], c["ani"], -c["breakpoints"]),
            )
            assigned = best["lineage"]
    return assigned, metrics


def categorize_host(source):
    """Categorize isolation_source into human / bovine / other/unknown."""
    if not source:
        return "other/unknown"
    s = source.lower()
    human_terms = ("human", "clinical", "patient", "stool", "blood", "urine")
    bovine_terms = ("bovine", "cattle", "cow", "beef", "fecal", "faeces", "faecal")
    if any(t in s for t in human_terms):
        return "human"
    if any(t in s for t in bovine_terms):
        return "bovine"
    return "other/unknown"


def main():
    pairs = load_triangle()
    meta = load_metadata()

    output_rows = []
    for acc in sorted(meta.keys()):
        assigned, metrics = assign_lineage(acc, pairs)
        row = {
            "nucleotide_acc": acc,
            "strain": meta[acc].get("strain", ""),
            "assigned_lineage": assigned,
        }
        row.update(metrics)
        output_rows.append(row)

    fieldnames = [
        "nucleotide_acc",
        "strain",
        "assigned_lineage",
        "ref_sakai_ani",
        "ref_9000_ani",
        "ref_tw14359_ani",
        "ref_180_ani",
        "ref_sakai_breakpoints",
        "ref_9000_breakpoints",
        "ref_tw14359_breakpoints",
        "ref_180_breakpoints",
        "ref_sakai_anchor_adjacency",
        "ref_9000_anchor_adjacency",
        "ref_tw14359_anchor_adjacency",
        "ref_180_anchor_adjacency",
    ]

    out_path = WORKDIR / "lineage_assignments.tsv"
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {out_path} with {len(output_rows)} rows")

    # Merge lineage + host_category into metadata and write metadata_with_lineage.tsv
    lineage_by_acc = {r["nucleotide_acc"]: r["assigned_lineage"] for r in output_rows}
    meta_out_path = WORKDIR / "metadata_with_lineage.tsv"
    meta_fieldnames = list(meta.values())[0].keys()
    if "assigned_lineage" not in meta_fieldnames:
        meta_fieldnames = list(meta_fieldnames) + ["assigned_lineage"]
    if "host_category" not in meta_fieldnames:
        meta_fieldnames = list(meta_fieldnames) + ["host_category"]

    with meta_out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=meta_fieldnames, delimiter="\t")
        writer.writeheader()
        for acc in sorted(meta.keys()):
            row = dict(meta[acc])
            row["assigned_lineage"] = lineage_by_acc.get(acc, "unassigned")
            row["host_category"] = categorize_host(row.get("isolation_source", ""))
            writer.writerow(row)
    print(f"Wrote {meta_out_path} with {len(meta)} rows")

    # Quick sanity report
    from collections import Counter
    counts = Counter(r["assigned_lineage"] for r in output_rows)
    print("Lineage counts:")
    for lineage, cnt in sorted(counts.items()):
        print(f"  {lineage}: {cnt}")


if __name__ == "__main__":
    main()
