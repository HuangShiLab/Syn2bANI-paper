#!/usr/bin/env python3
"""Parse syn2bani struct PAF outputs for top Syntracker discordant pairs.

Produces a summary TSV of chain-level structural variation statistics:
- chain count
- inversion count (adjacent chains on opposite strands)
- translocation count (adjacent chains on different reference contigs)
- indel count and total indel bases (gaps between collinear chains)
"""

import argparse
from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_STRUCT_DIR = ROOT / "results" / "syntracker_validation" / "struct_top_cases"
DEFAULT_OUT_TSV = ROOT / "results" / "syntracker_validation" / "struct_top_cases_summary.tsv"


def parse_paf(path):
    """Return list of chain dicts from a syn2bani struct PAF file."""
    chains = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            chains.append({
                "qname": cols[0],
                "qlen": int(cols[1]),
                "qstart": int(cols[2]),
                "qend": int(cols[3]),
                "strand": cols[4],
                "rname": cols[5],
                "rlen": int(cols[6]),
                "rstart": int(cols[7]),
                "rend": int(cols[8]),
                "matches": int(cols[9]),
                "block_len": int(cols[10]),
                "mapq": int(cols[11]),
            })
    return chains


def classify_svs(chains, min_gap=100):
    """Classify SVs from collinear chains sorted by query coordinate."""
    if not chains:
        return {}

    # Sort by query start; tie-break by query end for stable ordering.
    ordered = sorted(chains, key=lambda c: (c["qstart"], c["qend"]))

    n_chains = len(ordered)
    n_inversions = 0
    n_translocations = 0
    n_indels = 0
    total_indel_bp = 0
    indel_sizes = []

    for prev, nxt in zip(ordered, ordered[1:]):
        # Strand / orientation change -> inversion
        if prev["strand"] != nxt["strand"]:
            n_inversions += 1
            continue

        # Reference contig change -> translocation
        if prev["rname"] != nxt["rname"]:
            n_translocations += 1
            continue

        # Collinear pair: infer indel from gap difference
        query_gap = nxt["qstart"] - prev["qend"]

        if prev["strand"] == "+":
            ref_gap = nxt["rstart"] - prev["rend"]
        else:
            # Both are on minus strand; reference coordinates decrease along the chain order
            ref_gap = prev["rstart"] - nxt["rend"]

        size = abs(ref_gap - query_gap)
        if size >= min_gap:
            n_indels += 1
            total_indel_bp += size
            indel_sizes.append(size)

    # Genome-level summaries
    query_contigs = {c["qname"] for c in ordered}
    ref_contigs = {c["rname"] for c in ordered}
    query_len = sum({c["qname"]: c["qlen"] for c in ordered}.values())
    ref_len = sum({c["rname"]: c["rlen"] for c in ordered}.values())
    aligned_query_bp = sum(c["qend"] - c["qstart"] for c in ordered)
    aligned_ref_bp = sum(c["rend"] - c["rstart"] for c in ordered)

    return {
        "n_chains": n_chains,
        "n_query_contigs": len(query_contigs),
        "n_ref_contigs": len(ref_contigs),
        "query_len": query_len,
        "ref_len": ref_len,
        "aligned_query_bp": aligned_query_bp,
        "aligned_ref_bp": aligned_ref_bp,
        "n_inversions": n_inversions,
        "n_translocations": n_translocations,
        "n_indels": n_indels,
        "total_indel_bp": total_indel_bp,
        "mean_indel_bp": round(sum(indel_sizes) / len(indel_sizes), 1) if indel_sizes else 0.0,
        "max_indel_bp": max(indel_sizes) if indel_sizes else 0,
    }


def parse_ani(path):
    """Read the verbose ani TSV and return key metrics as a dict."""
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}
    row = df.iloc[0]
    return {
        "ani": row.get("ani", row.get("ani_uniform", None)),
        "synteny_score": row.get("synteny_score", None),
        "breakpoint_count": row.get("breakpoint_count", None),
        "af_query": row.get("af_query", None),
        "af_reference": row.get("af_reference", None),
        "n_anchors": row.get("n_anchors", None),
        "n_tags": row.get("n_tags", None),
    }


def main():
    parser = argparse.ArgumentParser(description="Parse syn2bani struct PAF outputs and summarize SV statistics.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_STRUCT_DIR,
                        help="Directory containing struct .tsv PAF outputs.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT_TSV,
                        help="Output summary TSV path.")
    args = parser.parse_args()

    rows = []
    for paf_path in sorted(args.input_dir.glob("*.tsv")):
        if paf_path.name.endswith("_ani.tsv"):
            continue
        case = paf_path.stem
        ani_path = paf_path.with_suffix("").with_name(case + "_ani.tsv")

        chains = parse_paf(paf_path)
        sv = classify_svs(chains)
        ani = parse_ani(ani_path) if ani_path.exists() else {}

        row = {
            "case": case,
            **sv,
            **ani,
        }
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary = summary[[
        "case",
        "ani",
        "synteny_score",
        "breakpoint_count",
        "n_chains",
        "n_query_contigs",
        "n_ref_contigs",
        "query_len",
        "ref_len",
        "aligned_query_bp",
        "aligned_ref_bp",
        "n_inversions",
        "n_translocations",
        "n_indels",
        "total_indel_bp",
        "mean_indel_bp",
        "max_indel_bp",
        "af_query",
        "af_reference",
        "n_anchors",
        "n_tags",
    ]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote summary for {len(summary)} cases to {args.output}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
