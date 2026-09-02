#!/usr/bin/env python3
"""Call cagPAI status (complete/partial/empty) per H. pylori genome.

Aligns each query genome against cagPAI marker sequences with minimap2,
computes per-marker query coverage and identity, and classifies the island
status using thresholds:
    present : coverage >= 0.8 AND identity >= 0.8
    complete: fraction_present >= 0.85
    empty   : fraction_present <= 0.15
    partial : otherwise
"""

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


MARKER_RE = re.compile(r"^(HP\d{4})\|?(.*)?$")


def parse_marker_name(raw: str):
    """Return (hp_tag, gene_or_extra) from a FASTA header like 'HP0520|cag1'."""
    m = MARKER_RE.match(raw)
    if m:
        return m.group(1), (m.group(2) or "")
    return raw, ""


def parse_paf_line(line: str):
    """Return dict of PAF fields and tags."""
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 12:
        return None
    record = {
        "qname": cols[0],
        "qlen": int(cols[1]),
        "qstart": int(cols[2]),
        "qend": int(cols[3]),
        "strand": cols[4],
        "tname": cols[5],
        "tlen": int(cols[6]),
        "tstart": int(cols[7]),
        "tend": int(cols[8]),
        "matches": int(cols[9]),
        "block_len": int(cols[10]),
        "mapq": int(cols[11]),
    }
    for tag in cols[12:]:
        if tag.startswith("de:f:"):
            try:
                record["de"] = float(tag[5:])
            except ValueError:
                pass
        elif tag.startswith("nm:i:"):
            try:
                record["nm"] = int(tag[5:])
            except ValueError:
                pass
    return record


def identity_from_record(rec: dict) -> float:
    """Best available identity estimate from a PAF record.

    minimap2's de:f tag is gap-compressed per-base sequence divergence
    (1 - identity), so identity = 1 - de.
    """
    if "de" in rec:
        return 1.0 - rec["de"]
    if "nm" in rec and rec["block_len"] > 0:
        return 1.0 - (rec["nm"] / rec["block_len"])
    if rec["block_len"] > 0:
        return rec["matches"] / rec["block_len"]
    return 0.0


def call_markers(paf_lines, marker_lengths: dict, cov_thr: float = 0.8, id_thr: float = 0.8):
    """Return dict marker -> {coverage, identity, present} from PAF records."""
    best = {}  # marker -> record with most matches
    for line in paf_lines:
        line = line.strip()
        if not line or line.startswith("@"):
            continue
        rec = parse_paf_line(line)
        if rec is None:
            continue
        marker = rec["tname"]
        if marker not in best or rec["matches"] > best[marker]["matches"]:
            best[marker] = rec

    result = {}
    for marker, tlen in marker_lengths.items():
        if marker in best:
            rec = best[marker]
            coverage = (rec["tend"] - rec["tstart"]) / tlen
            identity = identity_from_record(rec)
            result[marker] = {
                "coverage": coverage,
                "identity": identity,
                "present": coverage >= cov_thr and identity >= id_thr,
            }
        else:
            result[marker] = {"coverage": 0.0, "identity": 0.0, "present": False}
    return result


def load_marker_lengths(fasta: Path) -> dict:
    """Return {header: length} for a FASTA file."""
    lengths = {}
    name = None
    seq = []
    with fasta.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = len("".join(seq))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            lengths[name] = len("".join(seq))
    return lengths


def classify(fraction_present: float, complete_thr: float = 0.85, empty_thr: float = 0.15):
    if fraction_present >= complete_thr:
        return "complete"
    if fraction_present <= empty_thr:
        return "empty"
    return "partial"


def call_genome(query_fasta: Path, markers_fasta: Path, minimap2="minimap2",
                cov_thr=0.8, id_thr=0.8, complete_thr=0.85, empty_thr=0.15):
    """Run minimap2 and return status dict for one genome."""
    marker_lengths = load_marker_lengths(markers_fasta)
    marker_order = list(marker_lengths.keys())

    cmd = [
        minimap2,
        "-cx", "asm20",
        str(markers_fasta),
        str(query_fasta),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed for {query_fasta}: {proc.stderr}"
        )

    calls = call_markers(
        proc.stdout.splitlines(), marker_lengths, cov_thr=cov_thr, id_thr=id_thr
    )

    present = [m for m in marker_order if calls[m]["present"]]
    missing = [m for m in marker_order if not calls[m]["present"]]
    n_markers = len(marker_order)
    n_present = len(present)
    fraction = n_present / n_markers if n_markers else 0.0
    status = classify(fraction, complete_thr=complete_thr, empty_thr=empty_thr)

    return {
        "genome": query_fasta.stem,
        "n_markers": n_markers,
        "n_present": n_present,
        "fraction_present": fraction,
        "status": status,
        "missing_markers": ",".join(missing) if missing else "",
    }


def main():
    parser = argparse.ArgumentParser(description="Call cagPAI status per genome.")
    parser.add_argument(
        "queries",
        nargs="+",
        help="Query genome FASTA file(s) or directories containing .fna files.",
    )
    parser.add_argument(
        "--markers",
        default="cagpai_markers.fna",
        help="cagPAI marker FASTA.",
    )
    parser.add_argument(
        "--out",
        default="cagpai_states.tsv",
        help="Output TSV.",
    )
    parser.add_argument(
        "--minimap2",
        default="minimap2",
        help="Path to minimap2 executable.",
    )
    parser.add_argument(
        "--coverage-thr",
        type=float,
        default=0.8,
        help="Minimum marker coverage to count as present.",
    )
    parser.add_argument(
        "--identity-thr",
        type=float,
        default=0.8,
        help="Minimum marker identity to count as present.",
    )
    parser.add_argument(
        "--complete-thr",
        type=float,
        default=0.85,
        help="Fraction-present threshold for 'complete'.",
    )
    parser.add_argument(
        "--empty-thr",
        type=float,
        default=0.15,
        help="Fraction-present threshold for 'empty'.",
    )
    args = parser.parse_args()

    markers_fasta = Path(args.markers)
    if not markers_fasta.exists():
        sys.exit(f"Marker FASTA not found: {markers_fasta}")

    # Resolve inputs: files are used directly; directories are scanned for .fna.
    query_files = []
    for q in args.queries:
        p = Path(q)
        if p.is_dir():
            query_files.extend(sorted(p.glob("*.fna")))
        else:
            query_files.append(p)

    if not query_files:
        sys.exit("No query FASTA files found.")

    rows = []
    for qf in query_files:
        row = call_genome(
            qf,
            markers_fasta,
            minimap2=args.minimap2,
            cov_thr=args.coverage_thr,
            id_thr=args.identity_thr,
            complete_thr=args.complete_thr,
            empty_thr=args.empty_thr,
        )
        rows.append(row)
        print(
            f"{row['genome']:40s} {row['status']:8s} "
            f"{row['n_present']:3d}/{row['n_markers']:3d} "
            f"({row['fraction_present']:.3f})"
        )

    out_path = Path(args.out)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "genome",
                "n_markers",
                "n_present",
                "fraction_present",
                "status",
                "missing_markers",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} records to {out_path}")


if __name__ == "__main__":
    main()
