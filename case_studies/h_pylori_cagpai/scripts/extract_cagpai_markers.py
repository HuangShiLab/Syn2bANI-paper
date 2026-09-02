#!/usr/bin/env python3
"""Extract cagPAI marker sequences from H. pylori 26695 GenBank annotation.

The cagPAI in H. pylori 26695 spans locus tags HP0520-HP0547. This script
parses the RefSeq GenBank file (full annotation is in hp26695.gbff; the
hp26695.gb file only contains source features) and writes one FASTA entry
per marker CDS.
"""

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO


def parse_locus_tag_range(tag: str):
    """Return numeric part of an old_locus_tag like 'HP0520' or 'HP_0520'."""
    m = re.match(r"HP[_-]?(\d+)", tag, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def in_cagpai_range(tag: str, low: int = 520, high: int = 547) -> bool:
    """Check whether old_locus_tag numeric part falls in HP0520-HP0547."""
    num = parse_locus_tag_range(tag)
    if num is None:
        return False
    return low <= num <= high


def main():
    parser = argparse.ArgumentParser(description="Extract cagPAI marker FASTA.")
    parser.add_argument(
        "--genbank",
        default="/Users/macstudio/Downloads/Syn2bANI-paper/data/cagpai_pilot/hp26695.gbff",
        help="Input GenBank file with CDS annotations (full annotation is in .gbff).",
    )
    parser.add_argument(
        "--out",
        default="cagpai_markers.fna",
        help="Output FASTA file.",
    )
    args = parser.parse_args()

    gb_path = Path(args.genbank)
    out_path = Path(args.out)

    record = SeqIO.read(gb_path, "genbank")

    markers = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        old_tags = feat.qualifiers.get("old_locus_tag", [])
        matched_tag = None
        for tag in old_tags:
            if in_cagpai_range(tag):
                matched_tag = tag
                break
        if matched_tag is None:
            continue

        gene = feat.qualifiers.get("gene", [""])[0]
        locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
        product = feat.qualifiers.get("product", [""])[0]
        seq = feat.extract(record.seq)

        # Normalise tag to HPXXXX form.
        num = parse_locus_tag_range(matched_tag)
        norm_tag = f"HP{num:04d}"
        header = f"{norm_tag}|{gene}" if gene else norm_tag
        markers.append((header, str(seq), locus_tag, product))

    # Sort by numeric locus tag to keep biological order.
    markers.sort(key=lambda x: int(x[0].split("|")[0][2:]))

    with out_path.open("w") as fh:
        for header, seq, locus_tag, product in markers:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")

    print(f"Wrote {len(markers)} markers to {out_path}")
    for header, seq, locus_tag, product in markers:
        print(f"  {header:15s} len={len(seq):5d}  ref_locus_tag={locus_tag:12s}  {product}")


if __name__ == "__main__":
    main()
