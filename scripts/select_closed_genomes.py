#!/usr/bin/env python3
"""Select closed genomes of target species from GTDB metadata, and emit the
accession list plus the pair file for an all-vs-all structural comparison.

The inverted-fraction channel needs contiguous assemblies (see
`annotate_closed_inversion_pairs.py` for the measurement). This picks genomes
that qualify, from metadata already on disk, so the accessions are verified
rather than recalled. It also reports which ones are already in the local genome
directory, so only the remainder has to be downloaded.

Species names must be GTDB `s__` names without the prefix, e.g.
"Streptococcus pneumoniae". GTDB sometimes splits an NCBI species into
`Foo bar` and `Foo bar_A`; pass --genus instead to catch every one of them, or
use --list-species to see what exists before choosing.

Usage:
    # what does GTDB actually call things?
    python3 scripts/select_closed_genomes.py --metadata M --list-species Streptococcus

    # select, and see what still needs downloading
    python3 scripts/select_closed_genomes.py \
        --metadata   /lustre1/g/aos_shihuang/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
        --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
        --species "Streptococcus pneumoniae" "Salmonella enterica" \
        --outdir results/closed_inversions
"""

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

WANT = [
    "accession",
    "gtdb_taxonomy",
    "contig_count",
    "checkm_completeness",
    "checkm_contamination",
    "genome_size",
    "n50_contigs",
    "ncbi_assembly_level",
    "ncbi_genbank_assembly_accession",
    "ncbi_strain_identifiers",
]


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    have = [c for c in WANT if c in df.columns]
    for req in ("accession", "gtdb_taxonomy", "contig_count"):
        if req not in have:
            sys.exit(f"metadata lacks required column '{req}'")
    df = df[have].copy()
    df["acc"] = df.accession.str.replace(r"^(RS_|GB_)", "", regex=True)
    tax = df.gtdb_taxonomy.str.split(";", expand=True)
    df["genus"] = tax[5].str.replace("^g__", "", regex=True) if tax.shape[1] > 5 else None
    df["species"] = tax[6].str.replace("^s__", "", regex=True) if tax.shape[1] > 6 else None
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--genome-dir", default=None,
                    help="if given, report which accessions are already present")
    ap.add_argument("--species", nargs="*", default=[])
    ap.add_argument("--genus", nargs="*", default=[])
    ap.add_argument("--list-species", metavar="GENUS",
                    help="print GTDB species names in this genus and exit")
    ap.add_argument("--max-contigs", type=int, default=3)
    ap.add_argument("--min-completeness", type=float, default=98.0)
    ap.add_argument("--max-contamination", type=float, default=2.0)
    ap.add_argument("--max-per-species", type=int, default=200,
                    help="cap per species; all-vs-all is quadratic (default 200 -> 19,900 pairs)")
    ap.add_argument("--outdir", default="results/closed_inversions")
    args = ap.parse_args()

    df = load(Path(args.metadata))

    if args.list_species:
        g = df[df.genus == args.list_species]
        if not len(g):
            sys.exit(f"no GTDB genus '{args.list_species}'")
        t = (g.groupby("species")
               .agg(total=("acc", "size"),
                    closed=("contig_count", lambda x: (x <= args.max_contigs).sum()))
               .sort_values("closed", ascending=False))
        print(t.to_string())
        return

    if not args.species and not args.genus:
        sys.exit("give --species and/or --genus (or --list-species to explore)")

    sel = df[df.species.isin(args.species) | df.genus.isin(args.genus)].copy()
    if not len(sel):
        sys.exit("nothing matched; check names with --list-species")

    keep = sel.contig_count <= args.max_contigs
    if "checkm_completeness" in sel:
        keep &= sel.checkm_completeness >= args.min_completeness
    if "checkm_contamination" in sel:
        keep &= sel.checkm_contamination <= args.max_contamination
    sel = sel[keep].sort_values(["species", "contig_count", "checkm_completeness"],
                               ascending=[True, True, False])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    present = set()
    if args.genome_dir:
        gd = Path(args.genome_dir)
        present = {p.stem for p in gd.glob("*.fna")} | {p.stem for p in gd.glob("*.fasta")}

    print(f"{'species':<46} {'closed':>7} {'kept':>6} {'local':>6} {'to fetch':>9}")
    frames = []
    for sp, g in sel.groupby("species"):
        g = g.head(args.max_per_species)
        loc = g.acc.isin(present).sum() if present else 0
        print(f"{str(sp):<46} {len(sel[sel.species == sp]):>7} {len(g):>6} "
              f"{loc:>6} {len(g) - loc:>9}")
        frames.append(g)
    sel = pd.concat(frames) if frames else sel

    cols = [c for c in ["acc", "species", "genus", "contig_count", "n50_contigs",
                        "genome_size", "checkm_completeness", "checkm_contamination",
                        "ncbi_assembly_level", "ncbi_genbank_assembly_accession",
                        "ncbi_strain_identifiers"] if c in sel]
    sel[cols].to_csv(outdir / "genomes.tsv", sep="\t", index=False)
    (outdir / "accessions.txt").write_text("\n".join(sel.acc) + "\n")
    if present:
        missing = sel[~sel.acc.isin(present)]
        (outdir / "accessions_to_download.txt").write_text(
            "\n".join(missing.acc) + ("\n" if len(missing) else "")
        )

    rows = []
    for sp, g in sel.groupby("species"):
        for q, r in itertools.combinations(sorted(g.acc), 2):
            rows.append({"pairid": f"{q}__{r}", "q_acc": q, "r_acc": r, "species": sp})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(outdir / "pairs_all_vs_all.tsv", sep="\t", index=False)

    print(f"\n{len(sel)} genomes, {len(pairs)} within-species pairs")
    print(f"  {outdir}/genomes.tsv")
    print(f"  {outdir}/accessions.txt")
    if present:
        print(f"  {outdir}/accessions_to_download.txt  ({len(sel) - len(sel[sel.acc.isin(present)])} missing)")
    print(f"  {outdir}/pairs_all_vs_all.tsv")


if __name__ == "__main__":
    main()
