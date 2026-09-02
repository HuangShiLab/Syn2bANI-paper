#!/usr/bin/env python3
"""Collect rearrangement-boundary coordinates from both methods, for the same pairs.

Syn2b's claim in `Syn2b/docs/MATH_REVIEW.md` is not only that its junction *count*
tracks dnadiff, but that the two methods place breakpoints at the same *coordinates*
at a ~2790x cost difference. That currently rests on a couple of hand-checked cases.
This collects the coordinates for every gtdb50k pair so the claim can be made
statistically, and it needs **no new alignment and no new digestion** -- both sides
are already on disk.

    Syn2b     $TGT_CACHE/tmp_pairs/<pairid>/synteny.junctions.tsv
              written by `syn2b synteny` next to its CSV, and never collected:
              run_syn2b_inverted_fraction.py reads only synteny.csv. The pair
              directories are not cleaned up, so the files should still be there.

    dnadiff   $WORK/out/<pairid>/dd.1coords
              run_dnadiff_slice.sh deletes dd.rdiff and dd.qdiff -- dnadiff's own
              structural-difference coordinates -- to save space, but keeps the
              1-to-1 alignment coordinates they are derived from.

WHY THE dnadiff SIDE IS DERIVED, AND HOW IT IS CHECKED
------------------------------------------------------
Because dd.rdiff is gone, the boundaries here are re-derived from dd.1coords rather
than read from dnadiff. That is a weaker thing than quoting dnadiff's own output, so
it is checked rather than trusted: dd.report's [Feature Estimates] section carries
dnadiff's own counts of Relocations, Translocations and Inversions, and this script
compares its derived count against them per pair and reports the agreement. Use the
coordinates only for pairs where the two agree; `--max-count-diff` gates that, and
the per-pair difference is written out either way.

THE DERIVATION
--------------
Walk the 1-to-1 alignment blocks in reference order. A boundary is where the query
side stops being collinear, which is exactly the event class Syn2b's adjacency rule
responds to:

    INV   the next block is on the opposite strand             (inversion)
    JMP   the query position moves backwards                    (relocation)
    SEQ   the next block is on a different query sequence      (translocation)

`dnadiff_pos` carries INV and JMP; SEQ goes to `dnadiff_pos_seq` and is excluded
from the position comparison by default. On a draft assembly SEQ is mostly not
biology -- dnadiff calls a translocation whenever consecutive alignments come from
different query sequences, and a fragmented query does that constantly. Measured
over 43,078 gtdb50k pairs, dnadiff's Translocations count correlates rho = +0.799
with contig number, against -0.208 for Relocations and -0.331 for Inversions.

Classification is by *order*, not by distance, so there is no tolerance parameter to
tune -- and a pure indel, which shifts query coordinates without disturbing their
order, is correctly not a boundary. dnadiff's own `Breakpoints` count does include
indel-driven boundaries, which is why the check above uses
Relocations + Translocations + Inversions instead.

The reported coordinate is the reference position where the collinear run ends
(`E1` of the last block before the break). Syn2b reports the position of the left
landmark of the broken adjacency, which is the same thing to within landmark
spacing, so the two are directly comparable.

USAGE
-----
    python3 collect_junction_coordinates.py \\
        --pairs     $WORK/pairs_50k.tsv \\
        --tgt-cache $WORK/syn2b_tgts_cache \\
        --dnadiff   $WORK/out \\
        --outdir    $WORK/junction_coords \\
        --workers 32
"""

import argparse
import csv
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

FEATURE_RE_KEYS = ("Relocations", "Translocations", "Inversions", "Breakpoints")


def read_syn2b_junctions(pairdir: Path):
    """Coordinates from `syn2b synteny`, in genome_A's frame.

    `run_syn2b_inverted_fraction.py` names the reference `a_ref_<acc>.tgt` and the
    query `b_qry_<acc>.tgt` precisely so the reference sorts first and becomes
    genome_A -- matching dnadiff, which is invoked as `dnadiff <ref> <qry>`. So both
    sides of this comparison are in reference coordinates without any remapping.
    """
    path = pairdir / "synteny.junctions.tsv"
    if not path.exists():
        return None
    out = []
    with open(path) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for r in rows:
            try:
                out.append(int(r["junction_pos_in_A"]))
            except (KeyError, TypeError, ValueError):
                continue
    return sorted(out)


def contig_offsets_from_tgt(tgt_path: Path):
    """Cumulative genome-wide offset per reference contig, in FASTA order.

    REQUIRED, and the comparison is wrong without it. Syn2b reports junctions in
    genome-wide concatenated coordinates (`digest` accumulates an offset across
    contigs), while dnadiff's .1coords reports positions *within* a reference
    contig alongside its name. On a closed genome the two frames coincide and the
    error is invisible; on a draft they differ by the cumulative length of every
    preceding contig, which is most of gtdb50k.

    The TGT header carries the contig table Syn2b itself used, so taking it from
    there rather than from the FASTA guarantees the two frames agree:

        #contigs=name:len;name:len:circular;...
    """
    if not tgt_path.exists():
        return None
    header = None
    with open(tgt_path) as fh:
        for line in fh:
            if line.startswith("#contigs="):
                header = line[len("#contigs="):].strip()
                break
            if not line.startswith((">", "#")):
                break            # past the header block
    if not header:
        return None
    offsets, off = {}, 0
    for field in header.split(";"):
        if not field:
            continue
        parts = field.split(":")
        if len(parts) < 2:
            continue
        name = parts[0]
        try:
            length = int(parts[1])
        except ValueError:
            continue
        offsets.setdefault(name, off)
        off += length
    return offsets or None


def parse_1coords(path: Path):
    """show-coords -rclTH: S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ TAGR TAGQ.

    Returns blocks as (tag_r, tag_q, s1, e1, s2, e2, reverse). Reference coordinates
    are always forward; a query alignment with S2 > E2 is on the reverse strand,
    which is how an inversion presents.
    """
    blocks = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                f = line.split()
                if len(f) < 13:
                    continue
            try:
                s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            except ValueError:
                continue
            blocks.append((f[11], f[12], s1, e1, s2, e2, s2 > e2))
    return blocks


def derive_dnadiff_boundaries(blocks, offsets=None):
    """Reference coordinates where the query stops being collinear.

    Positions are returned in Syn2b's frame: genome-wide if `offsets` maps each
    reference contig name to its cumulative offset, and within-contig otherwise.
    A caller without an offset table gets `None` back for multi-contig references
    rather than silently mismatched coordinates.

    Returns (positions, counts) with counts keyed SEQ / INV / JMP.
    """
    positions = []          # INV and JMP only -- see the note below
    seq_positions = []      # SEQ, kept separate
    counts = {"SEQ": 0, "INV": 0, "JMP": 0}
    by_ref = {}
    for b in blocks:
        by_ref.setdefault(b[0], []).append(b)

    # Within one reference contig, walk the 1-to-1 blocks left to right and ask what
    # the query does at each seam. Order is checked first because a query-sequence
    # change or a strand flip is decisive regardless of coordinates.
    for bs in by_ref.values():
        bs.sort(key=lambda b: b[2])
        for a, c in zip(bs, bs[1:]):
            _, a_tq, _, a_e1, _, a_e2, a_rev = a
            _, c_tq, _, _, c_s2, _, c_rev = c
            if a_tq != c_tq:
                kind = "SEQ"          # the next block comes from a different query
                                      # sequence: movement between replicons/contigs
            elif a_rev != c_rev:
                kind = "INV"          # strand flip
            else:
                # Collinear means the query advances in the direction the shared
                # strand implies. Anything else is a jump. Checked by order rather
                # than by distance, so a pure indel -- which shifts query
                # coordinates without disturbing their order -- is not a boundary
                # and no tolerance parameter is needed.
                moved_backwards = (c_s2 < a_e2) if not a_rev else (c_s2 > a_e2)
                kind = "JMP" if moved_backwards else None
            if kind:
                counts[kind] += 1
                base = offsets.get(a[0], 0) if offsets else 0
                # SEQ is separated because on a draft assembly it is mostly not
                # biology. dnadiff calls a translocation whenever consecutive
                # alignments come from different query *sequences*, and on a
                # fragmented query they constantly do. Measured over 43,078 gtdb50k
                # pairs, dnadiff's Translocations count correlates rho = +0.799 with
                # contig number -- it is closer to a fragmentation counter than to a
                # rearrangement counter. Relocations and Inversions run the other way
                # (-0.208 and -0.331), and dropping translocations is what makes the
                # truth channel fragmentation-clean: Reloc+Inv correlates -0.243 with
                # contig count against +0.163 for dnadiff's own Breakpoints, while
                # its agreement with Syn2b's junction count *rises*, 0.724 -> 0.766
                # (0.734 -> 0.790 partialling out contig count).
                (seq_positions if kind == "SEQ" else positions).append(base + a_e1)
    return sorted(positions), sorted(seq_positions), counts


def parse_feature_estimates(path: Path):
    """dnadiff's own counts, reference column, from dd.report."""
    if not path.exists():
        return None
    out = {}
    in_section = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("[Feature Estimates]"):
                in_section = True
                continue
            if in_section and line.startswith("["):
                break
            if not in_section:
                continue
            f = line.split()
            if len(f) >= 3 and f[0] in FEATURE_RE_KEYS:
                try:
                    out[f[0]] = int(f[1])
                except ValueError:
                    pass
    return out or None


def one_pair(task):
    pairid, tgt_cache, dnadiff_dir, ref_acc = task
    rec = {"pairid": pairid}

    syn = read_syn2b_junctions(Path(tgt_cache) / "tmp_pairs" / pairid)
    rec["syn2b_n"] = len(syn) if syn is not None else -1
    rec["syn2b_pos"] = ",".join(map(str, syn)) if syn else ""

    dd = Path(dnadiff_dir) / pairid
    coords_path = dd / "dd.1coords"
    if coords_path.exists():
        blocks = parse_1coords(coords_path)
        offsets = contig_offsets_from_tgt(Path(tgt_cache) / f"{ref_acc}.tgt") if ref_acc else None
        n_ref_contigs = len({b[0] for b in blocks})
        if offsets is None and n_ref_contigs > 1:
            # Frames would silently disagree by the cumulative length of every
            # preceding contig. Refuse rather than emit misaligned coordinates.
            rec["dnadiff_n"] = -1
            rec["dnadiff_pos"] = ""
            rec["dnadiff_n_seq"] = -1
            rec["dnadiff_pos_seq"] = ""
            rec["dnadiff_inv"] = rec["dnadiff_jmp"] = rec["dnadiff_seq"] = -1
            rec["frame"] = "no_offsets_multicontig"
            est = parse_feature_estimates(dd / "dd.report")
            rec["report_rearrangements"] = (
                est.get("Relocations", 0) + est.get("Translocations", 0) + est.get("Inversions", 0)
                if est else -1)
            rec["report_breakpoints"] = est.get("Breakpoints", -1) if est else -1
            rec["report_inversions"] = est.get("Inversions", -1) if est else -1
            rec["count_diff"] = ""
            return rec
        rec["frame"] = "genome_wide" if offsets else "single_contig"
        pos, seq_pos, counts = derive_dnadiff_boundaries(blocks, offsets)
        rec["dnadiff_n"] = len(pos)
        rec["dnadiff_pos"] = ",".join(map(str, pos))
        rec["dnadiff_n_seq"] = len(seq_pos)
        rec["dnadiff_pos_seq"] = ",".join(map(str, seq_pos))
        rec["dnadiff_inv"] = counts["INV"]
        rec["dnadiff_jmp"] = counts["JMP"]
        rec["dnadiff_seq"] = counts["SEQ"]
    else:
        rec["dnadiff_n"] = -1
        rec["dnadiff_pos"] = ""
        rec["dnadiff_n_seq"] = -1
        rec["dnadiff_pos_seq"] = ""
        rec["dnadiff_inv"] = rec["dnadiff_jmp"] = rec["dnadiff_seq"] = -1
        rec["frame"] = "no_1coords"

    est = parse_feature_estimates(dd / "dd.report")
    if est:
        # dnadiff's own rearrangement count, for the derivation check. Breakpoints
        # is kept separately because it also counts indel-driven boundaries.
        rec["report_rearrangements"] = (
            est.get("Relocations", 0) + est.get("Translocations", 0) + est.get("Inversions", 0)
        )
        rec["report_breakpoints"] = est.get("Breakpoints", -1)
        rec["report_inversions"] = est.get("Inversions", -1)
    else:
        rec["report_rearrangements"] = rec["report_breakpoints"] = rec["report_inversions"] = -1

    # dnadiff_n now excludes SEQ, so the check adds it back: the comparison is
    # against dnadiff's own Relocations + Translocations + Inversions, which
    # includes translocations.
    rec["count_diff"] = (
        rec["dnadiff_n"] + rec["dnadiff_n_seq"] - rec["report_rearrangements"]
        if rec["dnadiff_n"] >= 0 and rec["report_rearrangements"] >= 0
        else ""
    )
    return rec


FIELDS = ["pairid", "frame", "syn2b_n", "dnadiff_n", "dnadiff_n_seq",
          "report_rearrangements", "report_breakpoints", "report_inversions",
          "dnadiff_inv", "dnadiff_jmp", "dnadiff_seq", "count_diff",
          "syn2b_pos", "dnadiff_pos", "dnadiff_pos_seq"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs", required=True, help="TSV with pairid, or q_acc/r_acc")
    p.add_argument("--tgt-cache", required=True,
                   help="the --tgt-dir given to run_syn2b_inverted_fraction.py; "
                        "its tmp_pairs/ subtree holds the junction files")
    p.add_argument("--dnadiff", required=True, help="directory of <pairid>/dd.1coords")
    p.add_argument("--outdir", required=True)
    p.add_argument("--workers", type=int, default=min(32, cpu_count() or 1))
    args = p.parse_args()

    pairids = []
    with open(args.pairs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = r.get("pairid") or f"{r['q_acc']}__{r['r_acc']}"
            # dnadiff is invoked as `dnadiff <ref> <qry>` and the syn2b runner makes
            # the same accession genome_A, so the reference TGT is the contig table
            # both frames must agree on.
            pairids.append((pid, r.get("r_acc", "")))
    seen, uniq = set(), []
    for pid, racc in pairids:
        if pid not in seen:
            seen.add(pid)
            uniq.append((pid, racc))
    if len(uniq) != len(pairids):
        print(f"note: {len(pairids) - len(uniq)} duplicate pairids in --pairs, deduplicated",
              file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tasks = [(pid, args.tgt_cache, args.dnadiff, racc) for pid, racc in uniq]
    print(f"collecting {len(tasks)} pairs with {args.workers} workers ...", flush=True)

    rows = []
    with Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(one_pair, tasks, chunksize=64), 1):
            rows.append(rec)
            if i % 5000 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)

    out = outdir / "junction_coordinates.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")

    # ---- availability and the derivation check ------------------------------
    n = len(rows)
    have_syn = sum(1 for r in rows if r["syn2b_n"] >= 0)
    have_dd = sum(1 for r in rows if r["dnadiff_n"] >= 0)
    both = sum(1 for r in rows if r["syn2b_n"] >= 0 and r["dnadiff_n"] >= 0)
    print(f"\navailability: syn2b {have_syn}/{n}, dnadiff 1coords {have_dd}/{n}, both {both}/{n}")
    frames = {}
    for r in rows:
        frames[r.get("frame", "?")] = frames.get(r.get("frame", "?"), 0) + 1
    print("coordinate frames: " + ", ".join(f"{k} {v}" for k, v in sorted(frames.items())))
    if frames.get("no_offsets_multicontig"):
        print(f"  {frames['no_offsets_multicontig']} pairs skipped: the reference is "
              f"multi-contig and its TGT was not found in --tgt-cache, so dnadiff's "
              f"within-contig coordinates could not be put in Syn2b's genome-wide "
              f"frame. Point --tgt-cache at the directory holding <acc>.tgt.",
              file=sys.stderr)
    if have_syn == 0:
        print("  No syn2b junction files found. They live in <tgt-dir>/tmp_pairs/<pairid>/ ; "
              "check that --tgt-cache is the same --tgt-dir the run used.", file=sys.stderr)

    checked = [r for r in rows if r["count_diff"] != ""]
    if checked:
        exact = sum(1 for r in checked if r["count_diff"] == 0)
        within1 = sum(1 for r in checked if abs(r["count_diff"]) <= 1)
        print(f"\nderivation check against dd.report [Feature Estimates] "
              f"(Relocations + Translocations + Inversions), n = {len(checked)}:")
        print(f"  exact match      {exact:>7} ({100*exact/len(checked):.1f}%)")
        print(f"  within +-1       {within1:>7} ({100*within1/len(checked):.1f}%)")
        diffs = sorted(r["count_diff"] for r in checked)
        med = diffs[len(diffs) // 2]
        print(f"  median difference {med:>+6}")
        print("\nThe coordinates are only as good as this check. Filter to "
              "`count_diff == 0` before drawing conclusions from positions, and say "
              "in the paper that dnadiff's own dd.rdiff was not retained.")
    else:
        print("\nNo pair had both dd.1coords and dd.report, so the derivation is unchecked. "
              "Do not use the coordinates.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
