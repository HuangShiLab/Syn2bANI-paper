#!/usr/bin/env python3
"""check_discordant_genomes.py — per-genome forensic review of the genomes
where the 28-marker cagPAI panel (cagpai_complete) disagrees with the
collaborator's independent cagA genotyping (cagA_detected / cagA_result).

For each discordant genome (FASTA on the read-only external drive) we:
  a. compute assembly stats (n_contigs, total_len, N50);
  b. reproduce the panel's own marker call: minimap2 -cx asm20, genome as
     query, markers as target, per-marker best record by matches, marker-side
     coverage, identity = 1 - de:f (gap-compressed), present if >= 0.8/0.8
     (exactly the method of cagpai_status/call_cagpai_status.py);
  c. cagA-specific: presence decided by a k-mer test — the fraction of
     HP0547|cagA marker positions covered by >=1 shared exact 19-mer with
     the genome (both strands). This is robust to allele divergence (a
     full-length East-Asian cagA at ~89% identity still spans the whole
     gene) and immune to the gappy forced alignments that inflate
     minimap2's gap-compressed de:f identity on degraded loci. present if
     pos_cov >= 0.55, or >= 0.30 with an asm20 alignment spanning >= 75%
     of the gene; remnant if >= 0.08; absent otherwise. Panel-style asm20
     metrics (cov/idn/contig/pos) and a strict -cx asm5 full-length pass
     are reported alongside;
  d. cagPAI locus coverage: the 26695 reference window +/-100 kb flanks
     chopped into 10 kb chunks, each mapped independently with -cx asm20
     (chunk as target, genome as query). Chunking prevents a high-scoring
     spurious gappy chain from suppressing genuine small alignments, and
     the nmatch/alen >= 0.7 per-alignment filter rejects forced
     alignments of non-colinear sequence (which can cover 100+ kb with a
     misleading de:f ~ 0.05). A window with zero coverage but covered
     flanks is a clean (biological) deletion; zero coverage of window AND
     flanks means the locus region is not assembled (draft gap);
  e. cross-check the cag_sv_list column against SV BED records overlapping
     the cagPAI window (struct_vs_26695_filtered/<GCA>.vs_hp26695.bed).

Verdicts (see REVIEW.md): assembly_gap | panel_missed |
      true_biological_cagA_loss | disrupted_cagA | other,
with the assembly_gap notes split into 'draft gap' vs 'clean deletion'.

The external drive is treated as read-only: all temp files go to a local
temp dir under --outdir; outputs stay inside the repo.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile

# cagPAI locus in H. pylori 26695 (NC_000915.1), 1-based inclusive
CAG_S, CAG_E = 547_328, 583_481
FLANK = 100_000
WIN_LEN = CAG_E - CAG_S + 1

MINIMAP2 = "/opt/homebrew/bin/minimap2"

MARKERS = [
    "HP0520|cag1", "HP0521", "HP0522|cag3", "HP0523|cag4", "HP0524|cag5",
    "HP0525|virB11", "HP0526|cagZ", "HP0527|cagY", "HP0528|cagX",
    "HP0529|cagW", "HP0530|cagV", "HP0531|cagU", "HP0532|cagT",
    "HP0534|cagS", "HP0535|cagQ", "HP0536|cagP", "HP0537|cagM",
    "HP0538|cagN", "HP0539|cagL", "HP0540|cagI", "HP0541", "HP0542|cagG",
    "HP0543|cagF", "HP0544|cagE", "HP0545|cagD", "HP0546|cagC",
    "HP0546|cagB", "HP0547|cagA",
]
CAGA = "HP0547|cagA"


# ---------------------------------------------------------------- fasta io --
def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def assembly_stats(seqs):
    lens = sorted((len(s) for s in seqs.values()), reverse=True)
    total = sum(lens)
    half, acc, n50 = total / 2, 0, 0
    for L in lens:
        acc += L
        if acc >= half:
            n50 = L
            break
    return len(lens), total, n50


# --------------------------------------------------------------- minimap2 --
def minimap_paf(target_fasta, query_fasta, preset="asm20"):
    """minimap2 -cx <preset>; returns PAF lines (base-level alignment done,
    so de:f tags are available — required to reproduce the panel)."""
    cmd = [MINIMAP2, "-cx", preset, "--secondary=yes", "-t", "4",
           target_fasta, query_fasta]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode not in (0, 1):
        raise RuntimeError(f"minimap2 failed: {' '.join(cmd)}\n{out.stderr[-500:]}")
    return [l for l in out.stdout.splitlines() if l and not l.startswith("@")]


def parse_paf(lines):
    recs = []
    for l in lines:
        f = l.split("\t")
        if len(f) < 12:
            continue
        de = next((t for t in f[12:] if t.startswith("de:f:")), None)
        recs.append(dict(
            qname=f[0], qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]),
            tname=f[5], tlen=int(f[6]), tstart=int(f[7]), tend=int(f[8]),
            nmatch=int(f[9]), alen=int(f[10]),
            idn=(1 - float(de[5:])) if de else
                (int(f[9]) / int(f[10]) if int(f[10]) else 0.0)))
    return recs


def union_cov(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    tot, cs, ce = 0, intervals[0][0], intervals[0][1]
    for s, e in intervals[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + ce - cs


def clip(iv, a, b):
    """clip intervals to [a,b), dropping empties"""
    out = []
    for s, e in iv:
        s, e = max(s, a), min(e, b)
        if e > s:
            out.append((s, e))
    return out


# ------------------------------------------------------------------ checks --
def marker_recheck_panel_exact(markers_fasta, genome_fasta):
    """Reproduce call_cagpai_status.py: genome = query, markers = target,
    -cx asm20, best record per marker by matches, coverage on the marker
    (target) side, identity = 1 - de:f, present at >= 0.8/0.8."""
    paf = parse_paf(minimap_paf(markers_fasta, genome_fasta, "asm20"))
    best = {}
    for r in paf:
        # markers are the target in this orientation
        if r["tname"] not in MARKERS:
            continue
        m = r["tname"]
        if m not in best or r["nmatch"] > best[m]["nmatch"]:
            best[m] = r
    stats = {}
    for m in MARKERS:
        r = best.get(m)
        if r is None:
            stats[m] = dict(cov=0.0, idn=0.0, present=False, contig="")
        else:
            cov = (r["tend"] - r["tstart"]) / r["tlen"]
            stats[m] = dict(cov=cov, idn=r["idn"],
                            present=cov >= 0.8 and r["idn"] >= 0.8,
                            contig=r["qname"],
                            tpos=f'{r["tstart"]}-{r["tend"]}')
    n_present = sum(1 for m in MARKERS if stats[m]["present"])
    return n_present, stats


def genome_kmer_set(genome_fasta, k=19):
    km = set()
    with open(genome_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            s = line.strip().upper()
            for i in range(len(s) - k + 1):
                x = s[i:i + k]
                km.add(x)
                km.add(x.translate(str.maketrans("ACGT", "TGCA"))[::-1])
    return km


def kmer_position_coverage(query_seq, genome_km, k=19):
    """Fraction of query positions lying within >=1 shared exact k-mer.
    Robust to divergence level (a full-length cagA at 89% identity still
    spans the whole gene) unlike raw k-mer counts, and immune to the
    gappy forced-alignment artifacts that inflate gap-compressed identity.
    Empirics (k=19): full Western cagA ~0.75; full East-Asian cagA
    ~0.34; ~50% remnant ~0.20; absent 0.00."""
    s = query_seq.upper()
    n = len(s)
    if n < k:
        return 0.0
    cov = bytearray(n)
    for i in range(n - k + 1):
        if s[i:i + k] in genome_km:
            cov[i:i + k] = b"\x01" * k
    return sum(cov) / n


def caga_check(caga_fasta, genome_fasta, caga_seq):
    """cagA presence decided by the 19-mer fraction of the marker sequence
    shared with the genome (k-mer sharing is robust to the heavily gappy
    forced alignments that inflate minimap2's gap-compressed de:f identity
    on degraded/remnant loci). Minimap2 -cx asm20 (marker as target) and a
    strict -cx asm5 full-length pass provide contig/position/coverage
    detail. Classes: present (>=0.40) / remnant (0.05-0.40) / absent."""
    out = dict(cls="absent", idn="", cov="", contig="", pos="",
               strict_idn="", strict_cov="", kmer="")
    pos_cov = kmer_position_coverage(caga_seq, genome_kmer_set(genome_fasta))
    out["kmer"] = f"{pos_cov:.3f}"
    # panel-style asm20 detail
    cov = 0.0
    paf = parse_paf(minimap_paf(caga_fasta, genome_fasta, "asm20"))
    hits = [r for r in paf if r["tname"] == CAGA]
    if hits:
        b = max(hits, key=lambda r: r["nmatch"])
        cov = (b["tend"] - b["tstart"]) / b["tlen"]
        out.update(idn=f"{b['idn']:.3f}", cov=f"{cov:.3f}",
                   contig=b["qname"],
                   pos=f'{min(b["qstart"], b["qend"])}-{max(b["qstart"], b["qend"])}')
    # strict full-length confirmation
    paf5 = parse_paf(minimap_paf(genome_fasta, caga_fasta, "asm5"))
    hits5 = [r for r in paf5 if r["qname"] == CAGA]
    if hits5:
        b5 = max(hits5, key=lambda r: r["nmatch"])
        cov5 = (b5["qend"] - b5["qstart"]) / b5["qlen"]
        out["strict_idn"] = f"{b5['nmatch'] / b5['alen']:.3f}"
        out["strict_cov"] = f"{cov5:.3f}"
    if pos_cov >= 0.55 or (pos_cov >= 0.30 and cov >= 0.75):
        out["cls"] = "present"
    elif pos_cov >= 0.08:
        out["cls"] = "remnant"
    return out


def window_coverage(ref_fasta, genome_fasta, tmpdir, chunk=10_000):
    """cagPAI window +/-FLANK chopped into 10 kb chunks, each mapped
    independently with -cx asm20 (chunk as target, genome as query).
    Chunking prevents a high-scoring spurious gappy chain (de:f looks fine
    but nmatch/alen ~0.3-0.6, seen in empty-island genomes) from
    suppressing genuine small alignments inside the window, and the
    nmatch/alen >= 0.7 per-alignment filter rejects such forced
    alignments outright. Returns (win_cov, left_flank_cov, right_flank_cov,
    n_alignments, win_identity_weighted)."""
    ref = read_fasta(ref_fasta)
    chrom = list(ref.keys())[0]
    seq = ref[chrom]
    s0 = max(0, CAG_S - FLANK - 1)
    e0 = min(len(seq), CAG_E + FLANK)
    region = seq[s0:e0]
    win_a, win_b = CAG_S - 1 - s0, CAG_E - s0
    win_cov_b = union_cov_win = 0
    lf_cov_b = 0
    rf_cov_b = 0
    n_aln = 0
    idn_num = 0.0
    dirpath = os.path.join(tmpdir, "chunks")
    os.makedirs(dirpath, exist_ok=True)
    for c0 in range(0, len(region), chunk):
        c1 = min(len(region), c0 + chunk)
        cf = os.path.join(dirpath, f"chunk_{c0}.fna")
        with open(cf, "w") as fh:
            fh.write(f">c{c0}\n{region[c0:c1]}\n")
        paf = parse_paf(minimap_paf(cf, genome_fasta, "asm20"))
        good = [r for r in paf
                if r["tname"] == f"c{c0}" and
                r["nmatch"] / max(1, r["alen"]) >= 0.7]
        if not good:
            continue
        n_aln += len(good)
        iv = clip([(r["tstart"], r["tend"]) for r in good], 0, c1 - c0)
        cb = union_cov(iv)
        # shift into region coordinates and accumulate per zone
        ivw = clip([(s + c0, e + c0) for s, e in iv], win_a, win_b)
        win_cov_b += union_cov(ivw)
        lf_cov_b += union_cov(clip([(s + c0, e + c0) for s, e in iv], 0, win_a))
        rf_cov_b += union_cov(clip([(s + c0, e + c0) for s, e in iv],
                                   win_b, len(region)))
        for r in good:
            ov = min(r["tend"] + c0, win_b) - max(r["tstart"] + c0, win_a)
            if ov > 0:
                idn_num += ov * r["idn"]
    return (win_cov_b / WIN_LEN,
            lf_cov_b / max(1, win_a),
            rf_cov_b / max(1, len(region) - win_b),
            n_aln,
            idn_num / win_cov_b if win_cov_b else 0.0)


def bed_window_svs(bed_path):
    out = []
    ok = os.path.exists(bed_path)
    if ok:
        with open(bed_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                s, e = int(f[1]), int(f[2])
                if e < CAG_S or s > CAG_E:
                    continue
                out.append(f"{f[3].split('_')[0]}:{s}-{e}")
    return out, ok


# ------------------------------------------------------------------ verdict --
def classify(panel_status, panel_missing, caga, win_cov, lf, rf, n_pres):
    """panel_status: panel call; caga['cls']: k-mer based present/remnant/
    absent; win_cov: strict (real-identity) cagPAI window coverage."""
    notes = []
    flank_cov = max(lf, rf)
    locus_absent = win_cov < 0.05

    if caga["cls"] == "present":
        notes.append(
            f"cagA present in assembly (k-mer pos_cov={caga['kmer']}; asm20 "
            f"cov={caga['cov']}, idn={caga['idn']}, "
            f"{caga['contig']}:{caga['pos']}; strict asm5 cov="
            f"{caga['strict_cov'] or 'none'})")
        if panel_status in ("empty", "partial"):
            if win_cov >= 0.5:
                notes.append("cagPAI window largely covered: panel "
                             f"'{panel_status}' is not supported by the "
                             "assembly")
                return "panel_missed", "; ".join(notes)
            notes.append("but the cagPAI window itself is essentially "
                         "absent: the island backbone is deleted while "
                         "cagA is retained; panel 'empty' and collaborator "
                         "'Detected' are both defensible")
            return "disrupted_cagA", "; ".join(notes)
        return "other", "; ".join(notes) +             " (panel complete + cagA present; check boolean columns)"

    if caga["cls"] == "remnant":
        notes.append(
            f"cagA remnant only (k-mer pos_cov={caga['kmer']}; asm20 "
            f"cov={caga['cov']}, idn={caga['idn']}, "
            f"{caga['contig']}:{caga['pos']}; strict asm5 cov="
            f"{caga['strict_cov'] or 'none'}); window_cov={win_cov:.2f}")
        if panel_status == "empty" and locus_absent:
            notes.append("island backbone deleted, only a cagA fragment "
                         "left: panel 'empty' is assembly-supported; the "
                         "collaborator's call reflects the fragment")
        return "disrupted_cagA", "; ".join(notes)

    # cagA absent at k-mer level
    if locus_absent:
        if flank_cov >= 0.2:
            kind = ("clean deletion: flanks are covered and the assembly "
                    "supports a genuinely empty locus")
        else:
            kind = ("draft assembly gap: locus region and its flanks are "
                    "not assembled (fragmented draft), so the true locus "
                    "state cannot be read from this assembly")
        notes.append(f"cagPAI window uncovered ({kind})")
        if panel_status == "empty":
            notes.append("panel 'empty' is consistent with the assembly")
        return "assembly_gap", "; ".join(notes)
    if win_cov >= 0.5:
        if panel_status == "complete":
            notes.append(f"cagPAI window covered ({win_cov:.2f}, "
                         f"n_present={n_pres}) and island markers present, "
                         "but cagA has neither a minimap2 hit (asm20) nor "
                         "shared 19-mers: the island is present and cagA "
                         "is specifically deleted (or too diverged to map "
                         "at gene level)")
            return "true_biological_cagA_loss", "; ".join(notes)
        notes.append(f"window partly covered ({win_cov:.2f}) but cagA "
                     "absent: partial island without cagA")
        return "other", "; ".join(notes)
    notes.append(f"win_cov={win_cov:.2f}, flanks={lf:.2f}/{rf:.2f}, "
                 "cagA absent: mixed/ambiguous picture")
    return "other", "; ".join(notes)


# --------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-table",
                    default="/Users/macstudio/Downloads/Syn2bANI-paper/case_studies/"
                            "h_pylori_cagpai/results/metadata_assoc/merged_genome_table.tsv")
    ap.add_argument("--genomes-dir",
                    default="/Volumes/MoneyCat/Data/song_2026_hpylori/genomes")
    ap.add_argument("--markers",
                    default="/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status/cagpai_markers.fna")
    ap.add_argument("--ref26695",
                    default="/Users/macstudio/Downloads/Syn2bANI-paper/data/cagpai_pilot/hp26695.fna")
    ap.add_argument("--bed-dir",
                    default="/Volumes/MoneyCat/Data/song_2026_hpylori/struct_vs_26695_filtered")
    ap.add_argument("--outdir",
                    default="/Users/macstudio/Downloads/Syn2bANI-paper/case_studies/"
                            "h_pylori_cagpai/results/metadata_assoc/discordant_review")
    ap.add_argument("--minimap2", default=MINIMAP2)
    args = ap.parse_args()
    globals()["MINIMAP2"] = args.minimap2

    os.makedirs(args.outdir, exist_ok=True)

    # ---- cohort stats for size context ------------------------------------
    with open(args.merged_table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # ---- 1. discordant set --------------------------------------------------
    discordant = [r for r in rows
                  if r["cagA_detected"] == "" or
                  r["cagpai_complete"] != r["cagA_detected"]]
    print(f"[1] {len(rows)} genomes in table, {len(discordant)} discordant",
          file=sys.stderr)

    keep = ["GCA_accession", "status", "status_extended", "cagA_result",
            "cagA_genotype", "fastbaps", "country", "group", "n_cag_sv",
            "cag_sv_types", "cag_sv_list", "missing_markers", "n_present",
            "cagpai_complete", "cagA_detected"]

    out_cols = keep + [
        "n_contigs", "total_len", "N50",
        "n_present_panel_exact", "panel_exact_disagree",
        "cagA_class", "cagA_kmer_pos_cov", "cagA_idn", "cagA_cov", "cagA_contig", "cagA_pos",
        "cagA_strict_asm5_idn", "cagA_strict_asm5_cov",
        "window_cov", "window_idn", "left_flank_cov", "right_flank_cov",
        "window_map_blocks", "bed_sv_in_window", "sv_list_matches_bed",
        "verdict", "notes",
    ]

    results = []
    with tempfile.TemporaryDirectory(dir=args.outdir) as tmpdir:
        markers = read_fasta(args.markers)
        if CAGA not in markers:
            raise SystemExit(f"cagA marker {CAGA!r} not found in {args.markers}")
        caga_fa = os.path.join(tmpdir, "cagA_only.fna")
        with open(caga_fa, "w") as fh:
            fh.write(f">{CAGA}\n{markers[CAGA]}\n")

        for i, r in enumerate(discordant, 1):
            gca = r["GCA_accession"]
            gfa = os.path.join(args.genomes_dir, f"{gca}.fna")
            row = {k: r.get(k, "") for k in keep}
            row.update(dict.fromkeys(out_cols[len(keep):], ""))
            if not os.path.exists(gfa):
                row["verdict"], row["notes"] = "other", "FASTA missing on drive"
                results.append(row)
                continue
            print(f"  [{i}/{len(discordant)}] {gca} ...", file=sys.stderr)

            seqs = read_fasta(gfa)
            n_ctg, tot, n50 = assembly_stats(seqs)

            n_pres, mstats = marker_recheck_panel_exact(args.markers, gfa)
            miss_tbl = set(m.strip() for m in r["missing_markers"].split(",") if m.strip())
            disagree = []
            for m in MARKERS:
                pe = "present" if mstats[m]["present"] else "absent"
                tbl = "absent" if m in miss_tbl else "present"
                if pe != tbl:
                    s = mstats[m]
                    disagree.append(f"{m}(tbl={tbl},exact={pe},"
                                    f"cov={s['cov']:.2f},idn={s['idn']:.2f})"

                                    )
            caga = caga_check(caga_fa, gfa, markers[CAGA])
            win_cov, lf, rf, nblocks, win_idn = window_coverage(
                args.ref26695, gfa, tmpdir)
            bed_svs, bed_ok = bed_window_svs(
                os.path.join(args.bed_dir, f"{gca}.vs_hp26695.bed"))
            tbl_svs = ([] if r["cag_sv_list"] in ("", "none")
                       else r["cag_sv_list"].split(";"))
            sv_match = "yes" if set(bed_svs) == set(tbl_svs) else \
                       ("bed_missing" if not bed_ok else "no")

            verdict, notes = classify(r["status"], miss_tbl, caga, win_cov,
                                      lf, rf, n_pres)

            row.update({
                "n_contigs": n_ctg, "total_len": tot, "N50": n50,
                "n_present_panel_exact": n_pres,
                "panel_exact_disagree": ";".join(disagree) or "none",
                "cagA_class": caga["cls"], "cagA_kmer_pos_cov": caga["kmer"], "cagA_idn": caga["idn"],
                "cagA_cov": caga["cov"], "cagA_contig": caga["contig"],
                "cagA_pos": caga["pos"],
                "cagA_strict_asm5_idn": caga["strict_idn"],
                "cagA_strict_asm5_cov": caga["strict_cov"],
                "window_cov": f"{win_cov:.4f}",
                "window_idn": f"{win_idn:.3f}",
                "left_flank_cov": f"{lf:.4f}",
                "right_flank_cov": f"{rf:.4f}",
                "window_map_blocks": nblocks,
                "bed_sv_in_window": ";".join(bed_svs) or "none",
                "sv_list_matches_bed": sv_match,
                "verdict": verdict, "notes": notes,
            })
            results.append(row)

    out_tsv = os.path.join(args.outdir, "discordant_genomes_review.tsv")
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"[2] wrote {out_tsv}", file=sys.stderr)

    from collections import Counter
    print("\n== verdict counts ==")
    for v, n in Counter(r["verdict"] for r in results).most_common():
        print(f"  {v}\t{n}")


if __name__ == "__main__":
    main()
