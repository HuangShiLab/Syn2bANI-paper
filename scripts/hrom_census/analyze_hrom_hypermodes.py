#!/usr/bin/env python3
"""SynTracker-style top-5% species-mode analysis for the HROM census.

The HROM census has two complementary similarity axes:

- skani ANI: high values indicate sequence-similar pairs;
- Syn2b ``structural`` score: high values indicate conserved landmark
  adjacencies (an APSS-like structural-similarity score).

For each axis, this script selects the top ``--top-fraction`` pairs. Species
enrichment in each selected set is tested with a one-sided hypergeometric test
and controlled by Benjamini-Hochberg FDR.  The script writes compact plotting
data; use ``plot_hrom_hypermodes.py`` for the final two-panel figure.
"""
from __future__ import annotations

import argparse
import array
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom, spearmanr


def bh_qvalues(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.minimum(ranked, 1.0)
    q = np.empty(n, dtype=float)
    q[order] = ranked
    return q


def clean_species_name(raw: str) -> str:
    name = raw.strip()
    if name.startswith("s__"):
        name = name[3:]
    # Representative species fields often append the cluster ID.
    if "/HROM_Genome_" in name:
        name = name.split("/HROM_Genome_", 1)[0]
    return name or "HROM species"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--census", required=True,
                    help="hrom_within_species_sv_ani.tsv.gz")
    ap.add_argument("--species-metadata", required=True,
                    help="HROM_representative_genome_metadata.tsv")
    ap.add_argument("--genome-cluster-metadata", required=True,
                    help="hrom_genome_metadata.tsv mapping every genome to "
                         "its HROM representative cluster")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-fraction", type=float, default=0.05)
    ap.add_argument("--min-shared-tags", type=int, default=250,
                    help="minimum shared landmarks (default: 250)")
    ap.add_argument("--min-observable-fraction", type=float, default=0.50,
                    help="minimum observable_fraction (default: 0.50)")
    ap.add_argument("--min-species-pairs", type=int, default=200,
                    help="minimum total species pairs for enrichment testing")
    ap.add_argument("--plot-background-points", type=int, default=100000)
    ap.add_argument("--plot-selected-points-per-group", type=int, default=50000)
    ap.add_argument("--correlation-sample", type=int, default=500000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if not 0 < args.top_fraction < 1:
        raise ValueError("--top-fraction must be between 0 and 1")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Representative metadata maps HROM_Genome_XXXX to a display species name.
    cluster_species = {}
    with open(args.species_metadata, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cluster = row["genome"]
            raw = row.get("species") or row.get("gtdb_taxonomy") or cluster
            if raw.startswith("s__"):
                raw = raw[3:]
            name = clean_species_name(raw)
            cluster_species[cluster] = name
    genome_cluster = {}
    with open(args.genome_cluster_metadata, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genome_cluster[row["genome"]] = row["cluster"]
    # Disambiguate repeated names while retaining compact plot labels.
    counts = Counter(cluster_species.values())
    species_labels = {
        cl: (name if counts[name] == 1 else f"{name} [{cl}]")
        for cl, name in cluster_species.items()
    }
    species_names = sorted(species_labels)
    species_code = {name: i for i, name in enumerate(species_names)}

    ani_arr = array.array("d")
    struct_arr = array.array("d")
    code_arr = array.array("I")
    n_seen = 0
    n_dropped_quality = 0

    with gzip.open(args.census, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ["genome_A", "skani_ani", "structural", "shared_tags",
                    "observable_fraction"]
        missing = [x for x in required if x not in idx]
        if missing:
            raise ValueError(f"census lacks columns: {missing}")
        i_ani, i_struct = idx["skani_ani"], idx["structural"]
        i_shared, i_obs = idx["shared_tags"], idx["observable_fraction"]
        i_a, i_b = idx["genome_A"], idx["genome_B"]

        for line in fh:
            n_seen += 1
            f = line.rstrip("\n").split("\t")
            try:
                shared = int(f[i_shared])
                observable = float(f[i_obs])
                ani = float(f[i_ani])
                structural = float(f[i_struct])
            except ValueError:
                n_dropped_quality += 1
                continue
            if shared < args.min_shared_tags:
                n_dropped_quality += 1
                continue
            if observable < args.min_observable_fraction:
                n_dropped_quality += 1
                continue
            cluster = genome_cluster.get(f[i_a])
            # Census rows should be within-cluster; verify rather than assume.
            if cluster is None or genome_cluster.get(f[i_b]) != cluster:
                n_dropped_quality += 1
                continue
            code = species_code.get(cluster)
            if code is None:
                n_dropped_quality += 1
                continue
            ani_arr.append(ani)
            struct_arr.append(structural)
            code_arr.append(code)

    ani = np.frombuffer(ani_arr, dtype=np.float64)
    structural = np.frombuffer(struct_arr, dtype=np.float64)
    codes = np.frombuffer(code_arr, dtype=np.uint32)
    n_kept = len(ani)
    if n_kept == 0:
        raise ValueError("no pairs passed quality filters")

    n_top = int(math.floor(args.top_fraction * n_kept))
    if n_top == 0:
        raise ValueError("top-fraction selects no pairs")

    rng = np.random.default_rng(args.seed)
    ani_order = np.argpartition(ani, n_kept - n_top)[n_kept - n_top:]
    ani_top = np.zeros(n_kept, dtype=bool)
    ani_top[ani_order] = True
    struct_order = np.argpartition(structural, n_kept - n_top)[n_kept - n_top:]
    struct_top = np.zeros(n_kept, dtype=bool)
    struct_top[struct_order] = True
    both_top = ani_top & struct_top

    total_by_species = np.bincount(codes, minlength=len(species_names))
    ani_by_species = np.bincount(codes[ani_top], minlength=len(species_names))
    struct_by_species = np.bincount(codes[struct_top],
                                    minlength=len(species_names))
    both_by_species = np.bincount(codes[both_top],
                                  minlength=len(species_names))

    # Correlation is computed on a fixed random sample to avoid repeated
    # sorting of all 62.9M values.
    corr_n = min(args.correlation_sample, n_kept)
    corr_idx = rng.choice(n_kept, size=corr_n, replace=False)
    rho, corr_p = spearmanr(ani[corr_idx], structural[corr_idx])

    rows = []
    eligible_codes = []
    for code, name in enumerate(species_names):
        K = int(total_by_species[code])
        if K < args.min_species_pairs:
            continue
        k_struct = int(struct_by_species[code])
        k_ani = int(ani_by_species[code])
        k_both = int(both_by_species[code])
        # logsf retains significance when the survival probability underflows
        # to zero for extreme HROM clusters.
        log_p_struct = hypergeom.logsf(k_struct - 1, n_kept, K, n_top)
        log_p_ani = hypergeom.logsf(k_ani - 1, n_kept, K, n_top)
        p_struct = math.exp(log_p_struct)
        p_ani = math.exp(log_p_ani)
        r_struct = (k_struct / n_top) / (K / n_kept) if K and n_top else math.nan
        r_ani = (k_ani / n_top) / (K / n_kept) if K and n_top else math.nan
        pseudocount = 0.5
        log2_ratio = math.log2(
            (k_struct + pseudocount) / (k_ani + pseudocount))
        rows.append({
            "cluster": species_names[code],
            "species": species_labels[species_names[code]],
            "total_pairs": K,
            "struct_top_pairs": k_struct,
            "ani_top_pairs": k_ani,
            "both_top_pairs": k_both,
            "struct_enrichment_ratio": r_struct,
            "ani_enrichment_ratio": r_ani,
            "struct_p": p_struct,
            "ani_p": p_ani,
            "log2_struct_over_ani": log2_ratio,
            "neg_log10_min_p": -min(log_p_struct, log_p_ani),
        })
        eligible_codes.append(code)

    q_struct = bh_qvalues([r["struct_p"] for r in rows])
    q_ani = bh_qvalues([r["ani_p"] for r in rows])
    for i, r in enumerate(rows):
        r["struct_q"] = float(q_struct[i])
        r["ani_q"] = float(q_ani[i])
        if r["struct_q"] < 0.05 and r["struct_q"] <= r["ani_q"] \
                and r["log2_struct_over_ani"] <= -1:
            mode = "structural_enriched"
        elif r["ani_q"] < 0.05 and r["ani_q"] <= r["struct_q"] \
                and r["log2_struct_over_ani"] >= 1:
            mode = "ani_enriched"
        elif r["struct_q"] < 0.05 and r["ani_q"] < 0.05:
            mode = "both_enriched"
        else:
            mode = "not_enriched"
        r["mode"] = mode
    rows.sort(key=lambda r: (-r["neg_log10_min_p"], r["species"]))

    enrichment_tsv = out / "hypermode_species_enrichment.tsv"
    fields = list(rows[0].keys()) if rows else []
    with enrichment_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # Compact plotting data. Selected points are sampled to keep PDFs small;
    # all enrichment statistics above use every qualifying pair.
    background_pool = np.flatnonzero(~(ani_top | struct_top))
    background_idx = rng.choice(background_pool,
                                size=min(args.plot_background_points,
                                         len(background_pool)), replace=False)
    group_indices = {}
    for name, mask in [("ani_top", ani_top), ("structural_top", struct_top),
                       ("both_top", both_top)]:
        pool = np.flatnonzero(mask)
        group_indices[name] = rng.choice(
            pool,
            size=min(args.plot_selected_points_per_group, len(pool)),
            replace=False)
    # plot_data_idx contains indices in the filtered arrays and in the full
    # census arrays; positions are needed for the sampled group labels.
    plot_data_idx = np.concatenate(
        [background_idx] + list(group_indices.values()))
    plot_groups = np.full(len(plot_data_idx), "background", dtype=object)
    for name, idx in group_indices.items():
        plot_groups[np.isin(plot_data_idx, idx)] = name
    plot_tsv = out / "hypermode_figure_points.tsv.gz"
    with gzip.open(plot_tsv, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["ani", "structural", "group"])
        for pos, data_idx in enumerate(plot_data_idx):
            writer.writerow([f"{ani[data_idx]:.6f}",
                             f"{structural[data_idx]:.6f}",
                             plot_groups[pos]])

    # Exact counts of categorical groups (not the plotted samples).
    counts = {
        "kept_pairs": int(n_kept),
        "ani_top_pairs": int(ani_top.sum()),
        "structural_top_pairs": int(struct_top.sum()),
        "both_top_pairs": int(both_top.sum()),
    }
    report = {
        "parameters": vars(args) | {"species_metadata": args.species_metadata},
        "pairs_seen": n_seen,
        "pairs_dropped_quality_or_unknown_cluster": n_dropped_quality,
        "counts": counts,
        "selection_thresholds": {
            "ani_min_in_top": float(ani[ani_top].min()),
            "structural_min_in_top": float(structural[struct_top].min()),
        },
        "spearman_ani_vs_structural": {
            "n": int(corr_n), "rho": float(rho), "p": float(corr_p),
        },
        "species_tested": len(rows),
        "modes": dict(Counter(r["mode"] for r in rows)),
        "plot_data": str(plot_tsv),
        "enrichment_table": str(enrichment_tsv),
    }
    with (out / "hypermode_analysis_report.json").open("w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
