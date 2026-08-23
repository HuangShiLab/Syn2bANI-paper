#!/usr/bin/env python3
"""Sample candidate high-ANI pairs from the full GTDB-R207 genome collection.

Produces two strata for the unified benchmark:
  95-97% : same-genus, different-species pairs (inter-species)
  97-100%: same-species pairs (intra-species)

Excludes:
  - genomes used in the v5 calibration training set (2,520 pairs)
  - genomes used in the 50k held-out set (43,334 pairs)
  - GTDB representative genomes (to keep the stratum orthogonal to rep-based sets)

Outputs (results/gtdb50k/):
  high_ani_candidates.tsv  pairid q_acc r_acc stratum genus species info
  high_ani_genomes.txt     unique accessions to download
"""
import os
import random
import itertools
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")
METADIR = "/lustre1/g/aos_shihuang/data/gtdb-r207/metadata"

TARGET_PER_STRATUM = 4000  # candidate pairs; will be train/test split later
RANDOM_SEED = 20240824


def strip_prefix(acc):
    if acc.startswith("GB_") or acc.startswith("RS_"):
        return acc[3:]
    return acc


def parse_tax(tax):
    parts = {k[0]: k for item in tax.split(";") for k in [item]}
    return parts.get("g", ""), parts.get("s", "")


def load_excluded():
    excl = set()
    # training set: eval_pairs 2,074 + hi95 467
    train1 = pd.read_csv(os.path.join(ROOT, "results", "panel_by_band", "eval_pairs.tsv"), sep="\t")
    train2 = pd.read_csv(os.path.join(ROOT, "results", "anim_truth_hi95.tsv"), sep="\t")
    for col in ["query_asm", "ref_asm", "query", "reference"]:
        if col in train1.columns:
            excl.update(train1[col].astype(str))
        if col in train2.columns:
            excl.update(train2[col].astype(str))
    # 50k held-out genomes
    pairs50 = pd.read_csv(os.path.join(RES, "pairs_50k.tsv"), sep="\t")
    excl.update(pairs50["q_acc"].astype(str))
    excl.update(pairs50["r_acc"].astype(str))
    # normalize: strip prefix and version? keep both
    return excl | {strip_prefix(a) for a in excl}


def sample_pairs(items, max_per_group):
    """Sample up to max_per_group unordered pairs from items without replacement."""
    if len(items) < 2:
        return []
    combos = list(itertools.combinations(sorted(items), 2))
    if len(combos) <= max_per_group:
        return combos
    return random.sample(combos, max_per_group)


def main():
    random.seed(RANDOM_SEED)
    excluded = load_excluded()

    dfs = []
    for f in ["bac120_metadata_r207.tsv", "ar53_metadata_r207.tsv"]:
        df = pd.read_csv(os.path.join(METADIR, f), sep="\t", usecols=["accession", "gtdb_representative", "gtdb_taxonomy"])
        dfs.append(df)
    meta = pd.concat(dfs, ignore_index=True)
    meta["acc"] = meta["accession"].apply(strip_prefix)
    meta["is_rep"] = meta["gtdb_representative"].astype(str).str.lower().isin(["t", "true"])
    meta["genus"], meta["species"] = zip(*meta["gtdb_taxonomy"].apply(parse_tax))

    # keep only non-excluded, non-rep genomes with species-level taxonomy
    keep = meta[(~meta["is_rep"]) & (~meta["acc"].isin(excluded)) & (meta["species"] != "")].copy()
    print(f"eligible non-rep genomes: {len(keep):,}")

    # 97-100: intra-species pairs
    sp_groups = keep.groupby("species")["acc"].apply(list).to_dict()
    sp_groups = {k: v for k, v in sp_groups.items() if len(v) >= 2}
    print(f"species with >=2 eligible genomes: {len(sp_groups):,}")
    # sample species first to spread across species, then pairs within each
    chosen_species = random.sample(sorted(sp_groups.keys()), min(len(sp_groups), 2000))
    intra = []
    for sp in chosen_species:
        for a, b in sample_pairs(sp_groups[sp], 3):
            intra.append((a, b, "97-100", keep[keep.acc == a]["genus"].iloc[0], sp))
            if len(intra) >= TARGET_PER_STRATUM:
                break
        if len(intra) >= TARGET_PER_STRATUM:
            break
    print(f"intra-species candidates: {len(intra):,}")

    # 95-97: inter-species same-genus pairs
    # build genus -> species -> list of genomes
    genus_df = keep[keep["genus"] != ""].copy()
    genus_groups = {}
    for g, sub in genus_df.groupby("genus"):
        sp_acc = {}
        for sp, accs in sub.groupby("species")["acc"].apply(list).items():
            sp_acc[sp] = accs
        if len(sp_acc) >= 2:
            genus_groups[g] = sp_acc
    print(f"genera with >=2 species eligible: {len(genus_groups):,}")

    chosen_genera = random.sample(sorted(genus_groups.keys()), min(len(genus_groups), 2000))
    inter = []
    for g in chosen_genera:
        sp_acc = genus_groups[g]
        # sample up to 2 genomes per species, then cross-species pairs
        items = []
        for sp, accs in sp_acc.items():
            items.extend(random.sample(accs, min(2, len(accs))))
        for a, b in sample_pairs(items, 3):
            inter.append((a, b, "95-97", g, ""))
            if len(inter) >= TARGET_PER_STRATUM:
                break
        if len(inter) >= TARGET_PER_STRATUM:
            break
    print(f"inter-species same-genus candidates: {len(inter):,}")

    out = pd.DataFrame(intra + inter, columns=["q_acc", "r_acc", "stratum", "genus", "species"])
    out["pairid"] = out["q_acc"] + "__" + out["r_acc"]
    out = out[["pairid", "q_acc", "r_acc", "stratum", "genus", "species"]]
    os.makedirs(RES, exist_ok=True)
    out.to_csv(os.path.join(RES, "high_ani_candidates.tsv"), sep="\t", index=False)

    genomes = sorted(set(out["q_acc"]) | set(out["r_acc"]))
    with open(os.path.join(RES, "high_ani_genomes.txt"), "w") as fh:
        for g in genomes:
            fh.write(g + "\n")
    print(f"unique genomes to download: {len(genomes):,}")
    print("wrote high_ani_candidates.tsv and high_ani_genomes.txt")


if __name__ == "__main__":
    main()
