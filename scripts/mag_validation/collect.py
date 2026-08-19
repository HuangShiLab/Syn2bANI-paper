#!/usr/bin/env python3
"""collect.py <work>
Merge per-pair tool outputs into collect/:
  bins.tsv           (copied through from cohort stage)
  truth_dnadiff.tsv  (bin, anim_ani, anim_af_ref, anim_af_mag, aligned_bases_mag, af_tier)
  ani_fast_tools.tsv (bin, ref_path, ref_id, role, syn2bani*, skani*, fastani*)
syn2bani labels each side by the fasta's FIRST RECORD ID, so references are
resolved via the firstid maps built in s0/s0b, and bins via '{binid}|ctg'.
"""
import glob
import os
import sys

WORK = sys.argv[1]
PP = f"{WORK}/fast/per_pair"
RLD = f"{WORK}/fast/rl"

def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            vals = line.rstrip("\n").split("\t")
            yield header, vals

# first-record-id -> ref path (source genomes + GTDB)
fid2ref = {}
for m in glob.glob(f"{WORK}/refs/*/firstid_map.tsv"):
    with open(m) as f:
        next(f)
        base = os.path.dirname(m)
        paths = {}
        gp = os.path.join(base, "genome_paths.txt")
        genomes_tsv = os.path.join(base, "genomes.tsv")
        if os.path.exists(genomes_tsv):
            with open(genomes_tsv) as g:
                for line in g:
                    if line.startswith("genome_id\t"):
                        continue
                    gid, p = line.rstrip("\n").split("\t")
                    paths[gid] = p
        elif os.path.exists(gp):
            with open(gp) as g:
                for line in g:
                    p = line.strip()
                    gid = p.rsplit("/", 1)[-1].replace(".fna", "")
                    paths[gid] = p
        for line in f:
            gid, fid = line.rstrip("\n").split("\t")
            if gid in paths:
                fid2ref[fid] = paths[gid]

# ref path -> ref id (basename w/o extension, for readability)
def ref_id(path):
    b = path.rsplit("/", 1)[-1]
    for ext in (".fasta", ".fna", ".fa"):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b

# roles per bin
roles = {}
for p in glob.glob(f"{RLD}/*.refs.tsv"):
    with open(p) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                binid = os.path.basename(p)[:-len(".refs.tsv")]
                roles[(binid, parts[0])] = parts[1]

# --- gather per-bin tool outputs keyed by (bin, ref_path) ---
pairs = {}

def get(binid, refpath):
    return pairs.setdefault((binid, refpath), {})

for p in glob.glob(f"{PP}/*.syn2bani.tsv"):
    binid = os.path.basename(p)[: -len(".syn2bani.tsv")]
    try:
        with open(p) as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) < len(header):
                    continue
                q = v[idx["query"]].split("|")[0]
                refpath = fid2ref.get(v[idx["reference"]], v[idx["reference"]])
                d = get(q, refpath)
                d["s2b_ani"] = v[idx.get("ani", -1)]
                d["s2b_ani_gated"] = v[idx.get("ani_gated", -1)] if "ani_gated" in idx else ""
                d["s2b_ani_upper95"] = v[idx.get("ani_upper95", -1)] if "ani_upper95" in idx else ""
                d["s2b_gate"] = v[idx.get("gate", -1)] if "gate" in idx else ""
                d["s2b_af_query"] = v[idx.get("af_query", -1)] if "af_query" in idx else ""
                d["s2b_af_ref"] = v[idx.get("af_reference", -1)] if "af_reference" in idx else ""
                d["s2b_n_anchors"] = v[idx.get("n_anchors", -1)] if "n_anchors" in idx else ""
                d["s2b_n_chains"] = v[idx.get("n_chains", -1)] if "n_chains" in idx else ""
                d["s2b_std_err"] = v[idx.get("std_err", -1)] if "std_err" in idx else ""
                d["s2b_flag"] = v[-1]  # trailing flag column (header lists 'flag' twice)
    except Exception as e:
        print(f"[collect] WARN syn2bani parse {p}: {e}", file=sys.stderr)

for p in glob.glob(f"{PP}/*.skani.tsv"):
    binid = os.path.basename(p)[: -len(".skani.tsv")]
    try:
        with open(p) as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) < 5:
                    continue
                q = v[idx["Query_file"]]
                q = q.rsplit("/", 1)[-1]
                q = q[:-3] if q.endswith(".fa") else q
                refpath = v[idx["Ref_file"]]
                d = get(q or binid, refpath)
                d["skani_ani"] = v[idx.get("ANI", 2)]
                d["skani_af_ref"] = v[idx.get("Align_fraction_ref", -1)] if "Align_fraction_ref" in idx else ""
                d["skani_af_query"] = v[idx.get("Align_fraction_query", -1)] if "Align_fraction_query" in idx else ""
    except Exception as e:
        print(f"[collect] WARN skani parse {p}: {e}", file=sys.stderr)

for p in glob.glob(f"{PP}/*.fastani.tsv"):
    binid = os.path.basename(p)[: -len(".fastani.tsv")]
    try:
        with open(p) as f:
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) < 5:
                    continue
                q = v[0].rsplit("/", 1)[-1]
                q = q[:-3] if q.endswith(".fa") else q
                refpath = v[1]
                d = get(q, refpath)
                d["fastani_ani"] = v[2]
                d["fastani_ortho"] = v[3]
                d["fastani_frags"] = v[4]
    except Exception as e:
        print(f"[collect] WARN fastani parse {p}: {e}", file=sys.stderr)

os.makedirs(f"{WORK}/collect", exist_ok=True)

# truth table
with open(f"{WORK}/collect/truth_dnadiff.tsv", "w") as o:
    o.write("bin\tanim_ani\tanim_af_ref\tanim_af_mag\taligned_bases_mag\taf_tier\n")
    for p in sorted(glob.glob(f"{WORK}/truth/rows/chunk_*.tsv")):
        with open(p) as f:
            for line in f:
                o.write(line)

# fast-tools table
cols = ["bin", "ref_path", "ref_id", "role",
        "s2b_ani", "s2b_ani_gated", "s2b_ani_upper95", "s2b_gate", "s2b_flag",
        "s2b_af_query", "s2b_af_ref", "s2b_n_anchors", "s2b_n_chains", "s2b_std_err",
        "skani_ani", "skani_af_ref", "skani_af_query",
        "fastani_ani", "fastani_ortho", "fastani_frags"]
with open(f"{WORK}/collect/ani_fast_tools.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for (binid, refpath), d in sorted(pairs.items()):
        row = [binid, refpath, ref_id(refpath), roles.get((binid, refpath), "?")]
        row += [d.get(c, "") for c in cols[4:]]
        o.write("\t".join(row) + "\n")

# rep map passthrough
with open(f"{WORK}/collect/rep_map.tsv", "w") as o:
    o.write("bin\trep_path\tskani_ani\tskani_af_query\n")
    for p in sorted(glob.glob(f"{WORK}/repsearch/*/map.tsv")):
        with open(p) as f:
            for line in f:
                if not line.startswith("bin\t"):
                    o.write(line)

print(f"[collect] {len(pairs)} (bin,ref) pairs in ani_fast_tools.tsv")
