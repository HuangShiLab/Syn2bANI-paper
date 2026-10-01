#!/usr/bin/env python3
"""Normalize legacy GTDB census rows and prepare globally deduplicated stream.

Early GTDB tasks used contig IDs and retained within-block rows in cross-block
outputs. This program maps contig IDs back to genome accessions, canonicalizes
unordered pairs, and emits rows for deterministic external-sort deduplication.
"""
import argparse, collections, json, sys
from pathlib import Path


def load_contig_map(root, included_by_cluster):
    m = {}
    conflicts = collections.defaultdict(set)
    for cluster, accs in included_by_cluster.items():
        for acc in accs:
            m[acc] = acc
            tgt = root / "tgt" / f"{acc}.tgt"
            if not tgt.is_file():
                continue
            with tgt.open(errors="replace") as fh:
                line = fh.readline().strip()
            if not line.startswith(">"):
                continue
            token = line[1:].split()[0]
            candidates = {token}
            if "|" in token:
                candidates.update(token.split("|"))
            for key in candidates:
                if key in m and m[key] != acc:
                    conflicts[key].update([m[key], acc])
                else:
                    m[key] = acc
    return m, conflicts


def resolve(x, m, unmapped):
    if x in m:
        return m[x]
    base = x.split("|", 1)[0]
    if base in m:
        return m[base]
    unmapped[x] += 1
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--stats", required=True)
    args = ap.parse_args()
    root = Path(args.workdir)
    included = {}
    for p in (root / "effective_cluster_accessions").glob("*.txt"):
        included[p.name[:-4]] = p.read_text().split()
    mapping, conflicts = load_contig_map(root, included)
    unmapped = collections.Counter()
    total = emitted = skipped_self = skipped_unmapped = 0
    ia = ib = None
    for p in sorted((root / "outputs").glob("*.tsv.gz")):
        import gzip
        with gzip.open(p, "rt") as fh:
            for li, line in enumerate(fh):
                if li == 0:
                    if ia is None:
                        fields = line.rstrip("\n").split("\t")
                        ia, ib = fields.index("genome_A"), fields.index("genome_B")
                        print("\t".join(fields), flush=True)
                    continue
                total += 1
                f = line.rstrip("\n").split("\t")
                a = resolve(f[ia], mapping, unmapped)
                b = resolve(f[ib], mapping, unmapped)
                if a is None or b is None:
                    skipped_unmapped += 1
                    continue
                if a == b:
                    skipped_self += 1
                    continue
                if b < a:
                    a, b = b, a
                f[ia], f[ib] = a, b
                sys.stdout.write("\t".join(f) + "\n")
                emitted += 1
    with open(args.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in [("input_rows", total), ("emitted_rows", emitted),
                     ("skipped_self", skipped_self),
                     ("skipped_unmapped", skipped_unmapped),
                     ("mapping_keys", len(mapping)),
                     ("mapping_conflict_keys", len(conflicts))]:
            fh.write(f"{k}\t{v}\n")
        if unmapped:
            fh.write("unmapped_examples\t" + ",".join(list(unmapped)[:50]) + "\n")
        if conflicts:
            fh.write("conflict_examples\t" + ",".join(list(conflicts)[:50]) + "\n")
    print(f"normalization: input={total} emitted={emitted} unmapped={skipped_unmapped}",
          file=sys.stderr)
    if skipped_unmapped or conflicts:
        raise SystemExit("normalization found unmapped/conflicting genome IDs")


if __name__ == "__main__":
    main()
