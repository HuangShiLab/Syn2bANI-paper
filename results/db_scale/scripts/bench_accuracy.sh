#!/bin/bash
# Accuracy cross-check: ~500 pairs stratified by ANI band, drawn from the
# validated 100k matrix (results/matrix_gtdb_r207_100k_v8_final.tsv, s2b_ani = ani path).
# Per pair: syn2bani ani (MLE, ground truth), dist --enzymes ENZ4 (chained-kmer output),
# dist --enzymes ENZ4 --mash-ani (GBRT v7 output).
source "$(dirname "$0")/common.sh"
NPAIRS=${1:-500}
TSV="$BASE/dist_vs_ani.tsv"
PAIRS="$BASE/lists/accuracy_pairs.tsv"
MATRIX=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/matrix_gtdb_r207_100k_v8_final.tsv

if [ ! -s "$PAIRS" ]; then
    python3 - "$MATRIX" "$GENOMES" "$NPAIRS" > "$PAIRS" <<'EOF'
import csv, os, random, sys
matrix, gdir, npairs = sys.argv[1], sys.argv[2], int(sys.argv[3])
random.seed(42)
bands = {"80-85": [], "85-90": [], "90-95": [], "95-100": []}
with open(matrix) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        try:
            ani = float(row["s2b_ani"])
        except (ValueError, KeyError):
            continue
        if ani < 80: continue
        q = os.path.join(gdir, row["query"] + ".fna")
        r = os.path.join(gdir, row["reference"] + ".fna")
        if not (os.path.exists(q) and os.path.exists(r)): continue
        b = "80-85" if ani < 85 else "85-90" if ani < 90 else "90-95" if ani < 95 else "95-100"
        bands[b].append((q, r, ani))
per = npairs // 4
out = []
for b, rows in bands.items():
    random.shuffle(rows)
    take = rows[:per]
    out.extend(take)
    print(f"# band {b}: {len(rows)} available, took {len(take)}", file=sys.stderr)
for q, r, ani in out:
    print(f"{q}\t{r}\t{ani:.4f}")
EOF
fi

printf 'query\tref\ttruth_ani_mle_v8\tani_mle_rerun\tani_gated\tdist_chainedkmer\tdist_gbrt_v7\tdist_af_q\tdist_af_r\tdist_shared_tags\tdist_status\n' > "$TSV"

run_one() {
    q=$1; r=$2; truth=$3
    ani_out=$("$S2B" ani "$q" "$r" 2>/dev/null | tail -1)
    ani_mle=$(echo "$ani_out" | cut -f3); ani_gated=$(echo "$ani_out" | awk -F'\t' '{print $(NF-1)}')
    d_out=$("$S2B" dist "$q" "$r" --enzymes "$ENZ4" -p -t 1 2>/dev/null | tail -1)
    if [ -n "$d_out" ] && [ "${d_out%%$'\t'*}" != "query_file" ]; then
        d_ck=$(echo "$d_out" | cut -f3); d_afq=$(echo "$d_out" | cut -f4); d_afr=$(echo "$d_out" | cut -f5); d_st=$(echo "$d_out" | cut -f8)
        status=ok
    else
        d_ck=NA; d_afq=NA; d_afr=NA; d_st=NA; status=screened_out
    fi
    d_g=$("$S2B" dist "$q" "$r" --enzymes "$ENZ4" --mash-ani -p -t 1 2>/dev/null | tail -1)
    if [ -z "$d_g" ] || [ "${d_g%%$'\t'*}" = "query_file" ]; then d_g=NA; else d_g=$(echo "$d_g" | cut -f3); fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(basename "$q" .fna)" "$(basename "$r" .fna)" \
        "$truth" "$ani_mle" "$ani_gated" "$d_ck" "$d_g" "$d_afq" "$d_afr" "$d_st" "$status" >> "$TSV"
}
export -f run_one
export TSV S2B ENZ4

echo "[$(date -u +%H:%M:%S)] accuracy cross-check: $(wc -l < "$PAIRS") pairs"
cat "$PAIRS" | xargs -P "$THREADS" -L 1 bash -c 'run_one "$@"' _
echo "[$(date -u +%H:%M:%S)] accuracy cross-check done: $(($(wc -l < "$TSV") - 1)) rows"
