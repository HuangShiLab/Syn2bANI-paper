#!/bin/bash
#SBATCH --job-name=closed_top100_junctions
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/closed_inversions/logs/top100_junctions_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/closed_inversions/logs/top100_junctions_%j.err
# Run dnadiff + junction-coordinate comparison on the top 100 inverted closed-genome pairs.
set -euo pipefail

ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=/lustre1/g/aos_shihuang/data/closed_inversions
GENOMES=$WORK/genomes
OUT=$WORK/top100_dnadiff
TGT_CACHE=$WORK/syn2b_tgts_cache
PY=/group/aos_shihuang/conda/bin/python3
ANVIO_BIN=/group/aos_shihuang/conda/envs/anvio/bin
export PATH="$ANVIO_BIN:$PATH"
export PYTHONNOUSERSITE=1

mkdir -p "$OUT" "$WORK/logs" "$WORK/top100_results"

PAIRS=$WORK/top100_inverted_pairs.tsv

# 1. Run dnadiff for each pair, keeping dd.report and dd.1coords
tail -n +2 "$PAIRS" | while IFS=$'\t' read -r PID QA RA; do
    DD=$OUT/$PID
    [ -s "$DD/dd.report" ] && continue
    TMP=$(mktemp -d "$WORK/tmp.${PID}.XXXXXX")
    for ACC in "$QA" "$RA"; do
        SRC=$GENOMES/${ACC}.fna
        if [ ! -s "$SRC" ]; then
            echo "[warning] missing genome $SRC" >&2
            rm -rf "$TMP"
            continue 2
        fi
        if grep -qm1 -P "\s" "$SRC"; then
            awk '/^>/{sub(/[ \t].*$/,"");print;next}{gsub(/[ \t]/,"");print}' \
                "$SRC" > "$TMP/${ACC}.fna"
        else
            ln -s "$SRC" "$TMP/${ACC}.fna"
        fi
    done
    mkdir -p "$DD"
    (cd "$DD" && dnadiff -p dd "$TMP/${RA}.fna" "$TMP/${QA}.fna" > /dev/null 2>&1) \
        || echo "[closed_top100] $PID: dnadiff failed"
    rm -rf "$TMP"
    # keep only report + 1coords (delta files are large)
    rm -f "$DD"/dd.delta "$DD"/dd.mdelta "$DD"/dd.1delta "$DD"/dd.mcoords \
          "$DD"/dd.qdiff "$DD"/dd.rdiff "$DD"/dd.snps 2>/dev/null || true
done

# 2. Collect junction coordinates from Syn2b and dnadiff
$PY "$ROOT/scripts/gtdb50k/collect_junction_coordinates.py" \
    --pairs     "$PAIRS" \
    --tgt-cache "$TGT_CACHE" \
    --dnadiff   "$OUT" \
    --outdir    "$WORK/top100_results" \
    --workers   16

# 3. Compare positions
$PY "$ROOT/scripts/gtdb50k/compare_junction_positions.py" \
    --coords "$WORK/top100_results/junction_coordinates.tsv" \
    --out    "$WORK/top100_results/position_agreement.tsv"

echo "[closed_top100] results in $WORK/top100_results/"
