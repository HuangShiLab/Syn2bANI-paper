#!/bin/bash
# run_struct_discordant.sh — run syn2bani struct on the top ANI-synteny
# discordant pairs from the held-out and high-ANI test sets.
set -euo pipefail

WORK=${GTDB50K_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${GTDB50K_GENOMES:-/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all}
TOOLS=${TOOLS:-/lustre1/g/aos_shihuang/tools}
S2B=${S2B:-$TOOLS/syn2bani/syn2bani}
PY=${PY:-$TOOLS/python3/bin/python3}
N=${N:-100}

OUT=$WORK/struct_discordant_top${N}
mkdir -p "$OUT" "$WORK/logs"

# Combine discordant pair lists
cat "$WORK/discordant_ani95_syn98.tsv" "$WORK/discordant_high_ani_test_syn98.tsv" 2>/dev/null | \
    awk -F'\t' 'NR==1{next}{print $1}' | sort -u | head -n "$N" > "$OUT/top_pairs.txt"

NP=$(wc -l < "$OUT/top_pairs.txt")
echo "Running syn2bani struct on $NP discordant pairs"

while IFS= read -r PID; do
    QA=${PID%%__*}
    RA=${PID##*__}
    [ -s "$OUT/${PID}.paf" ] && continue
    QF=$GENOMES/${QA}.fna
    RF=$GENOMES/${RA}.fna
    if [ ! -s "$QF" ] || [ ! -s "$RF" ]; then
        echo "[struct] $PID: missing FASTA"
        continue
    fi
    "$S2B" struct "$QF" "$RF" > "$OUT/${PID}.paf" 2> "$OUT/${PID}.err" \
        || echo "[struct] FAIL $PID"
done < "$OUT/top_pairs.txt"

echo "Done: $(ls "$OUT"/*.paf 2>/dev/null | wc -l) PAF files"
