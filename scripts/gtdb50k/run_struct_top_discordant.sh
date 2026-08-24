#!/bin/bash
# run_struct_top_discordant.sh — run syn2bani struct on top GTDB discordant high-ANI pairs.
# Input: results/gtdb50k/top_discordant_high_ani_pairs.tsv
# Output: $WORK/struct_top_discordant/{pairid}.tsv (PAF chains) and {pairid}_ani.tsv (verbose ani)
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
S2B=${SYN2BANI:-/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani}
PAIRS=${1:-$WORK/top_discordant_high_ani_pairs.tsv}
OUTDIR="$WORK/struct_top_discordant"
ENZYMES="BcgI,AlfI,AloI,FalI"

mkdir -p "$OUTDIR"

# Skip header
tail -n +2 "$PAIRS" | while IFS=$'\t' read -r RANK PID QA RA _; do
    PAF_OUT="$OUTDIR/${PID}.tsv"
    ANI_OUT="$OUTDIR/${PID}_ani.tsv"
    [ -s "$PAF_OUT" ] && [ -s "$ANI_OUT" ] && continue

    QF="$GENOMES/${QA}.fna"
    RF="$GENOMES/${RA}.fna"
    if [ ! -s "$QF" ] || [ ! -s "$RF" ]; then
        echo "[skip] missing genome for $PID" >&2
        continue
    fi

    "$S2B" struct "$QF" "$RF" --enzymes "$ENZYMES" --paf --rearrangement --indel -o "$PAF_OUT" >/dev/null 2>&1 || true
    "$S2B" ani "$QF" "$RF" --enzymes "$ENZYMES" --verbose > "$ANI_OUT" 2>/dev/null || true

    if [ -s "$PAF_OUT" ]; then
        echo "[done] $PID"
    else
        echo "[fail] $PID" >&2
    fi
done

echo "[struct_top_discordant] finished"
