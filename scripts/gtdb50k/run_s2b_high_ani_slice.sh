#!/bin/bash
# run_s2b_high_ani_slice.sh <slice_id> — syn2bani ani --calibrate for high-ANI pairs.
# Output: $WORK/s2b_high_ani/{pairid}.tsv
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
S2B=${SYN2BANI:-/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani}
NSLICES=${NSLICES:-80}

SLICE=$1
PAIRS=$WORK/high_ani_candidates.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/s2b_high_ani"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r PID QA RA _; do
    OUT="$WORK/s2b_high_ani/${PID}.tsv"
    [ -s "$OUT" ] && continue
    QF="$GENOMES/${QA}.fna"
    RF="$GENOMES/${RA}.fna"
    if [ ! -s "$QF" ] || [ ! -s "$RF" ]; then
        continue
    fi
    "$S2B" ani "$QF" "$RF" --verbose --calibrate > "$OUT" 2>/dev/null || true
    # ensure at least a header exists on failure
    if [ ! -s "$OUT" ]; then
        echo -e "query\treference\tani\tani_uniform\taf_query\taf_reference\tstd_err\tani_cal\tsynteny_blocks\tsynteny_score\tbreakpoint_count\thet_shape\tretention\tani_from_loss\tani_from_hist\tenzyme_spread\tenzyme_chi2\tper_enzyme\tn_anchors\tn_chains\tn_tags\tmax_block_anchors\tmean_block_anchors\tflag\tani_gated\tgate\tani_upper95" > "$OUT"
    fi
done

echo "[high_ani s2b] slice $SLICE done"
