#!/bin/bash
# Common environment for db_scale benchmarks. Source this; do not run directly.
set -euo pipefail
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani   # provides skani 0.3.2
export S2B=/lustre1/g/aos_shihuang/Syn2bANI-bench/target/release/syn2bani
export SKANI=skani
export GENOMES=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
export BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale
export ENZ4="BcgI,AlfI,AloI,FalI"
export THREADS=${THREADS:-16}
export TIME=/usr/bin/time

mkdir -p "$BASE"/{lists,logs,out,sketches,scripts}

# run_timed <tsv> <tool> <phase> <n> <rep> <logfile> <cmd...>
# Appends: tool phase n rep wall_s max_rss_kb exit_status ; full /usr/bin/time -v in logfile
run_timed() {
    local tsv=$1 tool=$2 phase=$3 n=$4 rep=$5 log=$6; shift 6
    local status=0
    "$TIME" -v -o "$log" "$@" || status=$?
    local wall rss
    wall=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "$log" | tail -1)
    rss=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$log" | tail -1)
    # convert h:mm:ss / m:ss to seconds
    local secs
    secs=$(echo "$wall" | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$phase" "$n" "$rep" "$secs" "$rss" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$tsv"
    return $status
}
