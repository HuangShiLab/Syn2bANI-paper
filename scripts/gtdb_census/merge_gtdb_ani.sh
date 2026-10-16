#!/bin/bash
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
OUT=$R/gtdb_r207_ani_merged_sorted.tsv.gz
TMP=$W/ani_merge_tmp
mkdir -p "$TMP" "$R"
trap 'rm -rf "$TMP"' EXIT
# Existing per-cluster ANI files plus chunked completion file.
find "$W/ani" -maxdepth 1 -type f -name '*.ani.tsv.gz' -print0 > "$TMP/inputs.list"
printf '%s\0' "$R/gtdb_ani_remaining_chunked.tsv.gz" >> "$TMP/inputs.list"
# Do not run concurrent zcat processes with a shared stdout: large stdio
# writes can interleave and corrupt the merged stream.
(
  printf 'genome_A\tgenome_B\tskani_ani\n'
  while IFS= read -r -d '' input; do zcat "$input"; done < "$TMP/inputs.list"
) | awk -F '\t' '$1!="" && $2!="" && $1!=$2 {print $1 "\t" $2 "\t" $3}' |
  LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2 -u -S 64G --parallel=16 -T "$TMP" |
  pigz -p 16 > "$OUT"
echo GTDB-ANI-MERGE-DONE
