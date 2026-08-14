#!/bin/bash
# Setup: genome lists (nested subsets + held-out queries) and CLI smoke tests.
source "$(dirname "$0")/common.sh"

# --- deterministic master permutation (fixed seed) ---
if [ ! -s "$BASE/lists/master_perm.txt" ]; then
    ls "$GENOMES"/*.fna | shuf --random-source=<(yes 42) > "$BASE/lists/master_perm.txt"
fi
head -100  "$BASE/lists/master_perm.txt" > "$BASE/lists/n100.txt"
head -500  "$BASE/lists/master_perm.txt" > "$BASE/lists/n500.txt"
head -2000 "$BASE/lists/master_perm.txt" > "$BASE/lists/n2000.txt"
head -5000 "$BASE/lists/master_perm.txt" > "$BASE/lists/n5000.txt"
sed -n '5001,5100p' "$BASE/lists/master_perm.txt" > "$BASE/lists/queries100.txt"
wc -l "$BASE"/lists/*.txt

# --- CLI surface dumps for the record ---
{
  echo "== syn2bani version =="; "$S2B" --version
  echo "== skani version =="; "$SKANI" --version
  for sub in ani dist sketch search triangle db; do
      echo "== syn2bani $sub --help =="; "$S2B" $sub --help 2>&1 || true
  done
} > "$BASE/logs/cli_help.txt" 2>&1

# --- smoke tests on 3 genomes ---
SMOKE="$BASE/out/smoke"; mkdir -p "$SMOKE"
head -3 "$BASE/lists/n100.txt" > "$SMOKE/g3.txt"
mapfile -t G < "$SMOKE/g3.txt"

echo "--- smoke: sketch 4-enzyme panel"
"$S2B" sketch "${G[@]}" -o "$SMOKE/s2ba4" --enzymes "$ENZ4" -p -t 4
ls "$SMOKE/s2ba4"
echo "--- smoke: sketch default (BcgI)"
"$S2B" sketch "${G[@]}" -o "$SMOKE/s2ba_bcgi" -p -t 4
echo "--- smoke: dist 4-enzyme (2 positional fastas)"
"$S2B" dist "${G[0]}" "${G[1]}" --enzymes "$ENZ4" -p -t 4 | tee "$SMOKE/dist_enz4.tsv"
echo "--- smoke: dist default enzyme AloI,BslFI"
"$S2B" dist "${G[0]}" "${G[1]}" -p -t 4 | tee "$SMOKE/dist_default.tsv"
echo "--- smoke: dist --mash-ani (GBRT v7 output)"
"$S2B" dist "${G[0]}" "${G[1]}" --enzymes "$ENZ4" --mash-ani -p -t 4 | tee "$SMOKE/dist_gbrt.tsv"
echo "--- smoke: dist positional split check (3 fastas: which is ref?)"
"$S2B" dist "${G[0]}" "${G[1]}" "${G[2]}" --enzymes "$ENZ4" -p -t 4 | tee "$SMOKE/dist_3pos.tsv"
echo "--- smoke: triangle edge list"
"$S2B" triangle "${G[@]}" --edge-list -p -t 4 | tee "$SMOKE/triangle.tsv"
echo "--- smoke: search (fasta query vs 4-enzyme db)"
"$S2B" search "${G[0]}" "$SMOKE/s2ba4" -p -t 4 --min-ani 0.5 | tee "$SMOKE/search.tsv"
echo "--- smoke: search vs BcgI db"
"$S2B" search "${G[0]}" "$SMOKE/s2ba_bcgi" -p -t 4 --min-ani 0.5 | tee "$SMOKE/search_bcgi.tsv"
echo "--- smoke: ani (validated path)"
"$S2B" ani "${G[0]}" "${G[1]}" | tee "$SMOKE/ani.tsv"
echo "--- smoke: skani sketch"
"$SKANI" sketch -l "$SMOKE/g3.txt" -o "$SMOKE/skani_db" -t 4
ls "$SMOKE/skani_db"
echo "SMOKE DONE"
