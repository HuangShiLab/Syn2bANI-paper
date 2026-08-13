#!/bin/bash
# Efficiency benchmark: syn2bani v8 vs skani vs FastANI
# Machine: Mac Studio, 16 threads, all tools at full core count.
set -u

OUT=/Users/macstudio/Downloads/Syn2bANI-paper/results/efficiency_v8
LOGS=$OUT/logs
SYN=/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani
SKANI=/Users/macstudio/.cargo/bin/skani
FASTANI=/opt/homebrew/bin/fastani
REAL=/Users/macstudio/Downloads/Syn2bANI/prototype/realbench/genomes
DRAFT=/Users/macstudio/Downloads/Syn2bANI/prototype/draftbench/drafts
THREADS=16
REPS=3

mkdir -p "$LOGS" "$OUT/lists"

# Fixed nested genome order: E. coli completes first, then other
# Enterobacteriaceae, then real draft E. coli assemblies.
POOL=(
  $REAL/Ecoli_K12_MG1655.fasta
  $REAL/Ecoli_O157H7_Sakai.fasta
  $REAL/Ecoli_K12_W3110.fasta
  $REAL/Ecoli_CFT073.fasta
  $REAL/Ecoli_UTI89.fasta
  $REAL/Ecoli_BL21_DE3.fasta
  $REAL/Shigella_flexneri_301.fasta
  $REAL/Shigella_sonnei_Ss046.fasta
  $REAL/Escherichia_fergusonii.fasta
  $REAL/Salmonella_Typhimurium_LT2.fasta
  $REAL/Salmonella_Typhi_CT18.fasta
  $REAL/Klebsiella_pneumoniae_MGH78578.fasta
  $REAL/Enterobacter_cloacae_13047.fasta
  $REAL/Citrobacter_rodentium.fasta
  $DRAFT/GCA_001075925.fasta
  $DRAFT/GCA_001077875.fasta
  $DRAFT/GCA_001283205.fasta
  $DRAFT/GCA_001283245.fasta
  $DRAFT/GCA_001283605.fasta
  $DRAFT/GCA_001283865.fasta
  $DRAFT/GCA_001284145.fasta
  $DRAFT/GCA_001284645.fasta
)
TOTAL=${#POOL[@]}

# Manifest
printf "n_genomes\tgenome\n" > "$OUT/genome_subsets.tsv"
for n in 2 5 10 15 $TOTAL; do
  for ((i=0;i<n;i++)); do
    printf "%s\t%s\n" "$n" "$(basename "${POOL[$i]}")" >> "$OUT/genome_subsets.tsv"
  done
done

# Result tables
printf "tool\tmode\tn_genomes\tn_pairs\trep\twall_s\tpeak_rss_mb\n" > "$OUT/runtime_scaling.tsv"
printf "tool\tn_genomes\trep\twall_s\ttotal_size_kb\tartifact_path\n" > "$OUT/sketch_benchmark.tsv"

# run_and_log <tool> <mode> <n> <rep> <logfile> <cmd...>
# Appends one row to runtime_scaling.tsv; tool stdout/stderr -> logfile.
run_and_log() {
  local tool=$1 mode=$2 n=$3 rep=$4 log=$5; shift 5
  /usr/bin/time -l "$@" > "$log" 2> "$log.time"
  local wall rss
  wall=$(awk '$2=="real" {print $1; exit}' "$log.time")
  rss=$(awk '/maximum resident set size/ {print $1; exit}' "$log.time")
  local rss_mb
  rss_mb=$(awk -v b="$rss" 'BEGIN{printf "%.1f", b/1048576}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tool" "$mode" "$n" "$((n*n))" "$rep" "$wall" "$rss_mb" >> "$OUT/runtime_scaling.tsv"
}

# run_sketch <tool> <n> <rep> <logfile> <artifact_path> <cmd...>
run_sketch() {
  local tool=$1 n=$2 rep=$3 log=$4 art=$5; shift 5
  /usr/bin/time -l "$@" > "$log" 2> "$log.time"
  local wall size
  wall=$(awk '$2=="real" {print $1; exit}' "$log.time")
  size=$(du -sk "$art" | cut -f1)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tool" "$n" "$rep" "$wall" "$size" "$art" >> "$OUT/sketch_benchmark.tsv"
}

for n in 2 5 10 15 $TOTAL; do
  SUB=("${POOL[@]:0:$n}")
  FLIST=$OUT/lists/fasta_n$n.txt
  printf "%s\n" "${SUB[@]}" > "$FLIST"

  # --- syn2bani sketch (3 reps; keep rep1 artifacts for reuse/size) ---
  for rep in 1 2 3; do
    S2DIR=$OUT/sketches/syn2bani_n$n
    if [ $rep -eq 1 ]; then rm -rf "$S2DIR"; else rm -rf "$S2DIR.tmp"; fi
    TGT=$S2DIR; [ $rep -ne 1 ] && TGT=$S2DIR.tmp
    run_sketch syn2bani "$n" "$rep" "$LOGS/syn2bani_sketch_n${n}_r${rep}.log" "$TGT" \
      "$SYN" sketch --enzymes BcgI,AlfI,AloI,FalI -t 0 -o "$TGT" "${SUB[@]}"
    [ $rep -ne 1 ] && rm -rf "$S2DIR.tmp"
  done
  S2LIST=$OUT/lists/s2ba_n$n.txt
  ls "$OUT/sketches/syn2bani_n$n"/*.s2ba > "$S2LIST"

  # --- skani sketch (3 reps; keep rep1 db) ---
  for rep in 1 2 3; do
    SKDIR=$OUT/sketches/skani_n$n
    TGT=$SKDIR; [ $rep -ne 1 ] && TGT=$SKDIR.tmp
    run_sketch skani "$n" "$rep" "$LOGS/skani_sketch_n${n}_r${rep}.log" "$TGT" \
      "$SKANI" sketch -t $THREADS -o "$TGT" "${SUB[@]}"
    [ $rep -ne 1 ] && rm -rf "$SKDIR.tmp"
  done
  SKLIST=$OUT/lists/skani_n$n.txt
  ls "$OUT/sketches/skani_n$n"/*.sketch > "$SKLIST"

  # --- ANI modes, 3 reps each ---
  for rep in 1 2 3; do
    run_and_log syn2bani ani_fasta "$n" "$rep" "$LOGS/syn2bani_ani_fasta_n${n}_r${rep}.log" \
      "$SYN" ani --ql "$FLIST" --rl "$FLIST" -t 0 -o "$LOGS/syn2bani_ani_fasta_n${n}_r${rep}.tsv"

    run_and_log syn2bani ani_sketches "$n" "$rep" "$LOGS/syn2bani_ani_sk_n${n}_r${rep}.log" \
      "$SYN" ani --ql "$S2LIST" --rl "$S2LIST" -t 0 -o "$LOGS/syn2bani_ani_sk_n${n}_r${rep}.tsv"

    run_and_log skani dist "$n" "$rep" "$LOGS/skani_dist_n${n}_r${rep}.log" \
      "$SKANI" dist -t $THREADS --ql "$SKLIST" --rl "$SKLIST" -o "$LOGS/skani_dist_n${n}_r${rep}.tsv"

    run_and_log fastani all_vs_all "$n" "$rep" "$LOGS/fastani_n${n}_r${rep}.log" \
      "$FASTANI" --ql "$FLIST" --rl "$FLIST" -t $THREADS -o "$LOGS/fastani_n${n}_r${rep}.tsv"
  done
  echo "done n=$n"
done
echo "ALL DONE"
