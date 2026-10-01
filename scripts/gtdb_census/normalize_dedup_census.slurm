#!/bin/bash
#SBATCH --job-name=gtdb-dedup
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_dedup_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_dedup_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
OUT=$R/gtdb_r207_effective_normalized_dedup.tsv.gz
TMP=$W/sort_tmp
cd $L/Syn2b-paper
mkdir -p "$TMP" "$R" "$W/logs"
trap 'rm -rf "$TMP"' EXIT
python3 scripts/gtdb_census/build_effective_accessions.py --workdir "$W"
(
  printf 'genome_A\tgenome_B\tbreakpoints\tscj_distance\tbreakpoint_density\tinverted_fraction\traw_inverted_fraction\torientation_mismatches\torientation_mismatches_raw\torientation_uninformative\tobservable_fraction\tobservable_adjacencies\tstructural\tshared_tags\trepeats_dropped\tlandmarks_collapsed\tcircular\tlegacy_adjacency\n'
  python3 scripts/gtdb_census/normalize_dedup_census.py --workdir "$W" --stats "$W/normalization_stats.tsv" \
    | LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2 -S 48G --parallel=16 -T "$TMP" \
    | awk -F '\t' 'BEGIN{OFS="\t"} NR==1 {a=$1;b=$2;print;next} $1!=a || $2!=b {a=$1;b=$2;print}'
) | pigz -p 8 > "$OUT"
python3 scripts/gtdb_census/qc_normalized_census.py --table "$OUT" --workdir "$W" --out "$R"
echo GTDB-EFFECTIVE-DEDUP-DONE
