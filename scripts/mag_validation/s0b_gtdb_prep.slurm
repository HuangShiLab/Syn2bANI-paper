#!/bin/bash
#SBATCH --job-name=mv_gtdbprep
#SBATCH --partition=intel
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/gtdbprep_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/gtdbprep_%A.err
# s0b: sketch all GTDB r207 genomes_all (65,703) into a skani search DB,
#      and build first-record-id map (for syn2bani output name resolution).
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"

OUT=$WORK/refs/gtdb
mkdir -p "$OUT" "$WORK/logs"
cd "$OUT"

if [ ! -s genome_paths.txt ]; then
    ls "$GTDB_DIR"/*.fna | sort > genome_paths.txt
fi
echo "[gtdbprep] $(wc -l < genome_paths.txt) GTDB genomes"

# skani search DB (idempotent via sentinel)
if [ ! -f "$OUT/db.DONE" ]; then
    rm -rf "$OUT/db"
    $SKANI sketch -l genome_paths.txt -o "$OUT/db" -t "$SLURM_CPUS_PER_TASK"
    touch "$OUT/db.DONE"
fi

# first record id map: genome path -> first fasta header token
if [ ! -s firstid_map.tsv ]; then
    $PY "$SCRIPTS/fasta_index.py" <(awk -F/ '{n=$NF; sub(/\.fna$/,"",n); print n"\t"$0}' genome_paths.txt) \
        genome_lengths.tsv firstid_map.tsv
fi
echo "[gtdbprep] done"
