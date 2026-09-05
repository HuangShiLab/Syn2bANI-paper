#!/bin/bash
#SBATCH --job-name=gtdb_build_sketches
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_build_sketches_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_build_sketches_%j.err

set -euo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
GENOMES=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
S2B=/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani
SKANI=/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani
LOGS=$BASE/logs
mkdir -p $LOGS

REFLIST=$BASE/gtdb_r207_references.txt

# Syn2bANI sketch database
S2B_DB=$BASE/gtdb_r207_s2b_sketches
mkdir -p $S2B_DB
echo "=== Building Syn2bANI sketch DB ==="
/usr/bin/time -v python3 $BASE/build_s2b_gtdb_sketches.py $REFLIST $S2B_DB 500 32 \
    > $LOGS/s2b_sketch_db.out 2> $LOGS/s2b_sketch_db.err

# skani sketch database
SKANI_DB=$BASE/gtdb_r207_skani_sketches
mkdir -p $SKANI_DB
echo "=== Building skani sketch DB ==="
/usr/bin/time -v $SKANI sketch -t 32 -o $SKANI_DB -l $REFLIST \
    > $LOGS/skani_sketch_db.out 2> $LOGS/skani_sketch_db.err

echo "Done"
