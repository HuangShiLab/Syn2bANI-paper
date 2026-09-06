#!/bin/bash
#SBATCH --job-name=gtdb_build_skani_sketches
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_build_skani_sketches_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_build_skani_sketches_%j.err

set -euo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
SKANI=/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani
LOGS=$BASE/logs
mkdir -p $LOGS

REFLIST=$BASE/gtdb_r207_references.txt
SKANI_DB=$BASE/gtdb_r207_skani_sketches

# skani sketch requires output dir to NOT exist
rm -rf $SKANI_DB

echo "=== Building skani sketch DB ==="
/usr/bin/time -v $SKANI sketch -t 32 --separate-sketches -o $SKANI_DB -l $REFLIST \
    > $LOGS/skani_sketch_db.out 2> $LOGS/skani_sketch_db.err

echo "Done"
