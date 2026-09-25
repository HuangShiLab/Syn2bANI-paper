#!/bin/bash
#SBATCH --job-name=s2b_packed_db
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/s2b_packed_db_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/s2b_packed_db_%j.err

set -euo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
S2B=/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani
REFS=$BASE/gtdb_r207_references.txt
DBD=$BASE/gtdb_r207_s2b_packed.db
TMP=/tmp/s2b_db_chunks_$$
mkdir -p $TMP
rm -rf $DBD

# first chunk: db build; remaining chunks: db add
split -l 2000 "$REFS" $TMP/chunk_
FIRST=$(ls $TMP/chunk_* | head -1)
echo "=== db build: $FIRST ==="
/usr/bin/time -v $S2B db build -o "$DBD" -t 32 $(cat "$FIRST") 2>> $BASE/logs/s2b_packed_db.time
for CHUNK in $(ls $TMP/chunk_* | tail -n +2); do
    echo "=== db add: $CHUNK ==="
    $S2B db add -d "$DBD" -t 32 -p $(cat "$CHUNK")
done
rm -rf $TMP

echo "=== verify ==="
$S2B db list -d "$DBD" | wc -l
du -sh "$DBD"
echo BUILD_DB_DONE
