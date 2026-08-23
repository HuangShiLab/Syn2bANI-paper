#!/bin/bash
# Download GTDB-R207 non-representative genomes needed for the 95-97/97-100
# strata of the unified benchmark. Splits the accession list into chunks of
# CHUNK_SIZE and uses NCBI datasets CLI. Idempotent: skips chunks whose
# extracted FASTAs are already present.
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOME_DIR="$WORK/genomes_high_ani"
LIST="$WORK/high_ani_genomes.txt"
CHUNK_SIZE=${CHUNK_SIZE:-500}
DATASETS=${DATASETS:-/group/aos_shihuang/conda/bin/datasets}

mkdir -p "$GENOME_DIR"

# total accessions
N=$(wc -l < "$LIST")
NCHUNKS=$(( (N + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "Downloading $N genomes in $NCHUNKS chunks of $CHUNK_SIZE to $GENOME_DIR"

for i in $(seq 0 $((NCHUNKS - 1))); do
    START=$(( i * CHUNK_SIZE + 1 ))
    END=$(( (i + 1) * CHUNK_SIZE ))
    CHUNK_FILE="$WORK/chunk_${i}.txt"
    sed -n "${START},${END}p" "$LIST" > "$CHUNK_FILE"

    # count how many of this chunk are already downloaded
    MISSING=0
    while read -r ACC; do
        [ -z "$ACC" ] && continue
        if [ ! -s "$GENOME_DIR/${ACC}.fna" ]; then
            MISSING=$((MISSING + 1))
        fi
    done < "$CHUNK_FILE"

    if [ "$MISSING" -eq 0 ]; then
        echo "Chunk $i/$NCHUNKS already complete, skipping"
        rm -f "$CHUNK_FILE"
        continue
    fi

    echo "Chunk $i/$NCHUNKS: downloading $MISSING missing genomes..."
    TMP=$(mktemp -d "$WORK/dl_chunk_${i}.XXXXXX")
    ZIP="$TMP/chunk.zip"
    if $DATASETS download genome accession --inputfile "$CHUNK_FILE" --filename "$ZIP" --include genome > "$TMP/dl.log" 2>&1; then
        unzip -q "$ZIP" -d "$TMP/extracted" || { echo "  unzip failed chunk $i"; rm -rf "$TMP" "$CHUNK_FILE"; continue; }
        find "$TMP/extracted/ncbi_dataset/data" -name "*.fna" | while read -r FNA; do
            BASE=$(basename "$FNA" | sed 's/_.*$//')
            cp "$FNA" "$GENOME_DIR/${BASE}.fna"
        done
        echo "  chunk $i done"
    else
        echo "  download failed chunk $i (see $TMP/dl.log)"
    fi
    rm -rf "$TMP" "$CHUNK_FILE"
done

echo "Download complete. Present genomes: $(ls $GENOME_DIR/*.fna 2>/dev/null | wc -l) / $N"
