#!/bin/bash
# download_gtdb_r207.sh
# Download GTDB-R207 representative genomes for Syn2bANI benchmarking
# Target: Mac Studio with 2TB storage

set -euo pipefail

DATA_ROOT="${GTDB_ROOT:-$HOME/data/gtdb-r207}"
GENOMES_DIR="$DATA_ROOT/genomes"
METADATA_DIR="$DATA_ROOT/metadata"
LOG_FILE="$DATA_ROOT/download.log"

mkdir -p "$GENOMES_DIR" "$METADATA_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "GTDB-R207 Download Script"
echo "Data root: $DATA_ROOT"
echo "Start time: $(date)"
echo "========================================"

# --- Step 1: Download metadata ---
echo "[1/5] Downloading GTDB-R207 metadata..."

METADATA_URL="https://data.gtdb.ecogenomic.org/releases/release207/207.0"
wget -q --show-progress -c "$METADATA_URL/bac120_metadata_r207.tsv" \
    -O "$METADATA_DIR/bac120_metadata_r207.tsv" || true
wget -q --show-progress -c "$METADATA_URL/ar53_metadata_r207.tsv" \
    -O "$METADATA_DIR/ar53_metadata_r207.tsv" || true

# Taxonomy tree (useful for validation)
wget -q --show-progress -c "$METADATA_URL/bac120_r207.tsv" \
    -O "$METADATA_DIR/bac120_taxonomy_r207.tsv" || true
wget -q --show-progress -c "$METADATA_URL/ar53_r207.tsv" \
    -O "$METADATA_DIR/ar53_taxonomy_r207.tsv" || true

echo "[1/5] Metadata download complete."

# --- Step 2: Download representative genomes ---
echo "[2/5] Downloading GTDB-R207 representative genomes..."

# GTDB provides genome URLs in the metadata file
# Alternative: use NCBI datasets CLI (faster, recommended)

if command -v datasets &> /dev/null; then
    echo "Using NCBI datasets CLI..."
    
    # Extract accession list from metadata
    python3 -c "
import pandas as pd
bac = pd.read_csv('$METADATA_DIR/bac120_metadata_r207.tsv', sep='\t', low_memory=False)
ar = pd.read_csv('$METADATA_DIR/ar53_metadata_r207.tsv', sep='\t', low_memory=False)
df = pd.concat([bac, ar])
# Extract NCBI accession from column (varies by GTDB version)
# Common column names: 'ncbi_genbank_assembly_accession', 'accession'
acc_col = None
for col in ['ncbi_genbank_assembly_accession', 'accession', 'genome_accession']:
    if col in df.columns:
        acc_col = col
        break
if acc_col:
    df[acc_col].dropna().to_csv('$METADATA_DIR/accessions.txt', index=False, header=False)
    print(f'Extracted {len(df)} accessions')
else:
    print('WARNING: Could not find accession column. Available columns:')
    print(list(df.columns))
"
    
    # Download in batches using datasets
    if [ -f "$METADATA_DIR/accessions.txt" ]; then
        total=$(wc -l < "$METADATA_DIR/accessions.txt")
        batch_size=500
        echo "Downloading $total genomes in batches of $batch_size..."
        
        split -l "$batch_size" "$METADATA_DIR/accessions.txt" "$METADATA_DIR/batch_"
        
        for batch in "$METADATA_DIR"/batch_*; do
            echo "Processing batch: $(basename $batch)"
            datasets download genome accession \
                --inputfile "$batch" \
                --filename "${batch}.zip" \
                --assembly-source RefSeq \
                --assembly-level complete,chromosome,scaffold \
                --dehydrated || true
            
            if [ -f "${batch}.zip" ]; then
                unzip -q "${batch}.zip" -d "${batch}_extracted" || true
                find "${batch}_extracted" -name "*.fna" -exec mv {} "$GENOMES_DIR/" \;
                rm -rf "${batch}.zip" "${batch}_extracted"
            fi
        done
    fi
else
    echo "NCBI datasets CLI not found. Installing..."
    echo "Run: conda install -c conda-forge ncbi-datasets-cli"
    echo "Or download from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
    
    # Fallback: direct tar download (slower, larger)
    echo "Falling back to direct GTDB tar download..."
    TAR_URL="https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz"
    wget -q --show-progress -c "$TAR_URL" -O "$DATA_ROOT/gtdb_genomes_reps_r207.tar.gz"
    
    echo "Extracting tar archive (this may take 30+ minutes)..."
    tar -xzf "$DATA_ROOT/gtdb_genomes_reps_r207.tar.gz" -C "$GENOMES_DIR/"
    echo "You may delete the tar archive after extraction to save space:"
    echo "  rm $DATA_ROOT/gtdb_genomes_reps_r207.tar.gz"
fi

echo "[2/5] Genome download complete."

# --- Step 3: Validate downloads ---
echo "[3/5] Validating downloaded genomes..."

genome_count=$(find "$GENOMES_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" | wc -l)
echo "Found $genome_count genome files"

# Check for empty or corrupted files
find "$GENOMES_DIR" -name "*.fna" -size 0 -delete
find "$GENOMES_DIR" -name "*.fasta" -size 0 -delete

echo "[3/5] Validation complete."

# --- Step 4: Generate manifest ---
echo "[4/5] Generating genome manifest..."

python3 -c "
import os
import pandas as pd
from pathlib import Path

genomes_dir = '$GENOMES_DIR'
metadata_dir = '$METADATA_DIR'

# Find all genome files
files = []
for ext in ['*.fna', '*.fasta', '*.fa']:
    files.extend(Path(genomes_dir).glob(ext))

print(f'Found {len(files)} genome files')

# Build manifest
manifest = []
for f in files:
    size = f.stat().st_size
    # Try to extract genome ID from filename
    # GTDB format: GB_GCA_000xyz.1_genomic.fna or RS_GCF_000xyz.1_genomic.fna
    stem = f.stem.replace('_genomic', '')
    manifest.append({
        'genome_id': stem,
        'path': str(f),
        'file_size': size,
    })

manifest_df = pd.DataFrame(manifest)
manifest_df.to_csv(f'{metadata_dir}/manifest.tsv', sep='\t', index=False)
print(f'Manifest written: {metadata_dir}/manifest.tsv')
print(f'Total genomes: {len(manifest_df)}')
"

echo "[4/5] Manifest generation complete."

# --- Step 5: Quick stats ---
echo "[5/5] Computing genome statistics..."

if command -v seqkit &> /dev/null; then
    seqkit stats \
        $(find "$GENOMES_DIR" -name "*.fna" | head -100) \
        > "$METADATA_DIR/quick_stats_100samples.txt"
    echo "Quick stats (100 samples) saved to: $METADATA_DIR/quick_stats_100samples.txt"
else
    echo "seqkit not found. Install with: conda install -c bioconda seqkit"
fi

echo "========================================"
echo "Download Summary"
echo "End time: $(date)"
echo "Genomes: $genome_count"
echo "Manifest: $METADATA_DIR/manifest.tsv"
echo "========================================"
echo ""
echo "Next steps:"
echo "  1. Run: python3 scripts/sample_gtdb_pairs.py"
echo "  2. Run: python3 scripts/run_benchmark_matrix.py"
