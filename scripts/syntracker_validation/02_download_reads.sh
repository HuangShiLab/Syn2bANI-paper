#!/usr/bin/env bash
# Download raw FASTQs for all SynTracker isolates via NCBI SRA Toolkit.
# Run this on an I/O node (hpc2021-io1/2) because it is network-bound
# and may take several hours.
set -uo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
MANIFEST_DIR="${BASE}/manifests"
READS_DIR="${BASE}/reads"
mkdir -p "${READS_DIR}"

COMPRESS="${COMPRESS_CMD:-gzip}"
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -p 4"
fi

for manifest in "${MANIFEST_DIR}"/fastq_manifest_*.tsv; do
  [ -e "${manifest}" ] || continue
  echo "[download] $(basename "${manifest}")"
  tail -n +2 "${manifest}" | while IFS=$'\t' read -r isolate species sra layout url1 url2 bytes1 bytes2; do
    [ -z "${sra}" ] && continue
    d="${READS_DIR}/${isolate}"
    mkdir -p "${d}"

    # Skip if paired FASTQs already present.
    if [ -s "${d}/${sra}_1.fastq.gz" ]; then
      echo "  ${isolate}/${sra}_1.fastq.gz exists, skipping"
      continue
    fi

    echo "  prefetch ${sra}"
    prefetch -O "${d}" "${sra}" >/dev/null 2>&1 || {
      echo "WARN: prefetch failed for ${sra}" >&2
      continue
    }

    sra_file="${d}/${sra}/${sra}.sra"
    if [ ! -s "${sra_file}" ]; then
      echo "WARN: missing SRA file for ${sra}" >&2
      continue
    fi

    echo "  fasterq-dump ${sra}"
    fasterq-dump --outdir "${d}" "${sra_file}" >/dev/null 2>&1 || {
      echo "WARN: fasterq-dump failed for ${sra}" >&2
      continue
    }

    # Compress and tidy up.
    for fq in "${d}/${sra}.fastq" "${d}/${sra}_1.fastq" "${d}/${sra}_2.fastq"; do
      [ -e "${fq}" ] || continue
      ${COMPRESS} -f "${fq}"
    done
    rm -rf "${d}/${sra}.sra" "${d}/${sra}"
    echo "  done ${isolate}"
  done
done

echo "Reads staged in ${READS_DIR}"
