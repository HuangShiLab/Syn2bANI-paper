#!/usr/bin/env bash
# Download raw FASTQs for all SynTracker isolates from ENA.
# Run this on an I/O node (hpc2021-io1/2) because it is network-bound
# and may take >15 min.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
MANIFEST_DIR="${BASE}/manifests"
READS_DIR="${BASE}/reads"
mkdir -p "${READS_DIR}"

for manifest in "${MANIFEST_DIR}"/fastq_manifest_*.tsv; do
  [ -e "${manifest}" ] || continue
  echo "[download] $(basename "${manifest}")"
  tail -n +2 "${manifest}" | while IFS=$'\t' read -r isolate species sra layout url1 url2 bytes1 bytes2; do
    [ -z "${sra}" ] && continue
    d="${READS_DIR}/${isolate}"
    mkdir -p "${d}"
    for url in "${url1}" "${url2}"; do
      [ -z "${url}" ] && continue
      fname=$(basename "${url}")
      target="${d}/${fname}"
      if [ -s "${target}" ]; then
        echo "  ${isolate}/${fname} exists, skipping"
        continue
      fi
      echo "  downloading ${fname}"
      wget -q -c "http://${url}" -O "${target}.tmp" && mv "${target}.tmp" "${target}"
    done
  done
done

echo "Reads staged in ${READS_DIR}"
