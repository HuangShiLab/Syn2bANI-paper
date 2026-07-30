#!/usr/bin/env bash
# Assemble one isolate from paired or single FASTQ files.
# Called by the SLURM array driver (03_assemble_array.sh).
set -euo pipefail

isolate="${1:-}"
reads_dir="${2:-}"
asm_dir="${3:-}"

if [ -z "${isolate}" ] || [ -z "${reads_dir}" ] || [ -z "${asm_dir}" ]; then
  echo "Usage: $0 <isolate> <reads_dir> <asm_outdir>"
  exit 1
fi

out_fna="${asm_dir}/${isolate}.fna"
if [ -s "${out_fna}" ]; then
  echo "[${isolate}] assembly exists, skipping"
  exit 0
fi

mkdir -p "${asm_dir}"
tmpdir="${asm_dir}/.tmp_${isolate}"
rm -rf "${tmpdir}"
mkdir -p "${tmpdir}"

r1=$(find "${reads_dir}/${isolate}" -maxdepth 1 -name '*_1.fastq.gz' | head -1)
r2=$(find "${reads_dir}/${isolate}" -maxdepth 1 -name '*_2.fastq.gz' | head -1)
se=$(find "${reads_dir}/${isolate}" -maxdepth 1 -name '*.fastq.gz' ! -name '*_1.fastq.gz' ! -name '*_2.fastq.gz' | head -1)

CPUS="${SLURM_CPUS_PER_TASK:-8}"
MEM_GB="${ASSEMBLE_MEM_GB:-32}"

if command -v shovill >/dev/null 2>&1; then
  echo "[${isolate}] assembling with shovill (${CPUS} cpus, ${MEM_GB}G mem)"
  if [ -n "${r1}" ] && [ -n "${r2}" ]; then
    shovill --R1 "${r1}" --R2 "${r2}" --outdir "${tmpdir}" --cpus "${CPUS}" --ram "${MEM_GB}" --force --minlen 200
  elif [ -n "${se}" ]; then
    shovill --SE "${se}" --outdir "${tmpdir}" --cpus "${CPUS}" --ram "${MEM_GB}" --force --minlen 200
  else
    echo "[${isolate}] ERROR: no FASTQ found in ${reads_dir}/${isolate}" >&2
    exit 1
  fi
  cp "${tmpdir}/contigs.fa" "${out_fna}"
elif command -v spades.py >/dev/null 2>&1; then
  echo "[${isolate}] assembling with SPAdes (${CPUS} cpus, ${MEM_GB}G mem)"
  if [ -n "${r1}" ] && [ -n "${r2}" ]; then
    spades.py -1 "${r1}" -2 "${r2}" -o "${tmpdir}" -t "${CPUS}" -m "${MEM_GB}" --careful
  elif [ -n "${se}" ]; then
    spades.py -s "${se}" -o "${tmpdir}" -t "${CPUS}" -m "${MEM_GB}" --careful
  else
    echo "[${isolate}] ERROR: no FASTQ found in ${reads_dir}/${isolate}" >&2
    exit 1
  fi
  # Keep contigs >= 200 bp, rename to isolate
  python3 - <<PY
from pathlib import Path
seqs, name, buf = [], None, []
for line in Path("${tmpdir}/contigs.fasta").read_text().splitlines():
    if line.startswith('>'):
        if name and sum(len(s) for s in buf) >= 200:
            seqs.append((name, ''.join(buf)))
        name = line[1:].split()[0]
        buf = []
    else:
        buf.append(line.strip())
if name and sum(len(s) for s in buf) >= 200:
    seqs.append((name, ''.join(buf)))
with open("${out_fna}", 'w') as fh:
    for i, (_, seq) in enumerate(seqs, 1):
        fh.write(f">${isolate}_contig_{i}\n{seq}\n")
PY
else
  echo "[${isolate}] ERROR: neither shovill nor spades.py found" >&2
  exit 1
fi

rm -rf "${tmpdir}"

# Normalize the first contig header to the isolate name so that Syn2bANI's
# genome_id output can be trivially mapped back to the isolate.
python3 - <<PY
from pathlib import Path
lines = Path("${out_fna}").read_text().splitlines()
out = []
first = True
for line in lines:
    if line.startswith('>'):
        if first:
            out.append('>${isolate}')
            first = False
        else:
            out.append(line)
    else:
        out.append(line)
Path("${out_fna}").write_text('\n'.join(out) + '\n')
PY

echo "[${isolate}] wrote ${out_fna}"
