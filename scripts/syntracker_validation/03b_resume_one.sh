#!/usr/bin/env bash
# Resume an interrupted SPAdes assembly from its checkpoint tmpdir.
# Used when 03_assemble_one.sh hits the SLURM time limit: the tmpdir keeps
# SPAdes checkpoints, so --continue resumes from the last completed stage
# instead of restarting the whole assembly.
set -euo pipefail
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate spades 2>/dev/null || true
fi
isolate="$1"; asm_dir="$2"
out_fna="${asm_dir}/${isolate}.fna"
tmpdir="${asm_dir}/.tmp_${isolate}"
[ -s "${out_fna}" ] && { echo "[${isolate}] exists, skip"; exit 0; }
[ -d "${tmpdir}" ] || { echo "[${isolate}] no tmpdir to resume"; exit 1; }
CPUS="${SLURM_CPUS_PER_TASK:-16}"; MEM_GB="${ASSEMBLE_MEM_GB:-64}"
spades.py --restart-from last -o "${tmpdir}" -t "${CPUS}" -m "${MEM_GB}"
python3 - <<PY
from pathlib import Path
seqs, name, buf = [], None, []
for line in Path("${tmpdir}/contigs.fasta").read_text().splitlines():
    if line.startswith(">"):
        if name and sum(len(s) for s in buf) >= 200:
            seqs.append((name, "".join(buf)))
        name = line[1:].split()[0]; buf = []
    else:
        buf.append(line.strip())
if name and sum(len(s) for s in buf) >= 200:
    seqs.append((name, "".join(buf)))
with open("${out_fna}", "w") as fh:
    for i, (_, seq) in enumerate(seqs, 1):
        fh.write(f">${isolate}_contig_{i}\n{seq}\n")
PY
rm -rf "${tmpdir}"
python3 - <<PY
from pathlib import Path
lines = Path("${out_fna}").read_text().splitlines()
out, first = [], True
for line in lines:
    if line.startswith(">"):
        out.append(">${isolate}" if first else line); first = False
    else:
        out.append(line)
Path("${out_fna}").write_text("\n".join(out) + "\n")
PY
echo "[${isolate}] wrote ${out_fna}"
