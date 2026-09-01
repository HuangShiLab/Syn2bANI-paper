#!/bin/bash
#SBATCH --job-name=g50k_s2b_wrapper
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_wrapper_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_wrapper_%j.err
# s14: single-job wrapper that runs the 190 s2b slices in parallel using
#      local workers, avoiding the account array-job submission limit.
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export GTDB50K_WORK=$WORK NSLICES=190
mkdir -p "$WORK/logs"

WORKERS=8
PYTHON=/group/aos_shihuang/conda/bin/python3

$PYTHON - <<PY
import os, subprocess, sys
from multiprocessing import Pool, cpu_count

work = os.environ["GTDB50K_WORK"]
script = f"{work}/scripts/run_s2b_slice.sh"
workers = int(os.environ.get("S2B_WORKERS", "$WORKERS"))
nslices = int(os.environ.get("NSLICES", "190"))

def run_slice(sid):
    env = os.environ.copy()
    env["SLURM_ARRAY_TASK_ID"] = str(sid)
    r = subprocess.run(["bash", script, str(sid)], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        return (sid, r.returncode, r.stderr[:500])
    return (sid, 0, "")

with Pool(workers) as pool:
    for sid, rc, err in pool.imap_unordered(run_slice, range(nslices)):
        if rc != 0:
            print(f"slice {sid} failed: {err}", file=sys.stderr, flush=True)
        else:
            print(f"slice {sid} ok", flush=True)
PY
