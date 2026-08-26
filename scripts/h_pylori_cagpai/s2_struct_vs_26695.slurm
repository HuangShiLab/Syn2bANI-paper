#!/bin/bash
# s2_struct_vs_26695.slurm — struct every cohort genome vs 26695 (cagPAI+)
# and vs the empty-site construct, to classify cagPAI state per strain.
#SBATCH --job-name=hp_struct
#SBATCH --output=hp_struct_%A_%a.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-600%50

# Env expected (edit to site):
#   GENOMES   dir with *.fna cohort genomes
#   REFDIR    dir with hp26695.fna (cagPAI+) and hp26695_delpai.fna (empty site)
#   OUT       output dir for per-genome struct TSVs
#   S2B       path to syn2bani binary (>= 4d679cf)
GENOMES=${GENOMES:-$PWD/genomes}
REFDIR=${REFDIR:-$PWD/refs}
OUT=${OUT:-$PWD/struct_vs_ref}
S2B=${S2B:-$TOOLS/syn2bani/syn2bani}
mkdir -p "$OUT"

f=$(ls "$GENOMES"/*.fna | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "$f" ] && exit 0
id=$(basename "$f" .fna)

for ref in hp26695 hp26695_delpai; do
    o="$OUT/${id}.vs_${ref}.struct.tsv"
    [ -s "$o" ] && continue
    # struct reports events from the query's perspective; run the cohort
    # genome as query against each reference.
    "$S2B" struct "$f" "$REFDIR/${ref}.fna" -o "$o" 2>/dev/null
done
echo "done $id"
