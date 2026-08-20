#!/bin/bash
#SBATCH --job-name=mv_srcprep
#SBATCH --partition=intel
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/srcprep_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/srcprep_%A_%a.err
# s0: build per-dataset source-genome catalog (concat fasta, lengths, first-id map)
#     + species clusters via skani triangle at >=95% ANI. array: 0=strain 1=marine
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"

DS=$( [ "$SLURM_ARRAY_TASK_ID" = "0" ] && echo strain || echo marine )
if [ "$DS" = "strain" ]; then SRCDIR=$STRAIN_DIR/source_genomes/short_read/source_genomes
else SRCDIR=$MARINE_DIR/genomes; fi
OUT=$WORK/refs/$DS
mkdir -p "$OUT" "$WORK/logs"

cd "$OUT"
ls "$SRCDIR"/*.fasta | sort > genome_paths.txt
{ echo -e "genome_id\tpath"; awk -F/ '{n=$NF; sub(/\.fasta$/,"",n); print n"\t"$0}' genome_paths.txt; } > genomes.tsv
N=$(wc -l < genome_paths.txt)
echo "[srcprep] $DS: $N source genomes"

# genome lengths + first record id (first record id is what syn2bani prints for a ref)
$PY "$SCRIPTS/fasta_index.py" genomes.tsv genome_lengths.tsv firstid_map.tsv

# concatenated reference for minimap2 (headers rewritten to genome_id|orig)
if [ ! -s source_all.fa ]; then
    $PY "$SCRIPTS/concat_sources.py" genomes.tsv source_all.fa
fi

# skani triangle -> species clusters (connected components at >=95% ANI)
if [ ! -s triangle.tsv ]; then
    $SKANI triangle -l genome_paths.txt -o triangle.tsv -t "$SLURM_CPUS_PER_TASK"
fi
$PY "$SCRIPTS/cluster_species.py" triangle.tsv genomes.tsv $SPECIES_ANI > species_clusters.tsv
echo "[srcprep] $DS done"
