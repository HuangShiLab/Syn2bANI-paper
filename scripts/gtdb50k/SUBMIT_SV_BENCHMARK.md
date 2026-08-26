# HPC SV benchmark submission guide

## Prerequisites

- GTDB-R207 genomes at `/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all`
- Tools at `/lustre1/g/aos_shihuang/tools` (minimap2, python3)
- dnadiff 1-to-1 results already at `results/gtdb50k/out/*/dd.1coords`
- Syn2bANI results at `results/gtdb50k/s2b_50k.tsv`

## Step 1: Run minimap2 on all 43,334 pairs

```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
sbatch scripts/gtdb50k/s11_minimap2.slurm
```

- 190 slices, ~230 pairs per slice.
- Each slice runs minimap2, parses SV metrics on the fly, and deletes PAF immediately.
- Output: `results/gtdb50k/minimap2_rows/slice_{0..189}.tsv`

## Step 2: Re-parse dnadiff with large-gap filters

```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
python3 scripts/gtdb50k/parse_dnadiff_sv_filtered.py 5000 results/gtdb50k/sv_truth_50k_min5000.tsv
python3 scripts/gtdb50k/parse_dnadiff_sv_filtered.py 10000 results/gtdb50k/sv_truth_50k_min10000.tsv
```

## Step 3: Merge and compare

```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
python3 scripts/gtdb50k/analyze_sv_comparison.py
```

Outputs:
- `results/gtdb50k/sv_comparison_merged.tsv`
- `results/gtdb50k/SV_COMPARISON_REPORT.md`

## Step 4: Copy results back to local

```bash
scp -r shihuang@hpc2021.hku.hk:/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/SV_COMPARISON_REPORT.md results/gtdb50k/
scp -r shihuang@hpc2021.hku.hk:/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/sv_comparison_merged.tsv results/gtdb50k/
```
