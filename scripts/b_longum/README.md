# B. longum abfA functional case study

## Data location

- Genomes: `/lustre1/g/aos_shihuang/Strain2b/data/JNU_genomes/genome2023/fna`
- Work directory: `/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA`

## Step 0: Prepare

```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
bash scripts/b_longum/prepare_b_longum.sh
```

Creates `results/b_longum_abfA/manifest.tsv` and `pairs_all_vs_all.tsv`.

## Step 1: Create metadata.tsv

Create `results/b_longum_abfA/metadata.tsv` with columns:

```
accession	abfA_status	phenotype
strain001	complete	effective
strain002	deleted	ineffective
...
```

- `abfA_status`: complete / deleted / partial
- `phenotype`: effective / ineffective / unknown

## Step 2: Run Syn2bANI all-vs-all

```bash
sbatch scripts/b_longum/s1_s2b_b_longum.slurm
```

## Step 3: Merge results

```bash
python3 scripts/b_longum/merge_b_longum_results.py
```

## Step 4: Analyze

```bash
python3 scripts/b_longum/analyze_abfA_discordance.py
```

Outputs:
- `results/b_longum_abfA/ABfA_ANALYSIS_REPORT.md`
- `results/b_longum_abfA/abfA_pair_metrics.tsv`

## Step 5: Copy results back

```bash
scp -r shihuang@hpc2021.hku.hk:/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA/*.tsv results/b_longum_abfA/
scp -r shihuang@hpc2021.hku.hk:/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA/*.md results/b_longum_abfA/
```
