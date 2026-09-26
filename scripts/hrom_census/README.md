# HROM within-species structural census

This is the HROM analogue of the GTDB-R207 census. It uses the same block-pair
worker and idempotent output/checkpoint format, but HROM already defines each
conspecific cluster by its representative genome, so no GTDB taxonomy inference
is needed.

## Source data

- Metadata: `/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_Conspecific-genomes-metadata.tsv`
- FASTA root: `/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_nonredundant_genomes`
- Cluster rule: `HROM_Genome_XXXX_N.fna` belongs to cluster `HROM_Genome_XXXX`
- Current inventory: **5,113 clusters**, **145,149 genomes**, **273.7 Gbp**
- Multi-genome clusters: **2,241**
- Census genomes: **142,277**
- Unique unordered pairs: **62,872,901**

HROM FASTA headers are contig IDs (`>HROM_Genome_0002_contig_1`), not genome
IDs. The HROM jobs therefore pass `--ensure-filename-genome-id`, which rewrites
headers to `>HROM_Genome_0002|HROM_Genome_0002_contig_1` before digestion. This
makes all contigs in one FASTA belong to one genome in the Syn2b TGT and
synteny output.

## Resource estimate

Measured/derived planning values are in
`/lustre1/g/aos_shihuang/databases/HROM/census/resource_estimate.tsv` after
running preparation. The pre-run estimate is:

| Item | Estimate |
|---|---:|
| Input FASTA (already present) | ~274 GB sequence / ~290–310 GB on disk |
| TGT store | ~30–45 GB |
| Structural task output (gzip) | ~4–7 GB |
| skani ANI output (gzip) | ~2–4 GB |
| Recommended shared-workdir ceiling | 80 GB hard, 60 GB soft |
| Local scratch per task | <2 GB at block size 250 |
| Structural compute | 785 core-h from the generated 250-genome block-pair plan; budget **800–1,200 core-h** with I/O variability |
| Digest compute | ~1.7 core-h raw; budget **50–150 core-h** with Lustre I/O |
| ANI compute | budget **20–80 core-h** |
| Recommended production run | 3 jobs × 16 cores × 24 h = 1,152 core-h |
| Expected wall time after digest | ~8–20 h if 3 nodes start promptly |

The HROM plan has 3,525 block-pair tasks. Although the unique-pair count is
~11× smaller than the GTDB-R207 census (62.9M vs 711M), block-pair overhead and
the top-heavy cluster-size distribution make the generated task estimate 785
core-h.

## Commands

Run from the repository root on HPC:

```bash
# 1. Build the HROM task plan and workdir.
python3 scripts/hrom_census/prepare_hrom_census.py \
    --metadata /lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_Conspecific-genomes-metadata.tsv \
    --genome-root /lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_nonredundant_genomes \
    --workdir /lustre1/g/aos_shihuang/databases/HROM/census

# Optional Lustre-heavy check:
#   add --verify-paths
# Smoke test:
#   add --limit-clusters 20 --workdir /tmp/hrom_census_smoke

# 2. Pre-digest all TGTs (idempotent).
sbatch scripts/hrom_census/hrom_digest_all.slurm

# 3. Run structural census and ANI pass in parallel.
sbatch scripts/hrom_census/hrom_census_long_job.slurm
sbatch scripts/hrom_census/hrom_ani_job.slurm

# 4. When worker drain reports zero unfinished tasks:
sbatch scripts/hrom_census/hrom_merge_job.slurm
```

## Important implementation note

`census_worker.py` now has `--ensure-filename-genome-id`. Without it, HROM TGTs
would be keyed by the first contig ID (`HROM_Genome_0002_contig_1`) rather than
the genome accession, so task row filters would fail or produce wrong genome
labels. This was verified with a two-genome Syn2b digest/synteny smoke test.
