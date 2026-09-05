# GTDB-R207 scale benchmark plan

## Data verification

- **Primary GTDB-R207 directory**: `/lustre1/g/aos_shihuang/data/gtdb-r207`
  - `genomes_all/`: 65,703 representative genome FASTAs
  - `metadata/`: GTDB R207 metadata (bac120 + ar53)
  - `pair_chunks/`: 26 chunks defining the GTDB50k held-out pairs
  - `matrix_v8_chunks/`: Syn2bANI `ani` results for the held-out pairs
- **Symlink directory**: `/lustre1/g/aos_shihuang/databases/GTDB/GTDBr207` points to the same data.
- **Consistency**: The manifest, pair chunks, and matrix chunks all use the same manifest-style accessions (e.g., `GCA_000007325.1`). Mapping from GTDB accession (`GB_GCA_...`/`RS_GCF_...`) to manifest accession is: remove the `GB_`/`RS_` prefix.

## What was already computed

- **GTDB50k held-out benchmark**: 45,967 pairs (25,000 low + 616 mid + 20,351 mid_high ANI bands), plus a high-ANI expansion of 727 test pairs. This is a stratified held-out subset, **not** the full all-vs-all of 65,703 representatives (~2.16 billion pairs).
- **ANI**: `g50k_s2b` job (192 slices) completed. Total CPU time ≈ 26.5 h, elapsed ≈ 3.45 h (parallel array).
- **ANIm truth**: `dnadiff` job (472 slices) completed. Total CPU time ≈ 479.5 h, elapsed ≈ 236.5 h.
- **SV orientation (`raw_inverted_fraction`)**: `g50k_bcgI_invfrac` completed. Total CPU time ≈ 8.5 h.
- **Timestamps**: no per-pair timing logs exist, but wall time and memory were recovered from SLURM `sacct`.

## Timing summary (GTDB50k scale)

| Workflow | Job | CPU time | Avg elapsed / slice | Max RSS |
|---|---|---|---|---|
| Syn2bANI `ani` | g50k_s2b | ~26.5 h | ~1.1 min | ~1.9 GB |
| dnadiff ANIm | dnadiff | ~479.5 h | ~30 min | ~2.0 GB |
| BcgI inverted fraction | g50k_bcgI_invfrac | ~8.5 h | ~32 min | — |

- Syn2bANI ANI is ~18× faster than dnadiff in CPU time on this dataset.
- BcgI inverted-fraction SV calling is ~56× faster than dnadiff in CPU time.

## Representative query genomes (one-to-all)

Selected GTDB-R207 representative genomes from species that appear in later case studies. Manifest accessions map to filenames in `genomes_all/`.

| Species | Manifest accession | NCBI accession | Level | Strain / note |
|---|---|---|---|---|
| *Escherichia coli* | GCF_003697165.2 | GCA_003697165.2 | Complete Genome | DSM 30083 = ATCC 11775 |
| *Bifidobacterium longum* | GCF_000196555.1 | GCA_000196555.1 | Complete Genome | JCM 1217 |
| *Helicobacter pylori* | GCF_900478295.1 | GCA_900478295.1 | Complete Genome | NCTC 11637 = ATCC 43504 |
| *Staphylococcus aureus* | GCF_001027105.1 | GCA_001027105.1 | Complete Genome | — |
| *Pseudomonas aeruginosa* | GCF_001457615.1 | GCA_001457615.1 | Complete Genome | — |
| *Streptomyces rimosus* | GCF_008704655.1 | GCA_008704655.1 | Complete Genome | — |
| *Neisseria gonorrhoeae* | GCF_003315235.1 | GCA_003315235.1 | Scaffold | — |

All seven FASTA files exist in `genomes_all/`.

## Full-scale one-to-all benchmark

- **Reference set**: all 65,703 GTDB-R207 representative genomes.
- **Queries**: the 7 representative genomes above.
- **Modes**:
  1. `syn2b_fasta`: `syn2bani ani --ql query --rl all_refs`
  2. `syn2b_search`: `syn2bani search query s2b_sketch_db`
  3. `skani_fasta`: `skani dist --ql query --rl all_refs`
  4. `skani_sketch`: `skani dist --ql query --rl ref_sketch_list`
  5. `fastani`: `fastani --ql query --rl all_refs`
- **Infrastructure**: HPC SLURM array of 35 tasks (7 queries × 5 modes).
- **Prerequisite**: sketch databases for Syn2bANI and skani are being built once (job 3993221).

## Files created

- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/gtdb_one_to_all_queries.tsv`
- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/gtdb_r207_references.txt`
- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/build_gtdb_sketches.slurm`
- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/build_s2b_gtdb_sketches.py`
- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/bench_gtdb_one_to_all.py`
- `/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/bench_gtdb_one_to_all.slurm`

## Next steps

1. Wait for sketch-database build (job 3993221).
2. Submit `bench_gtdb_one_to_all.slurm`.
3. Collect `/usr/bin/time` and `sacct` metrics into `results/gtdb_one_to_all_scaling.tsv`.
4. Update Fig 4 with full-database one-to-all panel.
