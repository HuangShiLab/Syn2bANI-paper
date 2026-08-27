# SynTracker Fig. 3 independent validation workflow

This workflow reproduces the four-species comparison from
[Enav et al., Nature Biotechnology 2024](https://doi.org/10.1038/s41587-024-02276-2)
(Fig. 3) using Syn2bANI's new `anchor_adjacency` output.

The isolates are taken from the paper's Supplementary Tables 2–5 and assembled
from raw SRA reads. We then compare all within-species pairs (and, for
_H. pylori_, within-host pairs) with Syn2bANI and skani, and plot
ANI vs anchor_adjacency.

## Expected evolutionary patterns

| Species | Expected pattern |
|---|---|
| _Neisseria gonorrhoeae_ | ANI and synteny correlated (both SNPs and SVs) |
| hypermutator _E. coli_ | ANI varies widely, synteny stays high (SNP-driven) |
| _Helicobacter pylori_ | Mixed: some pairs high-ANI/high-synteny, others high-ANI/low-anchor-adjacency |
| _Streptomyces rimosus_ | ANI high/clonal, synteny varies widely (SV-driven) |

## Files produced so far (local/MacStudio)

- `data/syntracker/41587_2024_2276_MOESM3_ESM.xlsx` – downloaded SI tables
- `data/syntracker/samples_*.tsv` – isolate metadata
- `data/syntracker/references.tsv` – reference genome accessions
- `data/syntracker/fastq_manifest_*.tsv` – ENA FASTQ URLs for all 132 isolates

## HPC execution steps (on hpc2021)

All commands assume the working root `/lustre1/g/aos_shihuang/data/syntracker_validation`.

1. **Copy metadata and scripts to HPC**
   Run on the login node:
   ```bash
   bash /lustre1/g/aos_shihuang/Syn2bANI-paper/scripts/syntracker_validation/00_setup.sh
   ```

2. **Download reference genomes** (light I/O task, can run on login or I/O node)
   ```bash
   bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/01_download_references.sh
   ```

3. **Download raw FASTQ reads** (network-bound, run on an I/O node; may take hours)
   ```bash
   ssh shihuang@hpc2021-io1.hku.hk
   bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/02_download_reads.sh
   ```

4. **Assemble isolates** (compute-intensive, submit via SLURM)
   On the login node:
   ```bash
   bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/03_assemble_array.sh
   ```
   Monitor with `squeue -u shihuang`.

5. **Run Syn2bANI comparisons** (compute, submit via SLURM)
   ```bash
   bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/04_submit_syn2bani.sh
   ```

6. **Run skani comparisons** (compute, submit via SLURM)
   ```bash
   bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/05_submit_skani.sh
   ```

7. **Plot results** (can run locally after copying results)
   ```bash
   python3 scripts/syntracker_validation/06_plot_results.py \
     --syn2bani-dir data/syntracker_validation/syn2bani \
     --skani-dir data/syntracker_validation/skani \
     --metadata-dir data/syntracker_validation/samples \
     --outdir figures/syntracker_validation
   ```

## Notes / requirements

- `shovill` is tried first for assembly; if it is not installed on the compute
  nodes the script falls back to `spades.py`.
- `datasets` (NCBI CLI) must be available for reference download.
- `syn2bani` is expected at `/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani`.
- `skani` must be in `$PATH` on the compute nodes.
