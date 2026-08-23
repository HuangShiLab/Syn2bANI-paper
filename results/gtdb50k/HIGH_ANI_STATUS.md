# High-ANI unified benchmark status

## What is ready

- `high_ani_truth.tsv` (7,695 rows): dnadiff/ANIm truth for the 8,000 candidate
  high-ANI pairs downloaded from GTDB-R207 non-representative genomes.
- `high_ani_pairs_final.tsv` (2,348 rows): genome-level 60/40 train/test split
  of all pairs with `anim_ani >= 95`.
  - 95–97%: 224 train / 98 test
  - 97–100%: 1,393 train / 633 test
- HPC runners for the final 2,348 pairs:
  - `s7_s2b_high_ani_final.slurm` + `run_s2b_high_ani_final_slice.sh`
  - `s8_skani_high_ani_final.slurm` + `run_skani_high_ani_final_slice.sh`
  - `s9_fastani_high_ani_final.slurm` + `run_fastani_high_ani_final_slice.sh`
- `merge_high_ani_results.py`: will merge the three tool outputs into
  `high_ani_results.tsv` once they exist.
- `calibration_v6.py`: will retrain the linear calibrator on v5 training rows
  plus the new high-ANI train split, and evaluate on the 43,334 held-out pairs
  plus the high-ANI test split.

## Next steps on HPC

1. Pull the latest repository on the HPC login node:
   ```bash
   cd /lustre1/g/aos_shihuang/Syn2bANI-paper
   git pull origin main
   ```
2. Submit the three arrays:
   ```bash
   cd results/gtdb50k
   mkdir -p logs
   sbatch ../../scripts/gtdb50k/s7_s2b_high_ani_final.slurm
   sbatch ../../scripts/gtdb50k/s8_skani_high_ani_final.slurm
   sbatch ../../scripts/gtdb50k/s9_fastani_high_ani_final.slurm
   ```
3. After they finish, merge outputs:
   ```bash
   python3 scripts/gtdb50k/merge_high_ani_results.py
   ```
4. Train v6 and evaluate:
   ```bash
   python3 scripts/calibration_v6.py
   ```
5. Pull `high_ani_results.tsv`, `linear_cal_v6.json`, and
   `calibration_v6_*.tsv` back to the local repo for manuscript updates.

## Notes

- The split drops "mixed" pairs (one genome in train, the other in test) to
  guarantee a genome-level holdout. 2,157 of the 4,505 ANI>=95 pairs are mixed
  and therefore excluded.
- skani path in the new runner is `/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani`.
