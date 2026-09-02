# Syn2b structural channel on the SynTracker isolate cohorts (raw assemblies)

**Run:** 2026-09-02 on HPC  
**Label:** `raw` — SPAdes assemblies as produced, no reference ordering  
**Syn2b binary:** `/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b` at `c10bfa3`  
**Command:**
```bash
python3 scripts/syntracker_validation/08_syn2b_structural.py \
    --assembly-dir /lustre1/g/aos_shihuang/data/syntracker_validation/assemblies \
    --samples-dir  /lustre1/g/aos_shihuang/data/syntracker_validation/samples \
    --out-dir      /lustre1/g/aos_shihuang/data/syntracker_validation/syn2b_structural_raw \
    --syn2b        /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --label        raw --workers 16
```

## Self-comparison control

A genome compared against itself must give 0 breakpoints, SCJ 0 and
inverted_fraction 0. `observable_fraction < 1` is allowed on draft assemblies;
it estimates contig count.

| cohort | n | passed | median bp | median SCJ | median obs_frac | K_est |
|---|---:|---:|---:|---:|---:|---:|
| Streptomyces_rimosus | 20 | PASS | 0.0 | 0.0 | 0.9752 | 424 |
| Escherichia_coli_hypermutator | 23 | PASS | 0.0 | 0.0 | 0.9982 | 11 |
| Neisseria_gonorrhoeae | 12 | PASS | 0.0 | 0.0 | 0.9958 | 15 |
| Helicobacter_pylori | 77 | PASS | 0.0 | 0.0 | 0.9977 | 4 |

All cohorts pass the self-comparison control, so the cohort-level signal is not
confounded by the assembly-fragmentation floor that affected the earlier
Syn2bANI pass.

## Cohort summary

| cohort | role | n_pairs | self_floor_bp | median_bp | excess_over_floor | expectation |
|---|---:|---:|---:|---:|---:|---|
| Streptomyces_rimosus | positive control | 190 | 0.0 | 10.0 | **10.0** | structural variation dominates; signal must exceed floor |
| Escherichia_coli_hypermutator | negative control | 253 | 0.0 | 0.0 | **0.0** | SNP variation dominates; signal must stay at floor |
| Neisseria_gonorrhoeae | both modes | 66 | 0.0 | 3.0 | **3.0** | both SNP and synteny scores move together |
| Helicobacter_pylori | mixed | 2926 | 0.0 | 7.0 | **7.0** | mixed, per participant |

## Interpretation

- **Positive control passes.** *S. rimosus* shows a clear structural signal above
the self-comparison floor (median 10 breakpoints vs 0), matching the published
finding that this clonal cohort varies structurally while popANI is pinned near
1.0.
- **Negative control passes.** The hypermutator *E. coli* cohort stays at the
floor (median 0 breakpoints), matching the published finding that variation
there is SNP-driven.
- **Orientation channel (`inverted_fraction`) is not interpretable on these
draft assemblies.** Median `min(f, 1-f)` drifts to 0.42–0.48, consistent with
the GTDB measurement that fragmentation drives `min(f, 1-f)` toward 0.5. Read
the breakpoint/SCJ columns, not orientation, on drafts.
- **SCJ carries a `+(K-1)` fragmentation term.** Median K_est ranges from 4
(*H. pylori*) to 424 (*S. rimosus*), so `scj_distance` is dominated by contig
count; `breakpoints` is the fragmentation-immune metric.

## Files

- `self_control_raw.tsv` — per-cohort self-comparison control
- `syn2b_structural_pairs_raw.tsv` — one row per pair, full structural channel
- `syn2b_structural_summary_raw.tsv` — per-cohort medians and expectations
