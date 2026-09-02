# SynTracker cohort re-analysis with the post-fix Syn2b structural channel

**Date:** 2026-09-02

The structural channel was re-run on the four SynTracker isolate collections with
the reference-side inflation fix (`Syn2bANI c974f5f`). Self-comparison controls now
read zero breakpoints on every cohort, so the cohort-level statistics are no
longer dominated by an assembly-fragmentation floor.

## Inputs

- Assemblies: `data/syntracker_validation/assemblies/` (132 isolates, present on HPC)
- Structural output: `data/syntracker_validation/syn2b_structural_raw/syn2b_structural_pairs_raw.tsv`
- skani ANI: `data/syntracker_validation/skani/`
- Metadata: `data/syntracker_validation/samples/`
- Analysis script: `scripts/syntracker_validation/09_analyze_structural_vs_ani.py`

```bash
python3 scripts/syntracker_validation/09_analyze_structural_vs_ani.py \
  --structural data/syntracker_validation/syn2b_structural_raw/syn2b_structural_pairs_raw.tsv \
  --skani-dir data/syntracker_validation/skani \
  --metadata-dir data/syntracker_validation/samples \
  --outdir figures/syntracker_validation
```

## Self-comparison control

| Cohort | n | median breakpoints | median SCJ | median inverted_fraction | control |
|---|---:|---:|---:|---:|:---:|
| *Streptomyces rimosus* | 20 | 0.0 | 0.0 | 0.0 | PASS |
| hypermutator *E. coli* | 23 | 0.0 | 0.0 | 0.0 | PASS |
| *Neisseria gonorrhoeae* | 12 | 0.0 | 0.0 | 0.0 | PASS |
| *Helicobacter pylori* | 77 | 0.0 | 0.0 | 0.0 | PASS |

## Cohort summary (post-fix)

| species | n_pairs | median breakpoints | excess over floor | median inverted_fraction | median SCJ | expectation |
|---|---:|---:|---:|---:|---:|:---|
| *S. rimosus* | 190 | 10.0 | 10.0 | 0.4895 | 760.0 | positive control: structural signal present |
| hypermutator *E. coli* | 253 | 0.0 | 0.0 | 0.4946 | 24.0 | negative control: structural signal at floor |
| *N. gonorrhoeae* | 66 | 3.0 | 3.0 | 0.4873 | 41.0 | both modes |
| *H. pylori* | 2926 | 7.0 | 7.0 | 0.4947 | 35.0 | mixed, participant-dependent |

The positive/negative control ordering now matches the published SynTracker
interpretation: *S. rimosus* shows excess breakpoints, while the hypermutator
*E. coli* set stays at the floor.

## ANI vs structural metrics

Spearman correlations between skani ANI and the Syn2b structural channel:

| species | rho(ANI, breakpoints) | rho(ANI, inverted_fraction) | rho(ANI, SCJ) |
|---|---:|---:|---:|
| hypermutator *E. coli* | −0.122 | −0.072 | −0.278 |
| *H. pylori* | **−0.873** | −0.091 | −0.746 |
| *N. gonorrhoeae* | −0.095 | 0.093 | −0.792 |
| *S. rimosus* | −0.455 | −0.154 | −0.637 |

The strongest signal is in *H. pylori*: lower ANI is associated with more
breakpoints, consistent with a mixed population where some pairs are
substantially diverged. The orientation channel (`inverted_fraction`) is near
0.5 across all cohorts, as expected for short-read draft assemblies where most
pairs are too fragmented for a reliable strand call.

## Outputs

- `figures/syntracker_validation/syntracker_ani_vs_breakpoints.png`
- `figures/syntracker_validation/syntracker_ani_vs_inverted_fraction.png`
- `figures/syntracker_validation/syntracker_h_pylori_ani_vs_breakpoints_by_host.png`
- `figures/syntracker_validation/correlation_summary.tsv`
- `figures/syntracker_validation/merged_ani_synteny.tsv`
- `results/syntracker_validation/syntracker_summary.tsv` (post-fix)
- `results/syntracker_validation/syntracker_summary_pre_fix.tsv` (old, kept for comparison)

The pre-fix `syntracker_summary.tsv` and the old figures derived from it should
not be cited.
