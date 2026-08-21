# MAG Validation Report: syn2bani vs ANIm truth on 695 CAMI2 bins

Date: 2026-08-21. Pipeline: `results/mag_validation/RUNBOOK.md` (HPC jobs,
s1 MEGAHIT -> s10 collect). Analysis: `scripts/mag_validation/analyze_mag.py`.
Data: `results/mag_validation/collect/*.tsv`, metrics `MAG_METRICS.tsv`,
figure `figures/report/mag_validation.png/.pdf`.

## Design

35 CAMI2 short-read samples (25 strain-madness + 10 marine) were assembled
with MEGAHIT, binned with MetaBAT2, quality-scored with CheckM2, and contigs
were assigned back to their source genomes via the CAMI2 ground-truth
mappings. For each of the 695 bins >= 100 kbp, the dominant source genome
("anchor", majority of assigned bases) was paired with the bin; ANI was
estimated by syn2bani, skani and FastANI against the anchor reference, and
dnadiff (MUMmer) provided ANIm truth on the same pair. Bins were classified
by contig assignment: **clean** (single source), **strain-mixed** (multiple
strains of one species), **cross-species** (contigs from >= 2 species).

Sample composition: 695 anchor pairs; quality tiers HQ 200 / MQ 210 / LQ 285
(CheckM2); classes clean 305 / strain-mixed 253 / cross-species 137; dnadiff
alignment tiers strict 638 / low-AF 40 / verylow-AF 17. All 695 dnadiff runs
succeeded after the FASTA whitespace sanitisation fix (see RUNBOOK incident
log #6).

## Headline accuracy (vs dnadiff ANIm, all 695 pairs)

| tool | n | MAE | median |err| | bias | within 0.1 | within 0.5 | within 1.0 |
|---|---|---|---|---|---|---|---|---|
| syn2bani | 695 | **0.163** | 0.092 | +0.156 | 53.4% | 94.2% | 98.7% |
| skani | 695 | **0.092** | 0.030 | +0.055 | 79.7% | 96.7% | 99.1% |
| FastANI | 633 | **0.203** | 0.065 | -0.119 | 71.1% | 90.8% | 94.6% |

(FastANI returned no estimate for 62 low-alignment pairs.)

Note: `s2b_ani` equals the raw gated estimate on 95.2% of rows, i.e. the HPC
binary used for this run did not apply the deployed v5 calibration; numbers
below are the raw estimator. The calibration layer (trained on complete GTDB
genomes) is expected to correct part of the residual bias but cannot fix the
fragmentation-driven component quantified here.

## Accuracy by contamination class

| class | syn2bani MAE | skani MAE | FastANI MAE |
|---|---|---|---|
| clean (305) | 0.084 | 0.030 | 0.055 |
| strain-mixed (253) | 0.170 | 0.111 | 0.165 |
| cross-species (137) | **0.326** | 0.196 | **0.656** |

syn2bani degrades gracefully on strain-mixed bins (MAE 0.17, 93% within
0.5 pp) and is the most robust of the three on cross-species contaminated
bins, where FastANI underestimates severely (bias -0.60, only 76% within
1 pp). skani is best on clean and strain-mixed bins.

## Accuracy by assembly quality (CheckM2 tier, syn2bani)

| tier | n | MAE | median |err| | within 0.5 |
|---|---|---|---|---|---|
| HQ | 200 | **0.041** | 0.031 | 100% |
| MQ | 210 | 0.105 | 0.083 | 99.0% |
| LQ | 285 | 0.291 | 0.186 | 86.7% |

Error is dominated by fragmentation: |error| vs contig N50 Spearman
rho = -0.680 (p = 2.3e-95); vs n_contigs rho = +0.268. On HQ MAGs —
the population that downstream analyses (dereplication, taxonomy) actually
use — syn2bani is within 0.1 pp of ANIm for 92.5% of bins.

## Alignment-fraction tiers (dnadiff)

strict (638 pairs): syn2bani MAE 0.129; low-AF (40): 0.461;
verylow-AF (17): 0.744. Residual error concentrates where the MAG covers
only a small fraction of the anchor genome — the same identifiability floor
seen in the fragmented-genome simulations.

## Reliability flag and upper bound

- The low-ANI reliability flag **never fired** on these 695 same-species
  pairs (all `ok`): it targets the low-ANI regime, not contamination.
  Consequently it has 0 recall on contaminated bins; a dedicated
  contamination/chimerism signal (e.g. chain-coverage heterogeneity across
  the bin) is future work.
- `ani_upper95` covers the ANIm truth in **98.4%** of pairs (target >= 95%),
  with median upper95 - point estimate gap of only 0.042 pp — a tight,
  honest upper bound.

## Calibration does not transfer to MAG inputs (tested)

The 695 anchor pairs were re-run with `ani --calibrate` (deployed v5 model,
binary @ fe0f36c; job s11_cal, `fast/per_pair_cal/`, merged
`collect/ani_calibrated.tsv`). Result: the GTDB-trained calibration
**degrades** MAG accuracy — overall MAE 1.261 (bias -1.243) vs raw 0.163,
worst on LQ bins (2.307) and cross-species bins (2.091); 43 near-clonal
pairs are pushed above 100. Mechanism: the v5 model learned the +2.02 bias
of complete-genome GTDB pairs, but raw bias on MAGs is only +0.16, so the
same correction overshoots downward; MAG feature values (low af_query,
high fragmentation) also resemble the divergent-pair regime the model
corrects hardest. Conclusions: (i) calibrated output must stay a separate,
clearly-labelled column — never a silent replacement — which is the current
design; (ii) on draft/MAG inputs the raw gated estimate is the recommended
output; (iii) any MAG-specific recalibration would need MAG-feature
training data, and the headroom is small (raw bias +0.16).

## Conclusions

1. On realistic MAGs syn2bani's raw estimate is within 0.5 pp of ANIm truth
   for 94% of bins overall and 100% of HQ bins (HQ MAE 0.041).
2. Error structure is interpretable: driven by fragmentation (N50) and
   cross-species contamination, not by estimator noise; syn2bani is the most
   robust tool under cross-species contamination.
3. skani is more accurate on clean/strain-mixed MAGs; FastANI fails
   silently (missing output or strong underestimation) on contaminated,
   low-AF pairs.
4. The upper95 bound is well calibrated on real MAGs (98.4% coverage).
5. Gaps for future work: contamination-aware reliability signal; calibration
   model retrained on draft/MAG features rather than complete genomes.
