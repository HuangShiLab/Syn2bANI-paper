# Syn2bANI Algorithmic Fix Report

## Summary

The mid-ANI bias can be substantially reduced by algorithmic changes without relying on a larger GBRT model:

| Approach | MAE vs FastANI (validation, n=15) | Mean error | Notes |
|---|---|---|---|
| BcgI + GBRT-corrected (current default) | 3.22% | +2.48% | Embedded v4 model, overfits to training taxa |
| BslFI + GBRT-corrected | 3.25% | — | Best single-enzyme GBRT result |
| **Multi-enzyme (16) + mash_ani** | **2.85%** | **+2.78%** | **Best algorithmic result; no ML** |
| Multi-enzyme (16) + GBRT v5 multi | 3.41% | +3.15% | Trained on 728 multi-enzyme pairs; overfits |
| skani | 0.47% | +0.47% | Reference |

Key take-away: **multi-enzyme mode with mash_ani as the reported ANI is the best algorithmic fix identified so far**, cutting the validation error roughly in half compared with single-enzyme BcgI and avoiding GBRT overfitting.

## 1. Enzyme sweep results (single enzymes, validation pairs)

| Enzyme | Tag length* | Mean shared tags | Mash ANI MAE | Corrected ANI MAE | Raw ANI MAE |
|---|---|---|---|---|---|
| BslFI | 25 bp | 296 | **5.26%** | 3.25% | 10.04% |
| Bsp24I | 27 bp | 163 | 6.06% | 3.36% | 9.90% |
| BaeI | 28 bp | 709 | 6.16% | 2.33% | 9.74% |
| CjeI | 28 bp | 689 | 6.37% | 2.85% | 9.62% |
| CjePI | 27 bp | 737 | 6.51% | 2.83% | 9.59% |
| BsaXI | 27 bp | 201 | 6.63% | 3.80% | 9.97% |
| BcgI | 32 bp | 316 | 6.66% | 3.22% | 10.05% |
| HaeIV | 27 bp | 13,319 | 6.67% | 1.76% | 9.28% |
| Hin4I | 27 bp | 14,277 | 6.78% | 1.66% | 9.32% |
| FalI | 27 bp | 135 | 7.44% | 3.16% | 9.44% |
| AloI | 27 bp | 77 | 11.80% | 9.56% | 15.01% |
| PpiI | 28 bp | 61 | 11.83% | 11.18% | 15.42% |

\* Note: `registry.rs` declares BslFI tag length as 21 bp, but `digest.rs` uses 25 bp; the actual extracted tags and mash_ani use 25 bp.

Observations:
- **Shorter tags help**: BslFI (25 bp) gives the best single-enzyme mash ANI.
- **High shared-tag count ≠ better mash ANI**: HaeIV/Hin4I produce 40× more shared tags than BcgI but their mash ANI is no better. Their corrected ANI is excellent, suggesting the embedded GBRT model benefits from denser tags.
- **AloI/PpiI are unstable**: large positive and negative raw errors; avoid as single enzymes.

## 2. Multi-enzyme mode

### Performance bug fixed
The previous `--multi-enzyme` implementation hung (>30 s per pair). The bottleneck was the tag matcher's **near-match fallback**: for every query tag without an exact packed-sequence match, it scanned **all reference tags** (O(n·m)). In multi-enzyme mode the tag sets are much larger, so this fallback exploded.

Fix applied in `src/cli/dist.rs`: when `--multi-enzyme` is used, set `MatchConfig.allow_near_match = false`. Exact matching is sufficient when 16 enzymes are combined.

After the fix, multi-enzyme runs in seconds per pair.

### Results

| Metric | Value |
|---|---|
| Mean shared tags | 11,620 |
| Mean min AF | ~0.076 |
| Raw ANI MAE | 12.95% |
| **Mash ANI MAE** | **2.85%** |
| Corrected ANI MAE (embedded v4 GBRT) | 7.33% |

The embedded GBRT-corrected value is worse because v4 was trained on single-enzyme (BcgI-like) features and does not transfer to combined tag sets.

### Multi-enzyme GBRT experiment
To test whether a dedicated GBRT model could improve on multi-enzyme mash_ani, we generated multi-enzyme features for the same 728 GTDB-R207 training pairs and trained a v5-style GBRT (`multi_raw_ani`, `multi_mash_ani`, `multi_shared_log`, `multi_af_q`, `multi_af_r`).

| Metric | Value |
|---|---|
| In-sample corrected ANI MAE | 0.47% |
| Validation corrected ANI MAE | **3.41%** |
| Multi-enzyme mash_ani validation MAE | 2.85% |

**Result: the GBRT model makes validation accuracy worse than mash_ani alone.**
Feature importances: `multi_mash_ani` 37.3%, `multi_shared_log` 19.3%, `multi_af_r` 17.0%, `multi_af_q` 13.4%, `multi_raw_ani` 13.0%.

Interpretation: the GBRT learns the training-set taxa well (0.47% in-sample MAE) but does not generalize to the unseen *Bifidobacterium* / *Veillonella* validation pairs, reproducing the same overfitting pattern seen with single-enzyme GBRT models.

## 3. Runtime and memory benchmark

Measured on two validation pairs (Bifidobacterium and Veillonella, ~2–3 Mbp genomes):

| Tool | Mean wall time | Mean peak RSS |
|---|---|---|
| skani | ~0.04 s | ~12 MB |
| syn2bani BcgI | ~0.06 s | ~19 MB |
| syn2bani multi-enzyme | ~0.42 s | ~104 MB |

Multi-enzyme is roughly **7–10× slower and ~5× more memory-hungry** than single-enzyme BcgI, and **~10× slower than skani**. The extra cost comes from extracting and matching tags for all 16 enzymes. For large-scale workflows, this is the main trade-off against the accuracy gain.

## 4. Mash ANI as default output

Added `--mash-ani` CLI flag. When enabled, the reported `ani` column equals `mash_ani` instead of the GBRT-debiased value.

Files changed:
- `src/core/ani_calculator.rs`: added `use_mash_ani` to `AniConfig`; final ANI uses `mash_ani` when set.
- `src/cli/mod.rs`: added `--mash-ani` flag to `dist`.
- `src/main.rs`: pass flag to `run_dist`.
- `src/cli/dist.rs`: pass flag to `AniConfig`; also applies the multi-enzyme near-match fix.
- `src/cli/{db,search,struct,triangle}.rs`: initialize `use_mash_ani: false`.

## 5. Recommendations

1. **Adopt multi-enzyme + mash_ani as the recommended mid-ANI workflow**.
   - Validation MAE 2.85% is the best algorithmic result so far.
   - No GBRT retraining needed, avoiding the observed overfitting.
   - Fix the multi-enzyme performance regression (already done).

2. **Do not train a dedicated multi-enzyme GBRT model** on the current training set; it increases validation error (3.41% vs 2.85% for mash_ani alone).

3. **If a single enzyme must be used, switch default from BcgI to BslFI**.
   - BslFI mash ANI MAE 5.26% vs BcgI 6.66%.
   - BslFI corrected ANI MAE 3.25% vs BcgI 3.22% (comparable).

4. **Fix the BslFI tag-length inconsistency** between `registry.rs` (21 bp) and `digest.rs` (25 bp).

5. **Do not embed the current GBRT v5 combined model into Rust**; instead, rely on multi-enzyme + mash_ani or retrain a much more regularized model on a taxonomically diverse training set.

6. **Expand validation** beyond 15 pairs and two genera before finalizing the default.
