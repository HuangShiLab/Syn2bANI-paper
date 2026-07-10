# Syn2bANI vs FastANI Head-to-Head Benchmark

## Dataset
- **Reference**: *E. coli* NZ_CP026351.1 (4,651,848 bp, complete chromosome)
- **Queries**: Derived with controlled SNP rates (0.05%–5%)
- **Fragmentation**: Exponential distribution, N50 from 500 bp to 100 kb
- **Completeness**: Sequential contig truncation, 30%–100%

## Key Results

### Divergence Accuracy

| Divergence | GT_ANI | Syn2bANI | S2b Error | FastANI | FA Error |
|-----------|--------|----------|-----------|---------|----------|
| 4.98% | 95.02% | 97.08% | 2.064% | 95.03% | 0.012% |
| 3.00% | 97.00% | 97.83% | 0.828% | 97.02% | 0.015% |
| 2.00% | 98.00% | 98.53% | 0.527% | 98.01% | 0.010% |
| 1.00% | 99.00% | 99.19% | 0.189% | 99.00% | 0.001% |
| 0.50% | 99.50% | 99.62% | 0.119% | 99.50% | 0.001% |
| 0.20% | 99.80% | 99.86% | 0.062% | 99.80% | 0.001% |
| 0.10% | 99.90% | 99.93% | 0.030% | 99.90% | 0.000% |
| 0.05% | 99.95% | 99.97% | 0.021% | 99.95% | 0.000% |

### Fragmentation Robustness (2% divergence baseline)

| N50 | Syn2bANI | FastANI | S2b Error | FA Error |
|-----|----------|---------|-----------|----------|
| 500 | 98.52% | 98.01% | 0.517% | 0.010% |
| 1,000 | 98.52% | 98.01% | 0.517% | 0.010% |
| 2,000 | 98.52% | 98.01% | 0.517% | 0.010% |
| 5,000 | 98.53% | 98.01% | 0.527% | 0.010% |
| 10,000 | 98.53% | 98.01% | 0.527% | 0.010% |
| 20,000 | 98.53% | 98.01% | 0.527% | 0.010% |
| 50,000 | 98.53% | 98.01% | 0.527% | 0.010% |
| 100,000 | 98.53% | 98.01% | 0.527% | 0.010% |

### Completeness Robustness (2% div, N50~10k)

| Completeness | Syn2bANI | FastANI | S2b Error | FA Error |
|-------------|----------|---------|-----------|----------|
| 30% | 98.57% | 98.01% | 0.569% | 0.011% |
| 50% | 98.57% | 98.02% | 0.562% | 0.012% |
| 60% | 98.55% | 98.02% | 0.541% | 0.011% |
| 80% | 98.55% | 98.02% | 0.543% | 0.012% |
| 100% | 98.53% | 98.01% | 0.527% | 0.010% |

## Interpretation

1. **FastANI (Python) is near-perfect on SNP-only data** because:
   - No structural variation = k-mers map perfectly
   - Full sequence available = all k-mers counted
   - This represents the "best case" for k-mer methods

2. **Syn2bANI has a small systematic overestimation** (~0.5% at 2% div, ~2% at 5% div):
   - Fixed-anchor tags that differ by >1-2 bp are excluded, biasing toward conserved regions
   - The debias model partially corrects this but needs refinement
   - **This is expected**: tag-based methods inherently sample a subset of the genome

3. **Both methods are equally robust to fragmentation/completeness** on this data:
   - Neither method is affected by N50 or completeness in this SNP-only scenario
   - **Real-world difference**: Syn2bANI would maintain accuracy with rearrangements/inversions, while FastANI's k-mer chaining would break

## Limitations

- **No structural variation**: Real MAGs have rearrangements; this test favors k-mer methods
- **Single enzyme**: Multi-enzyme mode may improve Syn2bANI accuracy
- **Python FastANI**: Simplified implementation; real FastANI has fragment-level alignment and regression correction

## Next Steps

1. Add structural variation (inversions, translocations) to test where Syn2bANI's fixed-anchor advantage manifests
2. Implement multi-enzyme consensus mode in Syn2bANI
3. Train a proper GBRT debias model
