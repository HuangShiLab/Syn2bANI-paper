# GTDB-R207 50k held-out benchmark report

43,334 same-genus pairs sampled from GTDB R207 representative genomes, disjoint from the 2,541-pair v5 calibration training set (both directions excluded). Truth: nucmer/dnadiff ANIm. Methods: syn2bani raw gated estimate, syn2bani v5-calibrated estimate (true extrapolation), skani.

## Overall

| method | n | MAE | median | bias | SD | r |
|---|---|---|---|---|---|---|
| syn2bani_raw | 43229 | 1.6701 | 1.3110 | +1.4346 | 1.6839 | 0.8494 |
| syn2bani_cal | 39903 | 0.6194 | 0.4977 | -0.1223 | 0.8179 | 0.9620 |
| skani | 43334 | 1.0324 | 0.8000 | -0.8829 | 0.9859 | 0.9752 |

## By ANI band (skani pre-screen band)

| band | method | n | MAE | bias | r |
|---|---|---|---|---|---|
| 80-85 | syn2bani_raw | 12094 | 1.9642 | +1.6328 | 0.1201 |
| 80-85 | syn2bani_cal | 8823 | 0.7155 | -0.2092 | 0.4136 |
| 80-85 | skani | 12172 | 1.9690 | -1.9674 | 0.6825 |
| 85-90 | syn2bani_raw | 15973 | 1.8918 | +1.6613 | 0.4917 |
| 85-90 | syn2bani_cal | 15918 | 0.6445 | -0.0954 | 0.7741 |
| 85-90 | skani | 16000 | 0.8930 | -0.7421 | 0.8312 |
| 90-95 | syn2bani_raw | 14758 | 1.2141 | +1.0478 | 0.7686 |
| 90-95 | syn2bani_cal | 14758 | 0.5397 | -0.1047 | 0.8997 |
| 90-95 | skani | 14758 | 0.4307 | -0.1663 | 0.9329 |
| 95-100 | syn2bani_raw | 404 | 0.7544 | +0.6653 | 0.5120 |
| 95-100 | syn2bani_cal | 404 | 0.4411 | +0.0701 | 0.6065 |
| 95-100 | skani | 404 | 0.3185 | +0.0335 | 0.7266 |

## Common subset (pairs with a calibrated score; BELOW_DETECTION excluded)

| band | method | n | MAE | bias | r |
|---|---|---|---|---|---|
| all | syn2bani_raw | 39903 | 1.6988 | +1.4988 | 0.8336 |
| all | syn2bani_cal | 39903 | 0.6194 | -0.1223 | 0.9620 |
| all | skani | 39903 | 0.9575 | -0.7962 | 0.9734 |
| 80-85 | syn2bani_raw | 8823 | 2.2012 | +1.9924 | -0.0615 |
| 80-85 | syn2bani_cal | 8823 | 0.7155 | -0.2092 | 0.4136 |
| 80-85 | skani | 8823 | 1.9843 | -1.9825 | 0.6571 |
| 85-90 | syn2bani_raw | 15918 | 1.8937 | +1.6646 | 0.4875 |
| 85-90 | syn2bani_cal | 15918 | 0.6445 | -0.0954 | 0.7741 |
| 85-90 | skani | 15918 | 0.8931 | -0.7438 | 0.8310 |
| 90-95 | syn2bani_raw | 14758 | 1.2141 | +1.0478 | 0.7686 |
| 90-95 | syn2bani_cal | 14758 | 0.5397 | -0.1047 | 0.8997 |
| 90-95 | skani | 14758 | 0.4307 | -0.1663 | 0.9329 |
| 95-100 | syn2bani_raw | 404 | 0.7544 | +0.6653 | 0.5120 |
| 95-100 | syn2bani_cal | 404 | 0.4411 | +0.0701 | 0.6065 |
| 95-100 | skani | 404 | 0.3185 | +0.0335 | 0.7266 |

## By phylum (n>=200)

| phylum | method | n | MAE | bias |
|---|---|---|---|---|
| Actinobacteriota | syn2bani_raw | 11450 | 1.3956 | +1.1884 |
| Actinobacteriota | syn2bani_cal | 10619 | 0.4964 | -0.0134 |
| Actinobacteriota | skani | 11450 | 1.0303 | -0.9086 |
| Bacteroidota | syn2bani_raw | 4348 | 1.8294 | +1.4770 |
| Bacteroidota | syn2bani_cal | 4074 | 0.6616 | -0.0506 |
| Bacteroidota | skani | 4359 | 1.1706 | -0.9900 |
| Campylobacterota | syn2bani_raw | 475 | 1.7555 | +1.5403 |
| Campylobacterota | syn2bani_cal | 440 | 0.5857 | +0.1998 |
| Campylobacterota | skani | 476 | 0.6616 | -0.3466 |
| Cyanobacteria | syn2bani_raw | 1704 | 1.7571 | +1.6341 |
| Cyanobacteria | syn2bani_cal | 1696 | 0.7136 | -0.4546 |
| Cyanobacteria | skani | 1705 | 0.8543 | -0.7341 |
| Firmicutes | syn2bani_raw | 3798 | 1.3149 | +0.8813 |
| Firmicutes | syn2bani_cal | 3582 | 0.5649 | -0.2026 |
| Firmicutes | skani | 3811 | 1.0581 | -0.9293 |
| Firmicutes_A | syn2bani_raw | 2520 | 2.1806 | +1.9844 |
| Firmicutes_A | syn2bani_cal | 2505 | 0.6365 | +0.1743 |
| Firmicutes_A | skani | 2521 | 0.8971 | -0.5976 |
| Firmicutes_C | syn2bani_raw | 290 | 2.0474 | +1.7766 |
| Firmicutes_C | syn2bani_cal | 288 | 0.6061 | +0.2821 |
| Firmicutes_C | skani | 290 | 0.8142 | -0.4233 |
| Halobacteriota | syn2bani_raw | 600 | 1.6222 | +1.2997 |
| Halobacteriota | syn2bani_cal | 584 | 0.5604 | -0.2060 |
| Halobacteriota | skani | 600 | 1.1041 | -0.8630 |
| Proteobacteria | syn2bani_raw | 14506 | 1.6431 | +1.4408 |
| Proteobacteria | syn2bani_cal | 12786 | 0.6952 | -0.3287 |
| Proteobacteria | skani | 14555 | 1.0852 | -0.9930 |
| Spirochaetota | syn2bani_raw | 382 | 1.7757 | +1.4554 |
| Spirochaetota | syn2bani_cal | 303 | 0.6679 | +0.4161 |
| Spirochaetota | skani | 382 | 0.7563 | -0.4303 |
| Thermoplasmatota | syn2bani_raw | 228 | 2.8279 | +2.7964 |
| Thermoplasmatota | syn2bani_cal | 220 | 0.6626 | +0.1572 |
| Thermoplasmatota | skani | 231 | 0.6923 | -0.4641 |
| Thermoproteota | syn2bani_raw | 284 | 3.4897 | +3.4327 |
| Thermoproteota | syn2bani_cal | 276 | 0.7556 | +0.0062 |
| Thermoproteota | skani | 304 | 0.7974 | -0.4861 |
| Verrucomicrobiota | syn2bani_raw | 358 | 2.6635 | +2.6312 |
| Verrucomicrobiota | syn2bani_cal | 357 | 0.6095 | +0.4043 |
| Verrucomicrobiota | skani | 358 | 0.8057 | -0.5854 |

## Diagnostics

- ani_cal > 100: 0 pairs (0.000%), by band: {}
- flag distribution: {'ok': 34761, 'INCONSISTENT': 5142, 'BELOW_DETECTION': 3431}
- gate distribution: {'gamma': 36016, 'uniform': 4213, 'uniform_fallback': 3000, 'none': 105}
- upper95 coverage (anim_ani <= ani_upper95): 99.39%
- AF tiers of truth set: {'strict': 22636, 'low-AF': 20639, 'verylow-AF': 59}
