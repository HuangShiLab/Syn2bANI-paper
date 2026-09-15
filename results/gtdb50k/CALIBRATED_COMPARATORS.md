# Comparators given the same calibration protocol

Band-holdout linear recalibration of skani and FastANI on their own outputs, evaluated on the GTDB-R207 held-out set against ANIm. Produced by `analysis/calibrated_comparators.py`.

## MAE (ANI points) by band, common subset (39,903 pairs with a Syn2bANI calibrated value)

| method                                             |   80-85 |   85-90 |   90-95 |   95-100 |   all |
|:---------------------------------------------------|--------:|--------:|--------:|---------:|------:|
| FastANI                                            |   1.828 |   1.004 |   0.455 |    0.361 | 0.977 |
| FastANI + linear(ani), band-holdout                |   0.456 |   0.531 |   0.621 |    0.545 | 0.548 |
| FastANI + linear(ani, mapped, total), band-holdout |   0.439 |   0.472 |   0.531 |    0.529 | 0.487 |
| skani                                              |   1.984 |   0.893 |   0.431 |    0.319 | 0.958 |
| skani + linear(ani), band-holdout                  |   0.828 |   0.652 |   0.823 |    0.561 | 0.753 |
| skani + linear(ani, AF), band-holdout              |   0.832 |   0.652 |   0.861 |    0.567 | 0.768 |
| syn2bani calibrated v5                             |   0.716 |   0.645 |   0.54  |    0.441 | 0.619 |
| syn2bani raw gated                                 |   2.201 |   1.894 |   1.214 |    0.754 | 1.699 |

## MAE by band, all 43,334 pairs (Syn2bANI calibrated is undefined on 3,431 of them)

| method                                             |   80-85 |   85-90 |   90-95 |   95-100 |   all |
|:---------------------------------------------------|--------:|--------:|--------:|---------:|------:|
| FastANI                                            |   1.835 |   1.005 |   0.455 |    0.361 | 1.045 |
| FastANI + linear(ani), band-holdout                |   0.444 |   0.532 |   0.621 |    0.545 | 0.538 |
| FastANI + linear(ani, mapped, total), band-holdout |   0.434 |   0.472 |   0.531 |    0.529 | 0.482 |
| skani                                              |   1.969 |   0.893 |   0.431 |    0.319 | 1.032 |
| skani + linear(ani), band-holdout                  |   0.803 |   0.654 |   0.823 |    0.561 | 0.753 |
| skani + linear(ani, AF), band-holdout              |   0.813 |   0.654 |   0.861 |    0.567 | 0.769 |
| syn2bani calibrated v5                             |   0.716 |   0.645 |   0.54  |    0.441 | 0.619 |
| syn2bani raw gated                                 |   1.964 |   1.892 |   1.214 |    0.754 | 1.67  |

## Pearson r by band, common subset

| method                                             |   80-85 |   85-90 |   90-95 |   95-100 |   all |
|:---------------------------------------------------|--------:|--------:|--------:|---------:|------:|
| FastANI                                            |   0.809 |   0.875 |   0.931 |    0.678 | 0.98  |
| FastANI + linear(ani), band-holdout                |   0.809 |   0.875 |   0.931 |    0.678 | 0.974 |
| FastANI + linear(ani, mapped, total), band-holdout |   0.815 |   0.894 |   0.943 |    0.608 | 0.979 |
| skani                                              |   0.657 |   0.831 |   0.933 |    0.727 | 0.973 |
| skani + linear(ani), band-holdout                  |   0.657 |   0.831 |   0.933 |    0.727 | 0.955 |
| skani + linear(ani, AF), band-holdout              |   0.669 |   0.834 |   0.932 |    0.624 | 0.954 |
| syn2bani calibrated v5                             |   0.414 |   0.774 |   0.9   |    0.607 | 0.962 |
| syn2bani raw gated                                 |  -0.061 |   0.488 |   0.769 |    0.512 | 0.834 |

## The 3,431 pairs without a calibrated Syn2bANI value

Bands: {'80-85': 3349, '85-90': 82}
Syn2bANI raw gated MAE on them: 1.326
skani MAE on them: 1.903; skani + linear: 0.747
FastANI MAE on them: 1.838; FastANI + linear: 0.417

