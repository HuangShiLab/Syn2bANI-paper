# Syn2bANI Algorithmic Bias Analysis

## 1. Alignment fraction vs ANI (GTDB-R207 100k matrix)

### ANI bins
      bin   n  mean_ani  median_skani_af  mean_skani_af  median_s2b_af_min  mean_s2b_af_min  median_s2b_shared
 (75, 80]  34    78.008            0.000          0.000              0.005            0.006              0.000
 (80, 85] 229    83.438           17.700         17.568              0.040            0.042             98.000
 (85, 90] 213    87.071           33.570         35.267              0.078            0.084            176.000
 (90, 95] 216    92.805           64.020         63.780              0.191            0.196            312.500
(95, 100]  36    95.991           77.880         70.894              0.283            0.268            423.000

Pairs with skani align frac < 50%: 56 / 728 (7.7%)
  ANI range: 76.58% - 83.35%, mean: 79.27%
Pairs with Syn2bANI min AF < 50%: 721 / 728 (99.0%)
  ANI range: 76.58% - 97.68%, mean: 87.58%

### Thresholds: fraction of pairs with AF < 50% below given ANI
 threshold  pct_skani_af_below_50  pct_s2b_af_min_below_50
        95                    8.1                     99.3
        90                   11.8                    100.0
        85                   21.3                    100.0
        80                  100.0                    100.0

## 2. BCGI training analysis

### Overall error
- train: raw MAE=15.055%, mean raw err=+1.054%; mash MAE=5.199%, mean mash err=+3.722%
- validation: raw MAE=10.051%, mean raw err=+10.051%; mash MAE=6.656%, mean mash err=+6.656%

### By ANI bin (validation)
ani_bin  n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  median_af_min  mean_af_min  mean_unmatched_ratio
  85-88 15             0.871        10.051   10.051          6.656     6.656      316.267          0.162        0.137                 0.857

### By ANI bin (training)
ani_bin   n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  median_af_min  mean_af_min  mean_unmatched_ratio
    <80  34             0.780       -33.135   49.248         -0.840    10.893       15.294          0.005        0.006                 0.990
  80-85 229             0.834         2.653   19.989          6.275     6.793      149.135          0.040        0.042                 0.951
  85-88 154             0.863         5.735   13.546          5.375     5.492      247.240          0.070        0.072                 0.922
  88-90  59             0.890         6.212    9.257          4.367     4.367      355.525          0.110        0.116                 0.877
  90-92  76             0.910        -1.605   12.765          1.283     4.978      505.618          0.154        0.152                 0.830
  92-95 140             0.938         1.925    5.956          0.931     2.292      551.607          0.215        0.221                 0.762
 95-100  36             0.960        -3.069    7.555         -0.343     1.565      760.250          0.283        0.268                 0.707

### Per-genus validation summary
          genus  n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  mean_af_min
Bifidobacterium  9             0.871        10.112   10.112          7.408     7.408      486.222        0.166
    Veillonella  6             0.869         9.961    9.961          5.528     5.528       61.333        0.094

Scatter plots saved: /lustre1/g/aos_shihuang/Syn2bANI-paper/results/algorithm_analysis/bcgi_bias_diagnostics.png

### Correlations with raw error (validation)
- bcgi_shared_tags: 0.236
- bcgi_af_min: 0.210
- bcgi_unmatched_ratio: -0.195
- bcgi_mash_ani: 0.218

## 2. CJEPI training analysis

### Overall error
- train: raw MAE=8.657%, mean raw err=+8.398%; mash MAE=4.374%, mean mash err=+3.752%
- validation: raw MAE=9.587%, mean raw err=+9.587%; mash MAE=6.507%, mean mash err=+6.507%

### By ANI bin (validation)
ani_bin  n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  median_af_min  mean_af_min  mean_unmatched_ratio
  85-88 15             0.871         9.587    9.587          6.507     6.507      736.933          0.184        0.175                 0.820

### By ANI bin (training)
ani_bin   n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  median_af_min  mean_af_min  mean_unmatched_ratio
    <80  34             0.780        16.494   16.494          5.728     5.728       60.588          0.010        0.012                 0.986
  80-85 229             0.834        11.941   11.941          6.082     6.132      343.795          0.057        0.059                 0.937
  85-88 154             0.863         9.493    9.493          4.719     4.994      545.195          0.088        0.091                 0.903
  88-90  59             0.890         7.243    7.243          3.822     3.822      858.390          0.139        0.140                 0.852
  90-92  76             0.910         5.596    5.596          2.635     3.017     1108.211          0.177        0.186                 0.803
  92-95 140             0.938         3.005    4.355          0.220     2.114     1393.871          0.233        0.239                 0.744
 95-100  36             0.960         2.302    2.302         -1.090     1.810     1614.444          0.293        0.275                 0.714

### Per-genus validation summary
          genus  n  mean_fastani_ani  mean_raw_err  mae_raw  mean_mash_err  mae_mash  mean_shared  mean_af_min
Bifidobacterium  9             0.871         9.777    9.777          6.829     6.829      915.111        0.189
    Veillonella  6             0.869         9.304    9.304          6.025     6.025      469.667        0.154

Scatter plots saved: /lustre1/g/aos_shihuang/Syn2bANI-paper/results/algorithm_analysis/cjepi_bias_diagnostics.png

### Correlations with raw error (validation)
- cjepi_shared_tags: 0.592
- cjepi_af_min: 0.301
- cjepi_unmatched_ratio: -0.335
- cjepi_mash_ani: 0.431

## 3. Mash-like estimator formula check

Mash ANI computed from geometric mean of AF: `1 + log(sqrt(af_q * af_r)) / tag_len`

- bcgi: max recompute diff vs stored mash_ani = 0.000000
- cjepi: max recompute diff vs stored mash_ani = 0.280442

## 4. Unmatched tags hypothesis

Estimated total tags = shared_tags / af; unmatched = total - shared.

### BCGI (training)
- mean shared tags: 325.2
- mean estimated total tags (query): 2658.4
- mean unmatched ratio: 0.875
- correlation raw_err vs unmatched_ratio: 0.814

### CJEPI (training)
- mean shared tags: 759.5
- mean estimated total tags (query): 5502.6
- mean unmatched ratio: 0.863
- correlation raw_err vs unmatched_ratio: 0.807
