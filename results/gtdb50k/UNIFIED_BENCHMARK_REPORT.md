# Unified GTDB-R207 benchmark report

Combined 43,334 held-out same-genus pairs (80-95% plus 95-100 split) with 727 high-ANI test pairs (95-97 / 97-100). Truth: dnadiff/ANIm.

## Overall (all pairs with a method estimate)

| method | n | MAE | median | bias | SD | r |
|---|---|---|---|---|---|---|
| syn2bani_raw | 43955 | 1.6474 | 1.2880 | +1.4152 | 1.6772 | 0.8672 |
| syn2bani_v6 | 40629 | 0.6180 | 0.4953 | -0.1274 | 0.8158 | 0.9672 |
| syn2bani_hybrid | 40629 | 0.6146 | 0.4908 | -0.1154 | 0.8269 | 0.9668 |
| skani | 44059 | 1.0181 | 0.7800 | -0.8686 | 0.9845 | 0.9778 |
| fastani | 44060 | 1.0321 | 0.8759 | -0.9699 | 0.8784 | 0.9834 |

## By ANI band

| band | method | n | MAE | bias | r |
|---|---|---|---|---|---|
| 80-85 | syn2bani_raw | 12094 | 1.9642 | +1.6328 | 0.1201 |
| 80-85 | syn2bani_v6 | 8823 | 0.7155 | -0.2092 | 0.4136 |
| 80-85 | syn2bani_hybrid | 8823 | 0.7193 | -0.2046 | 0.3914 |
| 80-85 | skani | 12172 | 1.9690 | -1.9674 | 0.6825 |
| 80-85 | fastani | 12172 | 1.8350 | -1.8344 | 0.8275 |
| 85-90 | syn2bani_raw | 15973 | 1.8918 | +1.6613 | 0.4917 |
| 85-90 | syn2bani_v6 | 15918 | 0.6445 | -0.0954 | 0.7741 |
| 85-90 | syn2bani_hybrid | 15918 | 0.6450 | -0.0950 | 0.7728 |
| 85-90 | skani | 16000 | 0.8930 | -0.7421 | 0.8312 |
| 85-90 | fastani | 16000 | 1.0048 | -0.9802 | 0.8762 |
| 90-95 | syn2bani_raw | 14758 | 1.2141 | +1.0478 | 0.7686 |
| 90-95 | syn2bani_v6 | 14758 | 0.5397 | -0.1047 | 0.8997 |
| 90-95 | syn2bani_hybrid | 14758 | 0.5411 | -0.1032 | 0.8987 |
| 90-95 | skani | 14758 | 0.4307 | -0.1663 | 0.9329 |
| 90-95 | fastani | 14758 | 0.4553 | -0.3042 | 0.9311 |
| 95-97 | syn2bani_raw | 496 | 0.7563 | +0.6655 | 0.6188 |
| 95-97 | syn2bani_v6 | 496 | 0.4487 | +0.0615 | 0.6932 |
| 95-97 | syn2bani_hybrid | 496 | 0.4514 | +0.0815 | 0.7092 |
| 95-97 | skani | 496 | 0.3178 | +0.0272 | 0.7919 |
| 95-97 | fastani | 496 | 0.3512 | -0.1841 | 0.7580 |
| 97-100 | syn2bani_raw | 634 | 0.2292 | +0.2022 | 0.9576 |
| 97-100 | syn2bani_v6 | 634 | 0.5474 | -0.4662 | 0.7729 |
| 97-100 | syn2bani_hybrid | 634 | 0.2324 | +0.1748 | 0.9433 |
| 97-100 | skani | 633 | 0.1386 | -0.0109 | 0.9697 |
| 97-100 | fastani | 634 | 0.2659 | -0.2220 | 0.7653 |

## Notes

- syn2bani_v6 is the ridge calibration trained on 2,520 v5 pairs + 1,614 high-ANI train pairs.

- syn2bani_raw is the gated MLE without calibration.

- syn2bani_hybrid uses ani_cal for ani_gated < 98.0% and ani_gated otherwise.

- skani/FastANI are shown where they returned a value.
