# skani 0.3.2 rerun for remaining small comparisons

Date: 2026-09-22.

This closes the two remaining skani 0.1.0 caveats in the Nature Methods
manuscript: the 12-pair simulation ladder and the 13-vs-1 complete
Enterobacteriaceae comparison. Both reruns used skani 0.3.2 with default
settings, matching the prior comparison policy. No low-ANI `-c` adjustment was
applied because the manuscript states that default parameters were used.

## Simulation ladder

- Input: 12 deterministic genomes in `../Syn2bANI/prototype/simindel/`.
- Command: `skani dist --ql <queries> --rl <reference> -t 8`.
- Raw output: `simindel_skani032.tsv`.
- Joined Syn2bANI/skani/FastANI table:
  `simindel_cross_tool_4e.tsv` (FastANI values are the original 1.33 values).

| comparator | n | MAE | mean bias | max abs error |
|---|---:|---:|---:|---:|
| skani 0.1.0 (withdrawn) | 12 | 0.3775 | -0.3775 | 1.38 |
| skani 0.3.2 | 12 | 0.4775 | -0.4775 | 1.60 |

## Enterobacteriaceae complete genomes

- Input: 13 complete chromosomes versus *E. coli* K-12 MG1655 in
  `../Syn2bANI/prototype/realbench/`.
- Command: `skani dist --ql perf_q13.txt --rl <reference list> -t 8`.
- Raw output: `enterobacteriaceae_skani032.tsv`.
- Per-pair summary: `enterobacteriaceae_skani032_report.txt`.

On the 8 reported pairs, the Syn2bANI gamma estimate has MAE 0.119 (bias
-0.073) against skani 0.3.2, and MAE 0.207 against FastANI. The previous
skani 0.1.0 comparison gave gamma MAE 0.094; the conclusion is unchanged.
The five very distant pairs remain below Syn2bANI's detection floor and are
not extrapolated.

## Manuscript change

The Tools and versions paragraph now states skani 0.3.2 for the simulation and
small-subset comparisons and removes the obsolete pre-release caveat. Figure 2
was regenerated from the updated cross-tool table.
