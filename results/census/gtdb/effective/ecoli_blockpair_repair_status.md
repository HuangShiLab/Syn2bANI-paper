# *E. coli* structural block-pair repair status

The final-table diagnostic (`ecoli_blockpair_coverage.json`) found no completely
missing or over-represented block pairs, but 100 of 1,431 *E. coli* block pairs
were partially short. The shortfall was 7,571,899 of 338,793,465 expected pairs
(97.77% coverage). There were no invalid within-cluster identifiers and no
unexpected extra pairs.

The affected outputs and checkpoints were moved to
`structural_repair_quarantine_4259011`, and exactly those 100 tasks were
resubmitted through the standard census worker. Downstream jobs are:

- structural repair: `4259080`;
- repaired normalized/deduplicated structural table and QC: `4259081`;
- repaired quality-sensitivity table: `4259083`.

The pre-repair tables remain in place until the repaired QC passes. The final
GTDB manuscript numbers should be taken from the repaired outputs, not from
`results/census/gtdb/effective/` until this handoff is complete.
