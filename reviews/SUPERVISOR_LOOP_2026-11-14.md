# Supervisor skill loop, 2026-11-14

## Frozen goal and completion rule

Deliver the current Syn2bANI Nature Methods manuscript only when:

1. every claim has a mapped current source and display;
2. all reviewer feedback has a fixed/answered/rejected disposition;
3. every captioned panel exists in the packaged rendering;
4. current figure generators run without private-path dependencies;
5. `scripts/audit_release.py` returns `RELEASE READY`;
6. DOCX/PDF package is regenerated and inspected.

## Skill packet

- Repository: `https://github.com/HuangShiLab/Supervisor-Skills`
- HEAD reviewed: `d3cfcf72a96237f907f884bb534f7a4ad004fe47`
- Skills used: `supervisor`, `supervisor-review`, `delivery-evidence-audit`,
  `computational-method-release-audit`, and the local generalized counterparts
  `benchmark-design-guard`, `paper-architecture`, and `claim-evidence-check`.

## Work ledger

| ID | Task | Acceptance | Status |
|---|---|---|---|
| T01 | Update/install current Supervisor-Skills | latest remote HEAD copied and descriptions verified | done |
| T02 | Scope and frame gate | `PROJECT_FRAME.md`, `PAPER_ARCHITECTURE.md`, `CLAIM_EVIDENCE_MATRIX.md` current | done |
| T03 | Feedback disposition gate | `reviews/FEEDBACK_DISPOSITION_LEDGER_2026-11-14.md`, no unresolved actionable row | done |
| T04 | Claim-evidence and denominator gate | recalculate simulation, S. aureus, OAPGC and GTDB source values | done |
| T05 | Provenance gate | manifest checksums and simulation manifests pass; no host paths | done |
| T06 | Rendered package gate | Fig. 2 labels match; S10 shows v0.1.1 result; main/ED numbering correct | done |
| T07 | Clean-checkout gate | fresh GitHub clone runs `scripts/audit_release.py` | done |
| T08 | Method-release audit | event-level, biological, aggregate and runtime claims separated; archive plan explicit | done |

## Round-2 findings and repairs

| Gate | Finding | Repair |
|---|---|---|
| PK04/RP | `simindel/manifest.tsv` contained `/Users/...` paths | replaced with source-repository-relative paths and updated checksum |
| EV04 | Text said 20–201 contigs while source/figure levels are 20–200 | corrected Results, Methods, caption, Table 3 and release docs |
| EV01/MS01 | Abstract calibration sentence could imply all complete-genome calibration was unavailable | changed to “Benchmark-only v6 calibration is not embedded” |
| CM09 | Code release lacked an explicit durable-archive condition | added tagged release plus Zenodo DOI before acceptance |
| EV01/EV04 | Fig. 5 did not distinguish 185 isolates from 183 locus-analyzed non-reference genomes; title said “loss” beyond curated labels | caption now gives 127/56 denominator and uses “abfA-cluster difference” |

## Final gate result

`scripts/audit_release.py` returned:

```json
{"verdict":"RELEASE READY","checks":42,"failures":[]}
```

A fresh clean clone check is recorded below in the loop log. The binary
supervisor verdict is **RELEASE READY for submission**, with acceptance-stage
Zenodo archives explicitly stated as release conditions.
