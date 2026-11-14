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

## Reviewer-round clean checkout

After the 2026-11-14 citation-order repair, a fresh GitHub clone at manuscript
commit `3e3716fe4b8526d265d7e632f25e96db83420e42` ran
`python3 scripts/audit_release.py`: **RELEASE READY, 48 checks, 0 failures**.

## Clean-checkout log

```text
git clone --depth 1 https://github.com/HuangShiLab/Syn2bANI-paper.git
HEAD: 957aadeef4f125213bbebaedc199514751785b7a
python3 scripts/audit_release.py
  verdict: RELEASE READY
  checks: 42
  failures: 0
python3 -m venv .venv
.venv/bin/python -m pip install pandas==3.0.6 numpy==2.5.4 matplotlib==3.11.2
.venv/bin/python scripts/generate_fig2_simulations.py          # pass
.venv/bin/python scripts/generate_fig_s6_simulation_mosaic.py # pass
.venv/bin/python scripts/generate_fig_s10_saureus_breakpoints.py
  n=2211, zero=1337, max=7, ANI=93.6-100.0
```

Rendered clean-checkout Figure 2 labels were rechecked as Syn2bANI 0.073,
skani 0.477 and FastANI 0.742; S10 showed the 60.5% zero-breakpoint result and
no obsolete 0.1.0 label. Regenerated PDF binaries may differ only in embedded
timestamps. The binary supervisor verdict is **RELEASE READY for
submission**, with acceptance-stage Zenodo archives explicitly stated as
release conditions.
