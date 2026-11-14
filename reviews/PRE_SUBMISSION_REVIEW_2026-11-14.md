# Pre-submission review and revision, 2026-11-14

## Skill packet

- Supervisor-Skills repository: `https://github.com/HuangShiLab/Supervisor-Skills`
- Reviewed HEAD: `2edb32733e4042a0408824441f0e832e2476b07b`
- Skills applied: `pre-submission-reviewer`, `supervisor-review`,
  `delivery-evidence-audit`, `computational-method-release-audit`,
  `project-experience-codifier`, plus the persisted Syn2bANI project frame,
  architecture, and claim-evidence matrix.

## Scope

Audited the current Nature Methods package at pre-repair commit `d89db7bba6bc2bc957b4e169bde29fcfc2da78b6`; repaired release commit `f049d038e392f55ce72fe9c7d246869b84602f40` before this
round, including the manuscript, Supplement, rendered DOCX/PDF, main and
Extended Data displays, source tables, scripts, checksums, release scripts and
clean-checkout path. The frozen 2026-10-21 reviewer reports were not edited.

## Findings before repair

| Severity | Finding | Repair |
|---|---|---|
| MAJOR | References were not in first-use order; old references 4/5/17/18 appeared after old reference 16, and several later dataset references were out of sequence. | Rebuilt the reference list in first-use order. |
| MAJOR | Thirteen listed references were not cited; GTDB, SynTracker, O157, FDA-ARGOS, HROM and OAPGC citations were missing or delayed. | Added first-use citations; removed unused references; renumbered all surviving citations in the manuscript and Supplement. |
| MINOR | “Robust to cluster scope” used unnecessary evaluative wording in Supplement. | Changed to “retained the same direction across cluster scope.” |
| MINOR | Several Online Methods sentences exceeded 60 words, especially rescue and structural-statistic definitions. | Retained because each sentence defines an algorithmic rule whose clause order is required; no meaning was obscured. No unsupported, filler or process-language finding remained. |

No CRITICAL issue was found in this round. The two MAJOR citation issues were
release-blocking for a Nature-style submission and were repaired.

## Changes made

1. Reordered references 1–33 by first use.
2. Removed unused references; all remaining references are cited.
3. Added explicit citations for GTDB-R207, SynTracker collections, O157:H7,
   FDA-ARGOS, *B. longum*, *H. pylori*, HROM and OAPGC.
4. Renumbered repeated citations in Online Methods, figure legends and
   Supplementary Note 8.
5. Added the comparator-aligner citations (minigraph-Cactus, SyRI and
   progressiveMauve) to the Discussion boundary statement.
6. Replaced the evaluative supplement phrase “robust to cluster scope.”
7. Updated the automated release gate to verify the 33-reference architecture
   and absence of obsolete citation strings.
8. Regenerated both DOCX files and re-rendered them to PDF.

## Review result after repair

- Abstract: 173 words.
- Main text excluding Abstract and Online Methods: approximately 2,184 words.
- Main displays: Figures 1–6 and Tables 1–2.
- Extended Data: Figures 1–5 and Tables 1–2.
- References: 33 consecutive entries, all cited and ordered by first use.
- Em dashes: 0 in manuscript and Supplement.
- Host-local paths in manuscript prose: 0.
- Banned-tone findings: none above the allowed threshold.
- Automated release gate: 48 checks, 0 failures.
- Rendered DOCX pages: manuscript 28; Supplement 23.

## Supervisor-review verdict

**READY.** No unresolved CRITICAL or MAJOR issue remains. The only non-blocking
acceptance-stage action is the already stated Zenodo deposition for the code
release and complete pair-level census tables.

## Final clean checkout

After committing the reviewer repairs, a fresh GitHub clone at manuscript commit
`3e3716fe4b8526d265d7e632f25e96db83420e42` ran
`python3 scripts/audit_release.py` successfully:

```json
{"verdict":"RELEASE READY","checks":48,"failures":[]}
```

## Post-review repair cycle (2026-11-14)

The final supervisor pass found additional release-gate violations not covered
by the initial prose-only review. Repairs were applied in checklist order:

1. **GC sweep provenance.** Removed the unsupported historical GC ladder from Methods, Supplement, Extended Data Table 1 and release scripts; the result is no longer claimed.
2. **ED Fig. 2 provenance.** Added `scripts/generate_extended_data_figure2.py`, which regenerates the packaged PDF/PNG from the unchanged SynTracker source figure and records the generator in the manifest.
3. **Reviewer minor comments.** Added `R1-m1–R1-m4`, `R2-m1–R2-m7` and `R3-m1–R3-m10` to the feedback disposition ledger with evidence for each disposition.
4. **Superseded artifacts.** Removed old manuscript variants, legacy captions, archive figures and archive result tables from the current release branch; Git history retains them.
5. **Process leakage.** Removed time-sensitive test-suite wording, caption history about withdrawn figures, and unnecessary “re-run” phrasing.
6. **Companion provenance.** Methods now specifies the Syn2b digest command, `.tgt` anchor format and batch `syn2b synteny` call.
7. **Supplement renumbering.** After GC removal, supplementary figures were renumbered S1–S14, all packaged files renamed, and every main-text reference updated.

The binary release gate now performs **55 checks** and reports:

```json
{"verdict":"RELEASE READY","checks":55,"failures":[]}
```

## Final clean checkout after repairs

A fresh clone at `f049d038e392f55ce72fe9c7d246869b84602f40` passed
`python3 scripts/audit_release.py` with **55 checks and 0 failures**. Running
`scripts/generate_extended_data_figure2.py` in that clone regenerated ED Fig. 2
and the subsequent audit remained **RELEASE READY**. Only PDF creation metadata
changed; the labeled panels remained deterministic.
