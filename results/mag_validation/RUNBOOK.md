# RUNBOOK: Task 1 — real-MAG validation pipeline (CAMI2 strain-madness + marine)

Status: **SUBMITTED 2026-08-19** — chain live on hpc2021, stage 1 (MEGAHIT) queued/running.
Design: see `DESIGN.md` (same directory). This file supersedes the design doc for
execution details (job IDs, actual paths, quota workarounds).

## Layout

- HPC work dir: `/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/`
  (referred to as `$WORK` below)
- Pipeline scripts: `$WORK/scripts/` — mirrored from this repo's
  `scripts/mag_validation/` (local copy is authoritative; push changes with
  cat-over-ssh, scp/sftp are DISABLED on hpc2021)
- Cohort: 25 strain-madness samples (seed 42: 3,4,11,13,14,17,25,27,28,29,31,35,
  54,64,69,71,75,77,81,84,86,88,89,94,97) + all 10 marine samples
  (`$WORK/lists/sample_list.tsv`, 35 rows)

## Submitted chain (job IDs)

| Stage | Job ID | Depends on | Array | What |
|---|---|---|---|---|
| s0_sourceprep | 3918192 | — | 0-1 (strain, marine) | source catalogs, first-id maps, skani triangle species clusters |
| s0b_gtdb_prep | 3918193 | — | — | skani sketch DB of 65,703 GTDB r207 genomes + first-id map |
| s1_assemble   | 3918145 | — | 0-17 (2 samples/task) | MEGAHIT --12, min-contig-len 1000 |
| s2_depth      | (submitted by controller after s1) | s1 | 0-17 | bowtie2 --interleaved + jgi depth |
| s3_bin        | controller | s2 | 0-11 (3/task) | MetaBAT2, min contig 2500 |
| s4_checkm2    | controller | s3 | — | stage bins ≥100 kbp into bins_all/ + CheckM2 |
| s5_assign     | controller | s3,s0 | 0-11 | minimap2 asm5 contig truth (≥80% cov, ≥95% id) |
| s6_cohort     | controller | s4,s5 | — | collect/bins.tsv + pairs manifests |
| s7_repsearch  | controller | s4,s0b | 0-11 | skani search vs GTDB DB → nearest rep |
| s8_fasttools  | controller | s6,s7 | 0-11 | syn2bani dist --verbose / skani dist / fastANI per pair + bin sketches |
| s9_truth      | controller | s6 | 0-(N/15) | dnadiff ANIm per MAG-vs-anchor, chunks of 15 |
| s10_collect   | controller | s8,s9 | — | collect/*.tsv |

First controller job: 3918176 (`after_s1`, dependency afterok:3918145).

**Why controllers?** hpc2021 enforces MaxSubmit≈50 / MaxJobs≈45 per user and
array elements count against it, so the full chain (6×35+ arrays) cannot be
submitted at once. `controller.sh` jobs (1 cpu, 10 min) submit each next stage
when its dependency completes, and self-retry (`--begin=now+20min`) if the
quota is still full. All job IDs are recorded in `$WORK/jobs/jobs.tsv`.

## How to check status

```bash
ssh -o BatchMode=yes shihuang@hpc2021.hku.hk
W=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation
squeue -u shihuang -o "%.18i %.9P %.30j %.8T %.20R"
cat $W/jobs/jobs.tsv                      # every submitted stage + job id
sacct -j <JOBID> --format=JobID,State,Elapsed,ExitCode -n
ls $W/logs/                               # per-stage stdout/stderr
tail -f $W/logs/megahit_<jobid>_0.out     # e.g. stage-1 task 0
```

## Where outputs land (on HPC, all under $WORK)

- `refs/{strain,marine}/` — genome lists, lengths, first-id maps, source_all.fa,
  skani triangle, species_clusters.tsv
- `refs/gtdb/db/` — skani DB of GTDB r207 (65,703 genomes); `refs/gtdb/firstid_map.tsv`
- `asm/{ds}__{sample}/final.contigs.fa` — MEGAHIT assemblies
- `depth/{ds}__{sample}/depth.txt` — jgi contig depths (BAM deleted after use)
- `bins/{ds}__{sample}/bin.*.fa` — raw MetaBAT2 bins
- `bins_all/` — staged bins ≥100 kbp, headers rewritten `{binid}|{contig}`
  (binid = `{ds}__{sample}__bin.N`); `bins_all/bin_stats.tsv` = size/N50/ncontigs
- `checkm2/quality_report.tsv`
- `assign/{ds}__{sample}/assign.tsv` — per-bin majority genome, completeness est,
  contamination bp split (same-species vs cross-species), class
  (clean / strain-mixed / cross-species / unassigned)
- `pairs/pairs.tsv` (anchor + top-3 same-species alternates),
  `pairs/pairs_anchor.tsv` (dnadiff list)
- `repsearch/{ds}__{sample}/map.tsv` — nearest GTDB r207 rep per bin
- `fast/per_pair/{binid}.{syn2bani,skani,fastani}.tsv`, `fast/rl/` ref lists,
  `sketches/{ds}__{sample}/` — reusable syn2bani .s2ba sketches
- `truth/out/{binid}/dd.report` + `truth/rows/chunk_*.tsv`
- `collect/` — **final deliverables**: `bins.tsv`, `truth_dnadiff.tsv`
  (anim_ani, anim_af_ref, anim_af_mag, af_tier with ≥60% strict),
  `ani_fast_tools.tsv` (per bin×ref: role, syn2bani ani/ani_gated/ani_upper95/
  gate/flag/AF/anchors, skani ANI/AF, fastANI ANI/ortho-frags), `rep_map.tsv`

## Pulling results back (scp/sftp disabled — cat-over-ssh)

```bash
H="ssh -o BatchMode=yes shihuang@hpc2021.hku.hk"
W=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation
mkdir -p results/mag_validation/collect
for f in bins.tsv truth_dnadiff.tsv ani_fast_tools.tsv rep_map.tsv; do
  $H "cat $W/collect/$f" > results/mag_validation/collect/$f
done
# or everything at once:
$H "cd $W && tar czf - collect/ pairs/ bins_all/bin_stats.tsv" | tar xzf - -C results/mag_validation/
```

## Expected wall time (from DESIGN.md §6)

s0 ~30 min · s0b 1–3 h · s1 6–12 h · s2 3–6 h · s3 <1 h · s4 1–2 h ·
s5 ~2 h · s6 minutes · s7 <1 h · s8 <1 h · s9 1–2 h · s10 minutes.
Total elapsed ~2–3 days including queueing (intel partition was 84/84
allocated at submit time).

## Pre-submission smoke tests (all passed, 2026-08-19)

- s0_sourceprep COMPLETED: strain (408 genomes → 23 skani-95% clusters, heavy
  strain redundancy as expected) and marine (977 → 756 clusters) catalogs built.
- MEGAHIT 1.2.9 `--12` on real interleaved gz reads (100k pairs of strain
  sample_3): exit 0, 378 contigs, "ALL DONE".
- s2 toolchain (bowtie2-build / bowtie2 --interleaved gz / samtools sort /
  jgi_summarize_bam_contig_depths): OK on the smoke assembly.
- MetaBAT2 runs and reads the depth file correctly (smoke assembly too small to
  bin — expected).
- minimap2 `-x asm5 -c` vs strain source_all.fa: OK; **bug found & fixed** —
  minimap2 writes plain-M CIGARs without `--eqx`, so parse_assign.py now
  (a) runs minimap2 with `--eqx` and (b) falls back to NM-tag identity for
  M-only CIGARs. Re-tested on real PAF: 359/378 contigs assigned (clean),
  19 unassigned.
- syn2bani dist / skani triangle / dnadiff output formats verified live; the
  collect-stage parsers are written against the verified formats.

## Status at handoff (2026-08-19 ~22:05 HKT)

s0 done; s0b (GTDB sketch DB) and s1 (MEGAHIT, 18 tasks) PENDING on intel
(partition 84/84 allocated at submission; reasons Resources/Priority — queued
cleanly, will start as nodes free); controller 3918176 holds the rest of the
chain. No action needed unless a controller dies (see Failure recovery).

## Failure recovery (all stages are idempotent — safe to resubmit)

- Every stage skips work whose outputs already exist (assemblies, depth.txt,
  bins, quality_report.tsv, assign.tsv, per-pair fast-tool TSVs, dnadiff
  reports). Resubmitting a stage only redoes what is missing.
- If the chain stalls (controller died without submitting the next stage):
  check `$W/jobs/jobs.tsv`, then resubmit the controller manually, e.g.
  `sbatch $W/scripts/controller.sh after_s3` (stage names: after_s1, after_s2,
  after_s3, after_s45, after_s67, after_s89). Controllers skip stages whose
  name already has a job ID in jobs.tsv; delete a line to force resubmission.
- If a stage's jobs fail: fix, then resubmit just that stage with the same
  dependency, e.g. `sbatch --dependency=afterok:<s3id> $W/scripts/s5_assign.slurm`,
  then `sbatch $W/scripts/controller.sh <next_stage>` to resume the chain.
- s9_truth array size is computed by the controller from the actual pair count
  (15 pairs/task), so any cohort size works.

## Deviations from DESIGN.md

1. **HPC work dir** is `/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/`
   (design doc said `/lustre1/g/aos_shihuang/mag_validation/`).
2. **CheckM2 DB** actual location: `/lustre1/g/aos_shihuang/databases/CheckM2_db/uniref100.KO.1.dmnd`
   (design doc said `/lustre1/g/aos_shihuang/db/CheckM2_db`).
3. **Marine reads layout**: only sample_0 has the dated subdir with plain
   `anonymous_reads.fq`; samples 1–9 are `marmgCAMI2_sample_N_reads/reads/anonymous_reads.fq.gz`.
   `common.sh:reads_path()` probes both.
4. **Sample packing**: arrays are 2 or 3 samples per task (18/12-task arrays)
   because of the per-user MaxSubmit/MaxJobs quota; total compute unchanged.
5. **Species clustering for source genomes** uses skani triangle at 95% ANI
   (self-contained; handles the unpublished MIKI-NS* / N_S*_contigs strains
   that have no GTDB taxonomy).
6. **syn2bani/skani/fastANI pair naming**: syn2bani labels sequences by first
   fasta record ID, so bin headers are rewritten to `{binid}|{contig}` at
   staging and ref first-record-id maps are built in s0/s0b for the collect join.
7. MetaBAT2 only (per user decision; no MaxBin2/CONCOCT ensemble).
8. MAG-vs-GTDB-rep axis included as cross-tool-agreement-only (s7/s8).
9. **Site env collision**: hpc2021 pre-sets `SCRIPTS` in the job environment even
   under `--export=NONE`, which broke `$SCRIPTS` resolution (two s0 failures).
   Fixed: all sbatch calls use `--export=NONE` and common.sh overrides use
   `MV_*`-prefixed names. Also: inside a batch job `$0` is the slurm spool path,
   so scripts source `common.sh` via the literal `SCRIPTS_HOME` path.

## CAMI2 strain-madness gold standard (cross-check download)

`strmgCAMI2_setup.tar.gz` (491 MB) from
`https://frl.publisso.de/data/frl:6425521/strain/short_read/` was downloaded to
`/lustre1/g/aos_shihuang/data/cami2_strain/gold_standard/` (login node, curl).
Per-sample GSA (`*_contigs.tar.gz`, ~200 MB each) and BAM read mappings
(~1.9 GB each) exist at the same URL if a deeper cross-check is needed later.
The pipeline does not depend on these files.

## Scale-up to all 100 strain samples

Regenerate `lists/sample_list.tsv` (all 100), raise the array bounds in
s1/s2/s3/s5/s7/s8 (PACK×tasks ≥ 100+10), push scripts, and re-run
`submit_all.sh` — completed samples are skipped automatically.

## Incident log

- 2026-08-20: s1_assemble 3918145 failed instantly on all tasks
  (`common.sh: No such file or directory` from the slurm spool dir) — the
  job had been submitted BEFORE the `$0`-relative-sourcing fix landed, and
  slurm runs the spooled copy from submit time. Lesson: after editing any
  pipeline script, resubmit affected pending jobs. Recovered per Failure
  recovery: resubmitted s1 (3918567) + controller after_s1 (3918568);
  job IDs appended to jobs.tsv (`jid` takes the last matching line).
