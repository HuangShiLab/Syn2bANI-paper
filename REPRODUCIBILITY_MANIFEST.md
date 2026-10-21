# Reproducibility manifest, 2026-10-21

## Software

- Syn2bANI 0.1.0 deployed snapshot: commit `fe0f36c`, embedded calibration v5.
- Syn2bANI fragmentation fix: commit `c974f5f`.
- Syn2bANI 0.1.1 event-count fix: commit `40b7795`.
- Syn2bANI HPC benchmark snapshot used for one-to-all scaling: commit `0b58795`.
- Companion Syn2b census/structural snapshot: commit `64717ab`.
- Benchmark-only calibration v6 is not embedded in the deployed 0.1.0 snapshot.
- Python analyses use the package families listed in `environment.yml`. This is
  a reproducible runtime specification, not a byte-for-byte lockfile or
  container image. Methods report exact Python, SciPy and statsmodels versions
  where those versions affect inference.

## Models

```text
0071b9306e0ca2a69fb86909087af0d01dc294a9518bfd799cbc3562c21e53e3  environment.yml
fe54a4b2f72ea2ccc024d46e893fcc5f436344c8702e31074fc7ee7084259ec5  results/panel_by_band/linear_cal_v5.json
569f6596d61d9e97486934cff5055a38a9c1e81bf61f9554c8177a45850c1e36  results/panel_by_band/linear_cal_v6.json
```

## Analysis artifact checksums

```text
63175236883b0145e59ca3b439d75753f49bc0a3c9f498909b84964ab2b47d4c  scripts/revision/dependence_aware_uncertainty.py
8e3d7502785fd0f475914fe5027067e13702a1afe1493f07af30af39c974a6a8  scripts/revision/review_response_sensitivity.py
0cbfa32bbebabc81e690f59c8b0a19a03fed701d182939af8b06f2946a51e2ce  results/revision/dependence_aware_uncertainty.tsv
419ca27a444be1073e2c57c702345774d621619b378fe5c9446275398846ecfb  results/revision/af_stratified_benchmark.tsv
e888376f5b633b7c824130c403095e485d891cba0bf62bdef2d64eab3f08f37a  results/revision/gtdb_cluster_scope_sensitivity.tsv
4182290d9103552c55e94798aca694d4fb6776d15b9d925b528db41cbb3e5ec1  results/census/gtdb/gtdb_species_ani_summary.tsv
4929b45429bfc62160c6d5466e02ecb692845a6400a77055fcf9f40fba414b64  results/census/gtdb/gtdb_species_ani_sv_joint.tsv
f5fca9dfbec550be164cb8f319c33ff416ccf480b238f5b4fdbe8932c48aa80f  results/census/hrom/quality_stratified/quality_stratified_summary.tsv
34fa755e5fa34155cdccb4b6141b32bd775712e78939e297dc25854b4682406f  results/census/oapgc/site_effects/species_effects_all_pairs.tsv
```

## Result-to-artifact mapping

| Result | Software/model | Primary script(s) | Archived output |
|---|---|---|---|
| GTDB v5 CV | Syn2bANI `fe0f36c`, v5 | `scripts/calibration_v5.py` | `results/panel_by_band/linear_cal_v5.json` |
| GTDB v6 hold-out | benchmark-only v6 | `scripts/calibration_v6.py`, `scripts/generate_fig3_gtdb_benchmark.py` | `results/gtdb50k/*`, `paper/figures/main/fig3*` |
| GTDB census | companion Syn2b `64717ab`, skani 0.3.2 | `scripts/gtdb_census/*` | `results/census/gtdb/*` |
| HROM census | companion Syn2b `64717ab`, skani 0.3.2 | `scripts/hrom_census/*` | `results/census/hrom/*` |
| OAPGC census | companion Syn2b `64717ab`, skani 0.3.2 | `scripts/oapgc_census/*` | `results/census/oapgc/*` |
| Dependence-aware uncertainty | this repository | `scripts/revision/dependence_aware_uncertainty.py` | `results/revision/*` |
| Review sensitivity audits | this repository | `scripts/revision/review_response_sensitivity.py` | `results/revision/af_stratified_benchmark.tsv`, `results/revision/gtdb_cluster_scope_sensitivity.tsv` |
