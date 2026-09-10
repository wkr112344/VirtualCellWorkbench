# E0–E11 experiment reproduction logs

This directory mirrors the 11-experiment reproduction suite described in
`docs/stage65_e0_calibration.md` and the bioRxiv preprint §3.

Each experiment is recorded as a single JSON or CSV file with the full
configuration, raw per-pair scores, and a one-paragraph interpretation.
These files are the source of truth for the paper's Table 1, Figures 1–4,
and the supplementary material.

| Stage | Experiment | Output file | Status |
|---|---|---|---|
| E0 | Reproduction fidelity vs. paper | `stage65_e0_calibration.md` | ✅ |
| E1 | ECFP4 + cell-embedding baseline | `stage67_e1_ecfp4_mlp.json` | ✅ PCC=0.137 |
| E3a | –cell ablation | `stage69_e3_ablations.json` | ✅ PCC=0.105 |
| E3b | +MoA ablation | `stage69_e3_ablations.json` | ✅ PCC=0.141 |
| E4 | 3-seed repeatability | `stage70_e4_3_seeds.json` | ✅ 0.136 ± 0.002 |
| E5 | Leave-3-cell-out (unseen cell) | `stage68_e5_unseen_cell.json` | ✅ PCC=0.222 |
| E6 | CPI strict (no target overlap) | `stage71_e6_cpi_strict.json` | ✅ +0.022 inflation |
| E7 | Effect-size stratification | (stage55, original session) | ✅ |
| E8 | 12 biological case studies | `stage73_e8_cases.json` | ✅ Vemurafenib A375/A549 cross-PCC = –0.27 |
| E9 | Per-cell-line difficulty | (stage60) | ✅ β_tumor = –0.063, p=0.036 |
| E10 | Nearest-neighbor baseline | (stage61) | ✅ PCC=0.072, 4.7× weaker |
| E11 | Drug novelty stratification | (stage62) | ✅ r = –0.965 |
