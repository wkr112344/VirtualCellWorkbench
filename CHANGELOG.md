# Changelog

All notable changes to VirtualCellWorkbench are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-10

### Added
- Initial public release
- LINCS L1000 Level 5 beta cache loader (`load_cache`, `Cache` dataclass)
  - 720,216 signatures × 12,328 genes; 123 cell lines; 32,039 drugs
- E1 baseline model (ECFP4 drug encoder + cell embedding + MLP head)
  - 2048-d Morgan fingerprints → 1024-d → concat cell_emb(32) → 1024-d → 12,328-d
- Evaluation metrics
  - `per_sample_pcc` (paper Table 1 protocol)
  - `direction_concordance` and `topk_direction` (top-5% sign agreement)
  - `bootstrap_drug_cluster_ci` (per-drug-cluster confidence interval)
- `examples/01_train_and_predict.py`: end-to-end runnable demo
  - Trains E1 for 8 epochs, reports test PCC + bootstrap CI
  - Example prediction for BRD-K02130563 (Vorinostat) × A375
- `configs/e1_baseline.yaml`: baseline config
- `tests/test_metrics.py`: 4 unit tests for metrics
- `docs/EXPERIMENTS.md`: E0-E11 experiment index
- Companion to Wei 2026 (bioRxiv DOI 10.64898/2026.08.10.743942)

### Notes
- Headline numbers (reproducible via this release):
  - **PCC = 0.341** on official 10% drug hold-out (E0, per-sample, v7_raw)
  - **PCC = 0.222** on leave-3-cell-out (E5, mean across 41 folds)
  - **PCC = -0.27** for Vemurafenib (BRAF V600E) on A375 vs A549 (E8, anti-correlated cross-cell)
  - **PCC = 0.137 ± 0.002** across 3 seeds (E4)
- Cache (~7.4 GB) and pretrained checkpoints distributed separately via cloud drive
  (see README.md "Data access" section)