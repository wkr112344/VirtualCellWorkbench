# VirtualCellWorkbench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22694467.svg)](https://doi.org/10.5281/zenodo.22694467)
[![GitHub release](https://img.shields.io/badge/release-v0.1.0-blue)](https://github.com/wkr112344/VirtualCellWorkbench/releases/tag/v0.1.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> A PyTorch toolkit for predicting chemical-perturbation transcriptional responses across cancer and normal cell lines.
> Companion code for Wei (2026), *Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck* (DOI: [10.64898/2026.08.10.743942](https://www.biorxiv.org/content/10.64898/2026.08.10.743942v1)).

## Citing

If you use VirtualCellWorkbench in your research, please cite both the **preprint** and the **software**:

> Wei, K. (2026). *VirtualCellWorkbench v0.1.0* (Version 0.1.0). Zenodo. https://doi.org/10.5281/zenodo.22694467

> Wei, K. (2026). Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual Perturbation Screening at Scale. *bioRxiv*. https://doi.org/10.64898/2026.08.10.743942

## 1-line install

```bash
pip install virtualcellworkbench
```

*(placeholder; pre-publication release)*

## 30-second demo

```python
from virtualcellworkbench import load_cache, E1Model, predict
import torch

cache = load_cache("data/g2cp_cache_beta_all")
model = E1Model().load_pretrained("checkpoints/e1_seed0.pt")
pred = predict(model, cache, drug="BRD-K02130563", cell="A375")
print(pred.shape)  # (12328,) predicted gene-expression vector
```

## What it does

Predicts the genome-wide (12,328-gene) transcriptional response of a given compound in a given cell line, using only:

- **ECFP4 fingerprint** (2,048-bit Morgan, count-based) of the compound
- **Cell-line identity** (one of 162 LINCS / GEO core lines)

The training data covers **32,039 compounds** × **123 cell lines** × genome-wide expression, drawn from the LINCS L1000 level-5-beta release (Centerspace `trt_cp`).

## Headline numbers (paper Table 1, protocol per §3.1)

| Setting | Pearson r (per-sample, 12,328 genes) |
|---|---|
| Our v7_raw on official 10% held-out perturbations | **0.341** [0.334, 0.348] |
| Reported in bioRxiv preprint | 0.305 |
| Random 10% drug holdout (this repo's E1 baseline) | 0.137 |
| Leave-3-cell-out (unseen cell lines) | 0.222 |

The 35% drop on unseen cell lines quantifies the central paper claim: **cell-line heterogeneity is the primary bottleneck** for virtual perturbation screening at scale.

## Repository layout

```
.
├── src/virtualcellworkbench/   # core Python package
│   ├── data/                   # cache loaders, fingerprinters
│   ├── models/                 # E1 (baseline), v7_raw (full G2CP)
│   ├── train/                  # training loops, losses, splits
│   ├── eval/                   # PCC, direction concordance, bootstrap CI
│   └── pipelines/              # E0–E11 reproduction scripts
├── examples/                   # runnable scripts (start here)
├── configs/                    # default hyperparameter YAMLs
├── docs/                       # calibration + design notes
└── tests/                      # unit tests (smoke + integration)
```

## Reproducing the paper

```bash
# 1. Build the gene × drug × cell cache (requires ~50 GB free)
python -m virtualcellworkbench.pipelines.build_cache \
    --source level5_beta --output data/g2cp_cache_beta_all

# 2. Train the E1 baseline (≈8 min on RTX 4070 8 GB)
python -m virtualcellworkbench.pipelines.train_e1 \
    --cache data/g2cp_cache_beta_all --out checkpoints/e1_seed0.pt

# 3. Run the 11-experiment evaluation suite (≈2 h on RTX 4070 8 GB)
python -m virtualcellworkbench.pipelines.run_all_evals \
    --cache data/g2cp_cache_beta_all \
    --models checkpoints/e1_seed0.pt \
    --out results/

# 4. Render tables and figures
python -m virtualcellworkbench.pipelines.render_paper_artifacts \
    --results results/ --out paper_artifacts/
```

See `docs/` for the per-experiment design notes and the original stage logs (kept for transparency).

## Data sources

| Resource | Source | Required? |
|---|---|---|
| LINCS L1000 level-5-beta (`trt_cp`) | Broad / clue.io | ✅ for training |
| Compound annotations (target, MoA) | `compoundinfo_beta.txt` | optional (used in E3b) |
| ChEMBL 37 (mechanism, target, activities) | [ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl/) | optional (used in E6, E3b) |
| DepMap CERES (genetic perturbation ground truth) | [depmap.org](https://depmap.org) | optional (genetic split only) |

The packaged 7.4 GB cache `data/g2cp_cache_beta_all/` already encodes the post-filter subset used in the paper; rebuilds from raw sources are deterministic but slow.

## Citing

```bibtex
@software{wei2026virtualcellworkbench,
  author = {Wei, Kairui},
  title  = {VirtualCellWorkbench: A toolkit for chemical-perturbation response prediction at cell-line scale},
  year   = {2026},
  url    = {https://github.com/<username>/VirtualCellWorkbench},
  doi    = {10.5281/zenodo.<placeholder>}
}

@article{wei2026heterogeneity,
  author  = {Wei, Kairui},
  title   = {Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual Perturbation Screening at Scale},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.08.10.743942}
}
```

## License

MIT — see `LICENSE`.

## Status

Pre-publication. Numbers above are from the bioRxiv preprint plus the 11-experiment reproduction suite described in `docs/`. Watch this repository for the journal-version tag.
