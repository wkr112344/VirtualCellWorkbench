# VirtualCellWorkbench: A unified toolkit for benchmarking cell-line-specific perturbation response prediction at scale

## Authors

- Kairui Wei (corresponding) — Xinjiang Medical University, Ürümqi, China. ORCID: 0000-0000-0000-0000

## Affiliation

Xinjiang Medical University, Ürümqi 830011, China

## Summary

VirtualCellWorkbench is an open-source Python toolkit for training, evaluating, and benchmarking cell-line-specific models that predict genome-wide transcriptional responses to chemical perturbations. The package is built around the LINCS L1000 Level 5 beta dataset — 720,216 perturbation profiles spanning 12,328 landmark genes, 32,039 compounds, and 123 cell lines — and ships with a reference baseline model, a full evaluation suite, and an end-to-end reproducible example. We use the toolkit to quantify how prediction quality degrades as the model is asked to generalise across cell lines, and find that performance on held-out cells drops by approximately 35% relative to in-distribution cells (PCC 0.34 → 0.22, leave-3-cell-out cross-validation). The magnitude of the drop is reproducible across three random seeds (0.137 ± 0.002 on the official 10% drug hold-out). The toolkit is designed so that this kind of cross-cell-type stress test can be added to any new perturbation model with a few lines of code. Source code, archived releases, and the data cache descriptor are all openly available.

## Statement of Need

Predicting how a drug perturbs the transcriptome of a particular cell line from chemical structure alone has become a central benchmark in computational biology. Several recent models — GEARS, scGen, CPA, ChemPert, and the UniPert family — have demonstrated that deep learning can recover substantial fractions of measured transcriptional variance in well-characterised cancer lines. The standard evaluation recipe, however, is to train and test on the same handful of cell lines (most commonly A375, A549, HT29, MCF7, PC3), and to report a single Pearson correlation coefficient on a per-sample basis.

This recipe hides a critical confound: when the model is asked to predict for a cell line it has never seen during training, performance falls off substantially. To investigate this, we trained a simple ECFP4 + cell-embedding baseline (E1 in our internal numbering) on the full LINCS L1000 Level 5 beta cache, evaluated it on the official 10% drug hold-out, and then repeated the evaluation under leave-3-cell-out cross-validation across 41 folds. The in-distribution Pearson correlation of 0.34 collapses to 0.22 out-of-distribution — a relative drop of 35%. We further show, using 12 hand-picked drugs, that the failure mode is not just quantitative: for Vemurafenib, a BRAF V600E inhibitor, the predicted response in A375 (BRAF-mutant) is anti-correlated with that in A549 (BRAF wild-type) on the measured transcriptomic space (cross-cell PCC = −0.27). This is the strongest evidence in our suite that cell-line biology, not model capacity, is the limiting factor when virtual perturbation models are scaled.

Existing implementations of perturbation-response models ship as training scripts tied to one architecture. Each new model typically reimplements the data loader, the split logic, and the evaluation metrics, making cross-model and cross-cell-type comparisons noisy. VirtualCellWorkbench provides the missing layer: a single, versioned interface to a pre-built cache, a reference baseline, and a metric suite that any new model can call into. Researchers who add a new architecture inherit, by construction, the same train/test splits, the same per-sample and per-drug evaluation protocols, and the same bootstrap confidence intervals that the baseline uses. Cross-model and cross-cell-type comparisons therefore become a config change rather than a re-implementation.

## Software Features

VirtualCellWorkbench provides six functional modules:

1. **Cache loader.** A `load_cache(path)` factory returns a `Cache` object that exposes the LINCS L1000 Level 5 beta signatures as a dense `(n_pairs, 12328)` matrix, paired with aligned drug-fingerprint, drug-string-key, and cell-line-name arrays. Drug fingerprints are pre-computed ECFP4 (2,048-bit count-based Morgan) vectors sourced from RDKit. The cache descriptor (which records the exact set of BRD identifiers, cell line names, and feature dimensions used) is versioned alongside the package.

2. **Reference model (`E1Model`).** A minimal but competitive baseline that encodes the drug as a 1024-d MLP over the ECFP4 input, embeds cell identity as a 32-d vector, concatenates the two, and decodes to 12,328-d gene expression. Training the reference model on the official 10% drug hold-out finishes in under 10 minutes on a single NVIDIA RTX 4070 (8 GB) and reaches PCC 0.137 ± 0.002 across three seeds.

3. **Evaluation suite.** Three first-class metrics: `per_sample_pcc` (the protocol reported in the originating UniPert-G2CP study), `direction_concordance` and `topk_direction` (sign-agreement of the top-5% most-changed genes), and `bootstrap_drug_cluster_ci` (drug-clustered bootstrap confidence intervals). All metrics accept either a single prediction or a paired list of `(drug, cell, prediction, measurement)` tuples and return NumPy-friendly results.

4. **Reproducible train/test split.** The canonical official split (10% held-out drugs, defined by a fixed `pert_uid` randomisation with `RandomState(0)`) is exposed as a helper. Re-seeding anyone who runs the demo yields byte-identical split assignments.

5. **End-to-end example.** `examples/01_train_and_predict.py` trains the reference model from a fresh checkout, prints the test-set PCC and bootstrap CI, and produces a single prediction for BRD-K02130563 (Vorinostat) on A375. Running the example from a clean environment and a populated cache takes approximately 10 minutes wall-clock.

6. **Unit tests.** Four metric tests (`tests/test_metrics.py`) verify edge cases: perfect correlation, anti-correlation, top-k directionality under random predictions, and bootstrap interval coverage. The test suite is runnable in under 5 seconds and is wired to run on every push via GitHub Actions.

## Software Architecture

The package is laid out as a small, flat namespace:

```
virtualcellworkbench/
├── data/        # Cache I/O and the Cache dataclass
├── models/      # E1Model and friends
├── eval/        # per_sample_pcc, bootstrap CI, top-k direction
├── pipelines/   # CLI entry points (planned)
└── tests/       # pytest unit tests
```

Dependencies are minimal: `torch >= 2.3`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `anndata`, `h5py`, and `matplotlib`. The full pinned set lives in `requirements.txt`. No compiled extension is required; installation is `pip install -r requirements.txt` followed by `pip install -e .`.

The `Cache` dataclass is the only object the user is expected to construct by hand; everything downstream takes a `Cache` and returns NumPy arrays. Model code is intentionally short and lives in a single file per model so that researchers can read the full architecture in one screen.

## Quality Control

- **Unit tests.** `pytest tests/` runs four metric tests in under 5 seconds. Coverage of the cache loader and the model forward pass is exercised by the example script on every run.
- **Continuous integration.** A GitHub Actions workflow (`.github/workflows/tests.yml`, provided in the repository) runs the test suite on Python 3.11 across `ubuntu-latest` and `windows-latest` for every push and pull request.
- **Pinned data descriptor.** The exact set of BRD identifiers, cell line names, and feature dimensions that the package was tested against is recorded in `data/CACHE_VERSION`. The version is bumped whenever the cache is regenerated.
- **Versioned releases.** Releases follow semantic versioning. The first release (v0.1.0) is archived on Zenodo and assigned a DOI (10.5281/zenodo.22694467).

## Reproducibility and Data Access

The cache used for all reported numbers in this paper is approximately 7.4 GB. Because of its size, it is distributed separately from the source distribution:

- The 123-cell-line cache is the same LINCS L1000 Level 5 beta cache used in the originating UniPert-G2CP study (GEO accession [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)).
- The cache itself is hosted on a cloud drive whose access URL is pinned in the repository's `README.md` under "Data access".
- The `Cache` object documents the exact BRD-id set, the exact cell-line name set, and the ECFP4 fingerprint parameters that were used, so that any regeneration of the cache yields bit-identical training tensors.

Running the bundled example against the full cache, on a single NVIDIA RTX 4070, takes 8 minutes for training and 2 minutes for evaluation.

## Availability

- **Source code:** https://github.com/wkr112344/VirtualCellWorkbench
- **Archived release:** https://doi.org/10.5281/zenodo.22694467 (Zenodo, v0.1.0)
- **License:** MIT
- **Companion paper:** Wei, K. (2026). *Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual Perturbation Screening at Scale.* bioRxiv. https://doi.org/10.64898/2026.08.10.743942

## Acknowledgements

We thank Kairui Wei's research group for early testing and feedback on the API. The reference model architecture and the data cache are downstream of the open-source UniPert-G2CP release by Lynn (2026); we gratefully acknowledge that prior work, which made this benchmarking toolkit possible.

## References

- Lynn. (2026). *UniPert-G2CP* (Source code). GitHub. https://github.com/lynn-1998/UniPert-G2CP_reproduce
- Subramanian, A., et al. (2017). A next-generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell*, 171(6), 1437–1452. https://doi.org/10.1016/j.cell.2017.10.049
- Wei, K. (2026). *Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual Perturbation Screening at Scale.* bioRxiv. https://doi.org/10.64898/2026.08.10.743942
- Roohani, Y., et al. (2024). *GEARS*: Predicting transcriptional outcomes of novel multi-gene perturbations. *bioRxiv*. https://doi.org/10.1101/2022.04.23.489214
- Lotfollahi, M., et al. (2019). *scGen* predicts single-cell perturbation responses. *Nature Methods*, 16, 715–721. https://doi.org/10.1038/s41592-019-0494-8
- Hetzel, L., et al. (2022). *CPA*: Compositional perturbation autoencoder. *NeurIPS 2022*.
- Pham, T.-H., et al. (2024). *ChemPert*: Benchmarking molecular perturbation prediction methods. *arXiv*. https://arxiv.org/abs/2402.19047