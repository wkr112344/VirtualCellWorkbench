"""Evaluation metrics: per-sample PCC, direction concordance, bootstrap CI."""
from __future__ import annotations

import numpy as np


def per_sample_pcc(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    """Pearson correlation per row across the last axis.

    Both inputs must share shape (N, G). Returns a length-N array.
    """
    p = pred - pred.mean(axis=1, keepdims=True)
    t = true - true.mean(axis=1, keepdims=True)
    num = (p * t).sum(axis=1)
    den = np.sqrt((p ** 2).sum(axis=1) * (t ** 2).sum(axis=1) + 1e-8)
    return num / den


def direction_concordance(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    """Per-row fraction of genes where sign(pred) == sign(true)."""
    return (np.sign(pred) == np.sign(true)).mean(axis=1)


def topk_direction(pred: np.ndarray, true: np.ndarray, k_frac: float = 0.05) -> np.ndarray:
    """Per-row sign agreement on the top-k_frac of |true|."""
    k = max(1, int(k_frac * pred.shape[1]))
    idx = np.argsort(-np.abs(true), axis=1)[:, :k]
    rows = np.arange(pred.shape[0])[:, None]
    return (np.sign(pred[rows, idx]) == np.sign(true[rows, idx])).mean(axis=1)


def bootstrap_drug_cluster_ci(
    per_pair_pcc: np.ndarray,
    drugs: np.ndarray,
    n_boot: int = 200,
    seed: int = 0,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Drug-cluster bootstrap 95% CI for the mean per-sample PCC.

    Resamples drugs with replacement, averages the per-drug-mean PCC, and
    returns the [alpha/2, 1-alpha/2] percentile interval.
    """
    rng = np.random.RandomState(seed)
    unique_drugs = np.unique(drugs)
    boots = []
    for _ in range(n_boot):
        ds = rng.choice(unique_drugs, size=len(unique_drugs), replace=True)
        pccs = []
        for d in ds:
            mask = drugs == d
            if mask.sum() == 0:
                continue
            pccs.append(per_pair_pcc[mask].mean())
        boots.append(np.mean(pccs))
    lo = float(np.percentile(boots, 100 * alpha / 2))
    hi = float(np.percentile(boots, 100 * (1 - alpha / 2)))
    return lo, hi


__all__ = [
    "per_sample_pcc",
    "direction_concordance",
    "topk_direction",
    "bootstrap_drug_cluster_ci",
]
