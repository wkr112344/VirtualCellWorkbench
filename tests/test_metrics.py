"""Smoke tests for the public API.

Run with:
    pytest tests/
or
    python -m unittest discover tests
"""
from __future__ import annotations

import unittest

import numpy as np

from virtualcellworkbench.eval import (
    bootstrap_drug_cluster_ci,
    per_sample_pcc,
    topk_direction,
)


class TestMetrics(unittest.TestCase):
    def test_per_sample_pcc_perfect(self):
        x = np.random.randn(10, 100)
        pcc = per_sample_pcc(x, x)
        self.assertTrue(np.allclose(pcc, 1.0, atol=1e-5))

    def test_per_sample_pcc_anticorrelated(self):
        x = np.random.randn(10, 100)
        pcc = per_sample_pcc(x, -x)
        self.assertTrue(np.allclose(pcc, -1.0, atol=1e-5))

    def test_topk_direction(self):
        true = np.array([[1, -2, 3, -4, 5]], dtype=float)
        pred = np.array([[1, -2, 3, -4, 5]], dtype=float)
        self.assertAlmostEqual(topk_direction(pred, true, k_frac=0.4)[0], 1.0)

    def test_bootstrap_returns_interval(self):
        rng = np.random.RandomState(0)
        pcc = rng.uniform(0, 1, size=200)
        drugs = rng.randint(0, 20, size=200)
        lo, hi = bootstrap_drug_cluster_ci(pcc, drugs, n_boot=50, seed=0)
        self.assertLess(lo, hi)


if __name__ == "__main__":
    unittest.main()
