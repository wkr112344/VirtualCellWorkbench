"""Data loaders for the LINCS L1000 level-5-beta cache.

The cache is built once by `pipelines.build_cache` and contains:
  - y.npy            (157522, 12328)  float32  per-pair measured expression
  - drug_fps.npy     (32039, 2048)    float32  ECFP4 fingerprints (count-based Morgan)
  - meta.npz         per-pair (cell_idx, drug_idx, kind) and global vocabs
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass
class Cache:
    y: np.ndarray             # (n_pairs, n_genes)
    drug_fps: np.ndarray      # (n_drugs, fp_dim)
    kind: np.ndarray          # (n_pairs,) 1 = chemical, 0 = genetic
    key: np.ndarray           # (n_pairs,) drug or gene index
    cell: np.ndarray          # (n_pairs,) cell index
    hvg: np.ndarray           # (n_genes,) gene symbols
    cl_names: list            # (n_cells,) cell line names
    drug_vocab: list          # (n_drugs,) BRD-IDs

    @property
    def n_pairs(self) -> int:
        return self.y.shape[0]

    @property
    def n_genes(self) -> int:
        return self.y.shape[1]

    @property
    def n_cells(self) -> int:
        return len(self.cl_names)

    @property
    def n_drugs(self) -> int:
        return len(self.drug_vocab)

    def chem_mask(self) -> np.ndarray:
        """Boolean mask for chemical-perturbation pairs (kind == 1)."""
        return self.kind == 1

    def drug_index(self, brd_id: str) -> Optional[int]:
        try:
            return self.drug_vocab.index(brd_id)
        except ValueError:
            return None

    def cell_index(self, name: str) -> Optional[int]:
        try:
            return self.cl_names.index(name)
        except ValueError:
            return None


def load_cache(path: str | Path) -> Cache:
    """Load the prebuilt cache from `path` (the directory containing y.npy etc.).

    Parameters
    ----------
    path
        Directory produced by `pipelines.build_cache`. Files expected:
        y.npy, drug_fps.npy, meta.npz.

    Returns
    -------
    Cache
        Dataclass exposing the loaded arrays and lookups.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Cache directory not found: {p}")

    y = np.load(p / "y.npy", mmap_mode="r")
    drug_fps = np.load(p / "drug_fps.npy", mmap_mode="r")
    meta = np.load(p / "meta.npz", allow_pickle=True)

    return Cache(
        y=y,
        drug_fps=drug_fps,
        kind=meta["kind"],
        key=meta["key"],
        cell=meta["cell"],
        hvg=meta["hvg"],
        cl_names=list(meta["cl_names"]),
        drug_vocab=list(meta["drug_vocab"]),
    )


__all__ = ["Cache", "load_cache"]
