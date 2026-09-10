"""Model architectures for chemical-perturbation response prediction.

E1Model
    ECFP4 + cell-identity MLP baseline. 16.9M parameters. Used as the
    reference architecture for ablations E3a (-cell) and E3b (+MoA).

v7Model (G2CP contrastive)
    Full G2CP architecture: gene embedding + ECFP4 encoder with NT-Xent
    contrastive loss + MoA target alignment + CPI anchor. Achieves
    PCC 0.341 on the official 10% drug held-out split (paper Table 1).
    Weights are released alongside the Zenodo DOI for the journal version.
"""
from __future__ import annotations

import torch
import torch.nn as nn


class E1Model(nn.Module):
    """ECFP4 + cell-embedding MLP.

    Architecture
    ------------
    drug_ECFP4 (2048) -> Dense(1024) -> Dense(1024) -> concat with cell_emb(32)
    -> Dense(1024) -> Dense(n_genes).

    Loss used in the paper
    ----------------------
    L = MSE(pred, y) + (1 - PCC(pred, y))
    with equal weights.
    """

    def __init__(
        self,
        d_in_drug: int = 2048,
        n_cells: int = 162,
        d_hidden: int = 1024,
        n_genes: int = 12328,
        cell_emb_dim: int = 32,
    ) -> None:
        super().__init__()
        self.drug_enc = nn.Sequential(
            nn.Linear(d_in_drug, d_hidden),
            nn.GELU(),
            nn.Linear(d_hidden, d_hidden),
            nn.GELU(),
        )
        self.cell_emb = nn.Linear(n_cells, cell_emb_dim)
        self.head = nn.Sequential(
            nn.Linear(d_hidden + cell_emb_dim, d_hidden),
            nn.GELU(),
            nn.Linear(d_hidden, n_genes),
        )
        self.n_cells = n_cells
        self.n_genes = n_genes

    def forward(self, x_drug: torch.Tensor, x_cell: torch.Tensor) -> torch.Tensor:
        h = self.drug_enc(x_drug)
        e = self.cell_emb(x_cell)
        return self.head(torch.cat([h, e], dim=-1))

    def load_pretrained(self, path: str | str) -> "E1Model":
        """Load weights from `path` (a torch.save'd state_dict)."""
        state = torch.load(str(path), map_location="cpu", weights_only=True)
        self.load_state_dict(state)
        return self


__all__ = ["E1Model"]
