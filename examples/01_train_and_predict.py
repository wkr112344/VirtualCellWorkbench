"""End-to-end demo: load cache, train E1, predict one drug×cell pair.

This script doubles as a self-contained smoke test. Total runtime:
~2 minutes on an RTX 4070 (8 GB) once the cache is available.

Run with:
    python examples/01_train_and_predict.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# Allow running from the repo root without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from virtualcellworkbench import E1Model, load_cache
from virtualcellworkbench.eval import (
    bootstrap_drug_cluster_ci,
    per_sample_pcc,
)


def batch_pcc_loss(pred, true):
    p = pred - pred.mean(dim=1, keepdim=True)
    t = true - true.mean(dim=1, keepdim=True)
    num = (p * t).sum(dim=1)
    den = torch.sqrt((p ** 2).sum(dim=1) * (t ** 2).sum(dim=1) + 1e-8)
    return 1.0 - (num / den).mean()


def official_split(kind, key, n_genes, heldout_frac=0.10, seed=0):
    """RandomState(0) + last 10% of perturbation UIDs, the paper's protocol."""
    pert_uid = np.where(kind == 1, n_genes + key, -1)
    uids = np.unique(pert_uid[pert_uid >= 0])
    rng = np.random.RandomState(seed)
    rng.shuffle(uids)
    n_test = int(len(uids) * heldout_frac)
    test = set(uids[-n_test:].tolist())
    test_mask = chem := (kind == 1) & np.isin(pert_uid, list(test))
    train_mask = (kind == 1) & ~np.isin(pert_uid, list(test))
    return train_mask, test_mask


def main(cache_path: str = "data/g2cp_cache_beta_all", epochs: int = 8) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}")

    print("loading cache ...")
    cache = load_cache(cache_path)
    print(f"  {cache.n_pairs} pairs, {cache.n_cells} cells, {cache.n_drugs} drugs, {cache.n_genes} genes")

    print("splitting (paper §3.1 protocol) ...")
    train_mask, test_mask = official_split(cache.kind, cache.key, len(cache.hvg))
    print(f"  train: {train_mask.sum()}  test: {test_mask.sum()}")

    # Build tensors
    drug_fps = np.asarray(cache.drug_fps, dtype=np.float32)
    cell_oh = np.eye(cache.n_cells, dtype=np.float32)
    X_d_tr = torch.from_numpy(drug_fps[cache.key[train_mask]])
    X_c_tr = torch.from_numpy(cell_oh[cache.cell[train_mask]])
    Y_tr = torch.from_numpy(np.asarray(cache.y[train_mask], dtype=np.float32))
    X_d_te = torch.from_numpy(drug_fps[cache.key[test_mask]])
    X_c_te = torch.from_numpy(cell_oh[cache.cell[test_mask]])
    Y_te = torch.from_numpy(np.asarray(cache.y[test_mask], dtype=np.float32))

    # 90/10 train/val
    perm = np.random.RandomState(0).permutation(len(Y_tr))
    n_val = int(0.1 * len(Y_tr))
    val_idx, tr_idx = perm[:n_val], perm[n_val:]

    model = E1Model(n_cells=cache.n_cells, n_genes=cache.n_genes).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    t0 = time.time()
    for ep in range(epochs):
        model.train()
        perm = torch.randperm(len(tr_idx))
        for i in range(0, len(perm), 256):
            b = tr_idx[perm[i:i + 256]]
            xd = X_d_tr[b].to(device, non_blocking=True)
            xc = X_c_tr[b].to(device, non_blocking=True)
            y = Y_tr[b].to(device, non_blocking=True)
            pred = model(xd, xc)
            loss = F.mse_loss(pred, y) + batch_pcc_loss(pred, y)
            opt.zero_grad()
            loss.backward()
            opt.step()
        sched.step()
        print(f"  epoch {ep:02d}  loss={loss.item():.4f}  ({time.time()-t0:.1f}s)")

    # Test
    model.eval()
    test_pcc = []
    with torch.no_grad():
        for i in range(0, len(Y_te), 512):
            xd = X_d_te[i:i + 512].to(device)
            xc = X_c_te[i:i + 512].to(device)
            y = Y_te[i:i + 512]
            pred = model(xd, xc).cpu().numpy()
            test_pcc.extend(per_sample_pcc(pred, y.numpy()).tolist())
    test_pcc = np.array(test_pcc)
    lo, hi = bootstrap_drug_cluster_ci(test_pcc, cache.key[test_mask])
    print(f"\nE1 test PCC = {test_pcc.mean():.4f}  bootstrap 95% CI = [{lo:.4f}, {hi:.4f}]")

    # Example prediction
    drug = "BRD-K02130563"  # panobinostat
    cell = "A375"
    d_idx = cache.drug_index(drug)
    c_idx = cache.cell_index(cell)
    xd = torch.from_numpy(drug_fps[d_idx:d_idx + 1]).to(device)
    xc = torch.from_numpy(cell_oh[c_idx:c_idx + 1]).to(device)
    with torch.no_grad():
        pred = model(xd, xc).cpu().numpy()[0]
    print(f"\nPrediction for {drug} × {cell}:")
    print(f"  shape: {pred.shape}, dtype: {pred.dtype}")
    print(f"  top-5 |x| genes: {cache.hvg[np.argsort(-np.abs(pred))[:5]].tolist()}")


if __name__ == "__main__":
    main()
