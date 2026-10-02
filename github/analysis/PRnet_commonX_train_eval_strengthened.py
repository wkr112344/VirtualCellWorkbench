#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRnet_commonX_train_eval_strengthened.py — 「完整A」受控 2×2 的 *加强版* 单 seed 探针。

目的：验证「把 PRnet 容量加大 + 决定性评测(n=0) 后，val_mse 能否打穿 1.07 的方差地板」。
如果能在单 seed 上把 val_mse 压到 < 1.0，说明 LINCS 的 X->Y 映射是可学的，再决定是否全量重跑 5-seed。
如果仍然卡在 ~1.07，则基本可判定是数据本身的真 null / 弱信号，加强模型也无济于事。

与 multiseed 脚本保持完全一致的两点：
  - 数据 split / 四格定义 / corr 计算逻辑 一字不改（同 PRnet_commonX_dualY.h5）
  - 超参对比基线：train-cap 默认 200000（与全量多 seed 一致，隔离变量）

相对 multiseed 的 *加强* 改动（全部集中在本文件顶部常量）：
  - 容量：HIDDEN [128] -> [512,256]，Z_DIM 64 -> 256，ADAPTOR [128] -> [512,256]，C_DIM 64 -> 128
  - 评测决定性：predict / val 用 n=0（原来是 torch.randn -> 注入噪声、压低 correlation）
  - 训练：去掉 ReduceLROnPlateau（val 卡住会把 LR 砍到 min_lr），改固定 LR + 更长 patience
  - 轻微正则：DR_RATE 0.05->0.1，WD 1e-8->1e-6（更大模型防过拟合）

须用含 torch 的 python 运行（dpb311）。
用法：
  python PRnet_commonX_train_eval_strengthened.py --seeds 1 --seed-base 2024 --train-cap 200000
"""
import argparse, os, random
import h5py
import numpy as np
import torch
import torch.nn as nn
from scipy.stats import rankdata

# ============================================================ PRnet 模型 (PGM) —— 与 multiseed 同结构
class PEncoder(nn.Module):
    def __init__(self, layer_sizes, z_dim, dr):
        super().__init__()
        self.FC = None
        if len(layer_sizes) > 1:
            self.FC = nn.Sequential()
            for i, (ins, outs) in enumerate(zip(layer_sizes[:-1], layer_sizes[1:])):
                if i == 0:
                    self.FC.add_module(f"L{i}", nn.Linear(ins, outs, bias=False))
                else:
                    self.FC.add_module(f"L{i}", nn.Linear(ins, outs))
                    self.FC.add_module(f"N{i}", nn.BatchNorm1d(outs))
                    self.FC.add_module(f"A{i}", nn.LeakyReLU(0.3))
                    self.FC.add_module(f"D{i}", nn.Dropout(p=dr))
        self.mean_encoder = nn.Linear(layer_sizes[-1], z_dim)
    def forward(self, x):
        if self.FC is not None:
            x = self.FC(x)
        return self.mean_encoder(x)

class PDecoder(nn.Module):
    def __init__(self, z_dim, layer_sizes, x_dim, dr):
        super().__init__()
        layer_sizes = [z_dim] + layer_sizes
        self.FirstL = nn.Sequential()
        self.FirstL.add_module("L0", nn.Linear(layer_sizes[0], layer_sizes[1], bias=False))
        self.FirstL.add_module("N0", nn.BatchNorm1d(layer_sizes[1]))
        self.FirstL.add_module("A0", nn.LeakyReLU(0.3))
        self.FirstL.add_module("D0", nn.Dropout(p=dr))
        if len(layer_sizes) > 2:
            self.HiddenL = nn.Sequential()
            for i, (ins, outs) in enumerate(zip(layer_sizes[1:-1], layer_sizes[2:])):
                if i + 3 < len(layer_sizes):
                    self.HiddenL.add_module(f"L{i+1}", nn.Linear(ins, outs, bias=False))
                    self.HiddenL.add_module(f"N{i+1}", nn.BatchNorm1d(outs, affine=True))
                    self.HiddenL.add_module(f"A{i+1}", nn.LeakyReLU(0.3))
                    self.HiddenL.add_module(f"D{i+1}", nn.Dropout(p=dr))
        else:
            self.HiddenL = None
        self.recon_decoder = nn.Sequential(nn.Linear(layer_sizes[-2], layer_sizes[-1]))
        self.relu = nn.ReLU()
    def forward(self, z):
        d = self.FirstL(z)
        x = self.HiddenL(d) if self.HiddenL is not None else d
        r = self.recon_decoder(x)
        dim = r.size(1) // 2
        return torch.cat((self.relu(r[:, :dim]), r[:, dim:]), 1)

class PAdaptor(nn.Module):
    def __init__(self, layer_sizes, comb_dim, dr):
        super().__init__()
        self.FC = None
        if len(layer_sizes) > 1:
            self.FC = nn.Sequential()
            for i, (ins, outs) in enumerate(zip(layer_sizes[:-1], layer_sizes[1:])):
                if i == 0:
                    self.FC.add_module(f"L{i}", nn.Linear(ins, outs, bias=False))
                else:
                    self.FC.add_module(f"L{i}", nn.Linear(ins, outs))
                    self.FC.add_module(f"N{i}", nn.BatchNorm1d(outs))
                    self.FC.add_module(f"A{i}", nn.LeakyReLU(0.3))
                    self.FC.add_module(f"D{i}", nn.Dropout(p=dr))
        self.comb_encoder = nn.Linear(layer_sizes[-1], comb_dim)
    def forward(self, x):
        if self.FC is not None:
            x = self.FC(x)
        return self.comb_encoder(x)

class PGM(nn.Module):
    def __init__(self, x_dim, c_dim, n_dim, hidden, z_dim, adaptor, comb_adapt_dim, dr):
        super().__init__()
        enc = hidden.copy(); enc.insert(0, x_dim + c_dim)
        dec = hidden.copy(); dec.reverse(); dec.append(x_dim * 2)
        adp = adaptor.copy(); adp.insert(0, comb_adapt_dim)
        self.encoder = PEncoder(enc, z_dim, dr)
        self.decoder = PDecoder(z_dim + c_dim + n_dim, dec, x_dim, dr)
        self.adaptor = PAdaptor(adp, c_dim, dr)
        self.n_dim = n_dim
    def forward(self, x, c, n):
        c = self.adaptor(c)
        z = self.encoder(torch.cat((x, c), 1))
        return self.decoder(torch.cat((z, c, n), 1))

# ============================================================ 加强版超参
X_DIM = 969
HIDDEN = [512, 256]      # 原 [128]
Z_DIM = 256              # 原 64
ADAPTOR = [512, 256]     # 原 [128]
C_DIM = 128              # 原 64
DRUG_DIM = 1026          # 1024 (Morgan) + dose_z + time_z
N_DIM = 10
DR_RATE = 0.1            # 原 0.05
LR = 1e-3
WD = 1e-6                # 原 1e-8
BATCH = 512
N_EPOCHS = 300
EARLY_PATIENCE = 40      # 原 20
CHUNK = 20000
# 评测决定性：True 时 n=0（原来训练/评测都用 randn）
EVAL_DETERMINISTIC = True

def make_model(device, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    return PGM(X_DIM, C_DIM, N_DIM, HIDDEN, Z_DIM, ADAPTOR, DRUG_DIM, DR_RATE).to(device)

def n_tensor(bs, device, deterministic):
    if deterministic:
        return torch.zeros(bs, N_DIM, device=device)
    return torch.randn(bs, N_DIM, device=device)

def weight_init(m):
    if isinstance(m, nn.Linear):
        m.weight.data.normal_(0, 0.02)
        if m.bias is not None: m.bias.data.zero_()
    elif isinstance(m, nn.BatchNorm1d):
        m.weight.data.normal_(1, 0.02)
        if m.bias is not None: m.bias.data.zero_()

# ============================================================ 流式 batch 生成器
def iter_batches(h5, idx, Yname, dose_z, time_z, batch, shuffle, device, with_y=True):
    sidx = np.sort(idx)
    for s in range(0, len(sidx), CHUNK):
        blk = sidx[s:s+CHUNK]
        Xb = h5["X_baseline"][blk].astype(np.float32)
        de = h5["drug_emb"][blk].astype(np.float32)
        dz = dose_z[blk][:, None]; tz = time_z[blk][:, None]
        cond = np.concatenate([de, dz, tz], 1).astype(np.float32)
        Y = h5[Yname][blk].astype(np.float32) if with_y else None
        perm = np.arange(len(blk))
        if shuffle:
            np.random.shuffle(perm)
        for b in range(0, len(blk), batch):
            pi = perm[b:b+batch]
            xb_t = torch.from_numpy(Xb[pi]).to(device)
            cond_t = torch.from_numpy(cond[pi]).to(device)
            if with_y:
                y_t = torch.from_numpy(Y[pi]).to(device)
                yield xb_t, cond_t, y_t
            else:
                yield xb_t, cond_t

# ============================================================ 向量化逐 sample 相关
def corr_rows(pred, true, method):
    pred = np.asarray(pred, dtype=np.float64)
    true = np.asarray(true, dtype=np.float64)
    if method == "pearson":
        pm = pred - pred.mean(1, keepdims=True)
        tm = true - true.mean(1, keepdims=True)
    else:  # spearman: rank then pearson
        pm = np.apply_along_axis(rankdata, 1, pred)
        tm = np.apply_along_axis(rankdata, 1, true)
        pm = pm - pm.mean(1, keepdims=True)
        tm = tm - tm.mean(1, keepdims=True)
    num = (pm * tm).sum(1)
    den = np.sqrt((pm**2).sum(1) * (tm**2).sum(1))
    return np.where(den > 1e-12, num / den, np.nan)

# ============================================================ 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=r"C:/Users/wkr20/WorkBuddy/Claw/gigascience_v5/data/prnet_commonX/PRnet_commonX_dualY.h5")
    ap.add_argument("--outdir", default=r"C:/Users/wkr20/WorkBuddy/Claw/gigascience_v5/data/prnet_commonX/strengthened")
    ap.add_argument("--epochs", type=int, default=N_EPOCHS)
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--seed-base", type=int, default=2024)
    ap.add_argument("--train-cap", type=int, default=200000, help="训练样本子采样上限 (0=全量)")
    ap.add_argument("--limit-eval", type=int, default=None)
    ap.add_argument("--tag", default="strengthened")
    ap.add_argument("--loss", default="nll", choices=["nll", "mse"],
                    help="nll=GaussianNLLLoss(含 logvar 头); mse=仅对 mean 做 MSE(逼模型学均值, 这才是评测用的量)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    nparam = 0
    print(f"[env] device={device} seeds={args.seeds} seed_base={args.seed_base} train_cap={args.train_cap} "
          f"HIDDEN={HIDDEN} Z_DIM={Z_DIM} C_DIM={C_DIM} ADAPTOR={ADAPTOR} "
          f"eval_deterministic={EVAL_DETERMINISTIC}")

    with h5py.File(args.data, "r") as f:
        split = f["split"][:].astype(int)
        dose = f["dose"][:].astype(np.float64)
        time = f["time"][:].astype(np.float64)
        n = len(split)
    dm, ds = dose.mean(), dose.std() or 1.0
    tm_, ts = time.mean(), time.std() or 1.0
    dose_z = ((dose - dm) / ds).astype(np.float32)
    time_z = ((time - tm_) / ts).astype(np.float32)

    train_idx = np.where(split == 0)[0]
    if args.train_cap and 0 < args.train_cap < len(train_idx):
        rng_t = np.random.default_rng(12345)
        train_idx = rng_t.choice(train_idx, args.train_cap, replace=False)
        print(f"[train] capped to {len(train_idx)} samples (fixed across seeds)")
    else:
        print(f"[train] full {len(train_idx)} samples")

    val_idx = np.where(split == 1)[0]
    if len(val_idx) > 60000:
        rng_v = np.random.default_rng(2024)
        val_idx = rng_v.choice(val_idx, 60000, replace=False)
    eval_idx = np.where(split == 3)[0]
    if args.limit_eval:
        eval_idx = eval_idx[:args.limit_eval]
    print(f"[split] n={n} train={len(train_idx)} val={len(val_idx)} eval(S25)={len(eval_idx)}")

    with h5py.File(args.data, "r") as f:
        Yb = f["Y_beta"][eval_idx]; Yd = f["Y_dcic"][eval_idx]
    M = len(eval_idx)

    criterion = nn.GaussianNLLLoss() if args.loss == "nll" else nn.MSELoss()
    print(f"[loss] using {args.loss.upper()} criterion")

    def train_arm(Yname, save_path, seed):
        print(f"\n===== TRAIN ARM  target={Yname}  seed={seed} =====")
        model = make_model(device, seed); model.apply(weight_init)
        if nparam == 0:
            global _np
            _np = sum(p.numel() for p in model.parameters())
            print(f"[model] params={_np:,}")
        opt = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=LR, weight_decay=WD)
        best_mse, patient, best_sd = np.inf, 0, None
        with h5py.File(args.data, "r") as h5:
            for ep in range(args.epochs):
                model.train(); tot = 0.0; cnt = 0
                for xb, cond, y in iter_batches(h5, train_idx, Yname, dose_z, time_z, BATCH, True, device, True):
                    opt.zero_grad()
                    out = model(xb, cond, n_tensor(xb.size(0), device, False))
                    dim = out.size(1)//2
                    mean = out[:, :dim]; logvar = torch.nn.functional.softplus(out[:, dim:])
                    if args.loss == "mse":
                        loss = criterion(mean, y)
                    else:
                        loss = criterion(mean, y, logvar)
                    loss.backward(); opt.step(); tot += loss.item()*xb.size(0); cnt += xb.size(0)
                model.eval(); vmse = 0.0; vc = 0
                with torch.no_grad():
                    for xb, cond, y in iter_batches(h5, val_idx, Yname, dose_z, time_z, BATCH, False, device, True):
                        out = model(xb, cond, n_tensor(xb.size(0), device, EVAL_DETERMINISTIC))
                        dim = out.size(1)//2; mean = out[:, :dim]
                        vmse += float(torch.nn.functional.mse_loss(mean, y).item())*xb.size(0); vc += xb.size(0)
                vmse /= max(vc, 1)
                print(f"  epoch {ep:3d} train_loss={tot/cnt:.4f} val_mse={vmse:.4f} best={best_mse:.4f}")
                if vmse < best_mse:
                    best_mse = vmse; patient = 0; best_sd = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                elif patient <= EARLY_PATIENCE:
                    patient += 1
                else:
                    print("  early stop"); break
        if best_sd is not None:
            torch.save(best_sd, save_path); print(f"  saved -> {save_path}")
        return best_sd

    def predict_arm(sd, idx):
        model = make_model(device, 0); model.load_state_dict(sd); model.eval()
        preds = []
        with h5py.File(args.data, "r") as h5:
            with torch.no_grad():
                for xb, cond in iter_batches(h5, idx, "Y_beta", dose_z, time_z, BATCH*2, False, device, False):
                    out = model(xb, cond, n_tensor(xb.size(0), device, EVAL_DETERMINISTIC))
                    preds.append(out[:, :out.size(1)//2].detach().cpu().numpy())
        return np.concatenate(preds, 0)

    # ---- 逐 seed 跑 2×2 ----
    BB_sp, BD_sp, DB_sp, DD_sp = [], [], [], []
    BB_pe, BD_pe, DB_pe, DD_pe = [], [], [], []
    seed_interactions = []
    for k in range(args.seeds):
        seed = args.seed_base + k
        tag = f"{args.tag}_seed{k}"
        sd_beta = train_arm("Y_beta", os.path.join(args.outdir, f"{tag}_arm_beta.pt"), seed)
        sd_dcic = train_arm("Y_dcic", os.path.join(args.outdir, f"{tag}_arm_dcic.pt"), seed)
        pred_beta = predict_arm(sd_beta, eval_idx)
        pred_dcic = predict_arm(sd_dcic, eval_idx)
        bb_sp = corr_rows(pred_beta, Yb, "spearman"); bd_sp = corr_rows(pred_beta, Yd, "spearman")
        db_sp = corr_rows(pred_dcic, Yb, "spearman"); dd_sp = corr_rows(pred_dcic, Yd, "spearman")
        bb_pe = corr_rows(pred_beta, Yb, "pearson");  bd_pe = corr_rows(pred_beta, Yd, "pearson")
        db_pe = corr_rows(pred_dcic, Yb, "pearson");  dd_pe = corr_rows(pred_dcic, Yd, "pearson")
        del pred_beta, pred_dcic
        BB_sp.append(bb_sp); BD_sp.append(bd_sp); DB_sp.append(db_sp); DD_sp.append(dd_sp)
        BB_pe.append(bb_pe); BD_pe.append(bd_pe); DB_pe.append(db_pe); DD_pe.append(dd_pe)
        int_sp = float(np.nanmean((bb_sp - bd_sp) - (db_sp - dd_sp)))
        int_pe = float(np.nanmean((bb_pe - bd_pe) - (db_pe - dd_pe)))
        seed_interactions.append((int_sp, int_pe))
        torch.cuda.empty_cache() if device.type == "cuda" else None
        print(f"[seed {seed}] BB_sp={np.nanmean(bb_sp):+.4f} BD_sp={np.nanmean(bd_sp):+.4f} "
              f"DB_sp={np.nanmean(db_sp):+.4f} DD_sp={np.nanmean(dd_sp):+.4f} interaction_sp={int_sp:+.4f}")

    rows = []
    for k in range(args.seeds):
        rows.append(dict(
            seed=args.seed_base + k,
            BB_sp=np.nanmean(BB_sp[k]), BD_sp=np.nanmean(BD_sp[k]),
            DB_sp=np.nanmean(DB_sp[k]), DD_sp=np.nanmean(DD_sp[k]),
            interaction_sp=seed_interactions[k][0],
            BB_pe=np.nanmean(BB_pe[k]), BD_pe=np.nanmean(BD_pe[k]),
            DB_pe=np.nanmean(DB_pe[k]), DD_pe=np.nanmean(DD_pe[k]),
            interaction_pe=seed_interactions[k][1],
        ))
    import csv
    sum_csv = os.path.join(args.outdir, f"{args.tag}_2x2_summary.csv")
    with open(sum_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    med_sp = float(np.nanmean([s[0] for s in seed_interactions]))
    med_pe = float(np.nanmean([s[1] for s in seed_interactions]))
    print("\n================ STRENGTHENED SEED RESULT (spearman) ================")
    for r in rows:
        print(f"  seed {r['seed']}: BB={r['BB_sp']:+.4f} BD={r['BD_sp']:+.4f} DB={r['DB_sp']:+.4f} "
              f"DD={r['DD_sp']:+.4f} int={r['interaction_sp']:+.4f}")
    print(f"  mean interaction_sp = {med_sp:+.4f}   interaction_pe = {med_pe:+.4f}")
    print(f"  (baseline 对比: 原 multiseed val_mse 卡 ~1.07 -> 看上面 train_arm 的 best val_mse 能否 < 1.0)")
    print(f"  -> {sum_csv}")

if __name__ == "__main__":
    main()
