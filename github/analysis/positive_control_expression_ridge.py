# DepMap positive control: expression-PCA -> multi-output Ridge
#   目的：在同一条 Figure 6 划分（train 636 / test 136）上，验证"真正依赖细胞状态特征"的
#   预测器在 residualized 评价下是否仍保留正信号，并与冻结模型对照。
#
# 指标定义（关键）：
#   by-cell raw        : 每个 test 细胞，跨基因 Pearson(pred, truth)
#   by-cell shift      : 同时减训练集 per-gene 均值 -> 【数学上与 raw 恒等】，仅作校验
#   by-cell zres       : (pred-mu)/sd vs (truth-mu)/sd，mu/sd 只用训练细胞 -> 真正的"去 gene-level 背景"
#   across-cell raw    : 每个基因，跨 test 细胞 Pearson(pred, truth)（对 per-gene 平移不变）
#   across-cell zres   : 同上但先做 per-gene z（等价于 raw，因为对 per-gene 仿射不变），做校验
#
# 对象：
#   mean_baseline (negative/基础对照)：预测 = 训练集 per-gene 均值（常数向量）
#   ridge_CERES / ridge_Chronos (positive control)：X_TPM -> PCA(100) -> Ridge，训练靶标分别为 CERES/Chronos
#   ridge_perm (label-permutation negative control)：打乱训练 Y 的行标签后同流程
#   frozen_G2CP (参照)：既有 depmap_analysis_ready 的 5-seed 平均预测，两臂各一份
#
# 输出：results/*.csv|json、figures/panel_ABC_positive_control.png
import os, sys, json, time
import numpy as np
from scipy.stats import rankdata   # noqa: F401
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AL = r"C:/Users/wkr20/Desktop/depmap_second_corpus/work/DeepDEP_training_data/aligned_21Q2_908_3omics"
FR = r"C:/Users/wkr20/WorkBuddy/Claw/depmap_analysis_ready_20260929"
R = lambda *a: os.path.join(BASE, *a)
N_PC = 100
ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0]
N_BOOT = 5000
SEED = 20260930
MIN_CELLS_PER_GENE = 30


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def load():
    X = np.load(os.path.join(AL, "X_TPM.npy")).astype(np.float32)
    Ycer = np.load(os.path.join(AL, "Y_CERES.npy")).astype(np.float32)
    Ychr = np.load(os.path.join(AL, "Y_Chronos.npy")).astype(np.float32)
    mask = np.load(os.path.join(AL, "common_observed_mask.npy"))
    split = json.load(open(os.path.join(AL, "split.json")))
    idx = json.load(open(os.path.join(FR, "index", "depmap_index.json")))
    return X, Ycer, Ychr, mask, split, idx


def corr_rows(P, T, M):
    """逐行（=逐细胞）跨列 Pearson。零方差预测按定义记 0（无线性信息）。"""
    out = np.full(P.shape[0], np.nan)
    for i in range(P.shape[0]):
        m = M[i] & np.isfinite(T[i]) & np.isfinite(P[i])
        if m.sum() < 30:
            continue
        p = P[i][m].astype(np.float64); t = T[i][m].astype(np.float64)
        p = p - p.mean(); t = t - t.mean()
        if (p * p).sum() == 0 or (t * t).sum() == 0:
            out[i] = 0.0
            continue
        d = np.sqrt((p * p).sum() * (t * t).sum())
        if d > 0:
            out[i] = (p * t).sum() / d
    return out


def sanitize_mu_sd(mu, sd, mtr, Ytr):
    """mu 的 NaN 用 0 代替、sd 的 NaN/0 用 1 代替 —— 保证 shift/zres 与 raw 用同一基因集。"""
    mu = np.nan_to_num(np.asarray(mu, dtype=np.float64), nan=0.0)
    sd = np.asarray(sd, dtype=np.float64)
    bad = (~np.isfinite(sd)) | (sd <= 0)
    sd = np.where(bad, 1.0, sd)
    return mu, sd


def corr_cols(P, T, M):
    """逐列（=逐基因）跨行 Pearson：对 per-gene 平移/缩放不变。返回 (G,) + 有效细胞数。"""
    n, G = P.shape
    out = np.full(G, np.nan); ns = np.zeros(G, np.int32)
    for g in range(G):
        m = M[:, g] & np.isfinite(T[:, g]) & np.isfinite(P[:, g])
        k = int(m.sum()); ns[g] = k
        if k < MIN_CELLS_PER_GENE:
            continue
        p = P[m, g].astype(np.float64); t = T[m, g].astype(np.float64)
        p = p - p.mean(); t = t - t.mean()
        d = np.sqrt((p * p).sum() * (t * t).sum())
        if d > 0:
            out[g] = (p * t).sum() / d
    return out, ns


def boot_median_ci(vals, cells, nboot=N_BOOT, seed=SEED):
    """按 test 细胞（=donor 单元）重采样，对 median 或 median-of-cell-medians 给 CI。"""
    rng = np.random.default_rng(seed)
    v = np.asarray(vals, dtype=np.float64)
    ok = np.isfinite(v)
    by = {}
    for c, x, o in zip(cells, v, ok):
        if o:
            by.setdefault(c, []).append(x)
    keys = sorted(by)
    if not keys:
        return float("nan"), float("nan"), float("nan"), False, np.array([])
    arrs = [np.array(by[k]) for k in keys]
    point = float(np.median(np.concatenate(arrs)))
    out = np.empty(nboot)
    for b in range(nboot):
        pick = rng.integers(0, len(keys), len(keys))
        out[b] = np.median(np.concatenate([arrs[i] for i in pick]))
    lo, hi = np.percentile(out, [2.5, 97.5])
    return point, float(lo), float(hi), bool(lo * hi > 0), out


def zres_stats(P, T, M, mu, sd):
    """三个口径：
       raw   : 原样
       shift : 同时减 per-gene 训练均值（用户原设计；注意 per-gene 常数不是全局常数，
               因此这会真实改变跨基因相关 —— 它剥掉的是跨基因共享的 gene-level background）
       z     : 再除以 per-gene 训练标准差（额外做 per-gene 尺度归一）
       返回 (r_raw, r_shift, r_z, corr(r_shift, r_z))"""
    Ps = P - mu; Ts = T - mu
    Pz = Ps / sd; Tz = Ts / sd
    r_raw = corr_rows(P, T, M)
    r_shift = corr_rows(Ps, Ts, M)
    r_z = corr_rows(Pz, Tz, M)
    o = np.isfinite(r_shift) & np.isfinite(r_z)
    agree = float(np.corrcoef(r_shift[o], r_z[o])[0, 1]) if o.sum() > 2 else float("nan")
    return r_raw, r_shift, r_z, agree


def main():
    t0 = time.time()
    X, Ycer, Ychr, mask, split, idx = load()
    log("载入: X%s Ycer%s mask%s" % (X.shape, Ycer.shape, mask.shape))

    tr = np.array(idx["train_rows_in_cell_order"], dtype=np.int64)
    te = np.array(idx["test_rows_in_cell_order"], dtype=np.int64)
    cells = [l.strip() for l in open(os.path.join(FR, "index", "cell_order.txt")) if l.strip()]
    te_cells = [cells[i] for i in te]
    assert [cells[i] for i in te] == [c for c in split["test_cells"]], "test 细胞顺序与 split.json 不一致"
    assert len(tr) == 636 and len(te) == 136
    log("train %d / test %d 细胞（与原 Figure 6 划分一致）" % (len(tr), len(te)))

    # --- 特征：训练集标准化 + PCA(100) ---
    sc = StandardScaler().fit(X[tr])
    Ztr = sc.transform(X[tr]); Zte = sc.transform(X[te])
    pca = PCA(n_components=N_PC, random_state=0).fit(Ztr)
    Ktr = pca.transform(Ztr).astype(np.float64); Kte = pca.transform(Zte).astype(np.float64)
    log("PCA100 训练集解释方差 %.3f" % pca.explained_variance_ratio_.sum())

    mtr = mask[tr]; mte = mask[te]
    preds, chosen = {}, {}

    def fit_ridge(Ytrain, tag):
        # alpha 只用训练集的 5-fold CV 选（test 全程不参与）
        kf = KFold(n_splits=5, shuffle=True, random_state=0)
        best, best_s = ALPHAS[0], -np.inf
        ytr = np.nan_to_num(Ytrain)
        for a in ALPHAS:
            s = []
            for i_tr, i_va in kf.split(Ktr):
                m = Ridge(alpha=a).fit(Ktr[i_tr], ytr[i_tr])
                ph = m.predict(Ktr[i_va])
                tv = ytr[i_va]
                ok = np.isfinite(tv)
                ss_res = ((ph - tv)[ok] ** 2).sum()
                ss_tot = ((tv[ok] - tv[ok].mean()) ** 2).sum()
                s.append(1 - ss_res / max(ss_tot, 1e-12))
            sc_ = float(np.mean(s))
            log("    %s alpha=%g CV-R2=%.4f" % (tag, a, sc_))
            if sc_ > best_s:
                best_s, best = sc_, a
        chosen[tag] = {"alpha": best, "cv_r2": best_s}
        mdl = Ridge(alpha=best).fit(Ktr, ytr)
        return mdl.predict(Kte), best

    log("训练 Ridge: CERES 臂")
    preds["ridge_CERES"], _ = fit_ridge(Ycer[tr], "ridge_CERES")
    log("训练 Ridge: Chronos 臂")
    preds["ridge_Chronos"], _ = fit_ridge(Ychr[tr], "ridge_Chronos")

    log("训练 Ridge: label-permutation 对照（打乱训练 Y 行标签）")
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(len(tr))
    preds["ridge_perm"], _ = fit_ridge(Ycer[tr][perm], "ridge_perm")

    # 冻结模型（既有 5-seed 平均，test 136 × 17393）
    for arm, Ytag in [("ceres", "CERES"), ("chronos", "Chronos")]:
        acc = None
        for k in range(5):
            f = os.path.join(FR, "predictions_test", "seed%d_pred_%s_arm_test.npy" % (k, arm))
            a = np.load(f).astype(np.float32)
            acc = a if acc is None else acc + a
        preds["frozen_%s" % Ytag] = (acc / 5.0).astype(np.float64)

    # --- 训练集 per-gene 均值/标准差（只用 train 细胞） ---
    def prep_mu_sd(Y):
        mu = np.nanmean(np.where(mtr, Y[tr], np.nan), axis=0)
        sd = np.sqrt(np.nanmean(np.where(mtr, (Y[tr] - np.nan_to_num(mu, nan=0.0)) ** 2, np.nan), axis=0))
        return sanitize_mu_sd(mu, sd, mtr, Y[tr])

    MU_SD = {"CERES": prep_mu_sd(Ycer), "Chronos": prep_mu_sd(Ychr)}
    mu_cer, _ = MU_SD["CERES"]
    mu_chr, _ = MU_SD["Chronos"]
    log("sanitize 后 mu/sd 全有限: %s" % all(np.isfinite(MU_SD[k][1]).all() for k in MU_SD))

    # 常量基线：预测 = 训练集 per-gene 均值（对每个细胞同一向量）
    preds["mean_baseline_CERES"] = np.tile(mu_cer, (len(te), 1))
    preds["mean_baseline_Chronos"] = np.tile(mu_chr, (len(te), 1))

    res_cell, res_gene, summary = [], [], []
    for tag, P in preds.items():
        for Ytag, Y in [("CERES", Ycer), ("Chronos", Ychr)]:
            T = Y[te]
            mu, sd = MU_SD[Ytag]
            r_raw, r_shift, r_z, agree = zres_stats(P, T, mte, mu, sd)
            rg, ns = corr_cols(P, T, mte)
            # 逐细胞记录
            for i, c in enumerate(te_cells):
                res_cell.append(dict(model=tag, eval=Ytag, cell=c, n_genes=int(mte[i].sum()),
                                     r_raw=r_raw[i], r_shift=r_shift[i], r_zres=r_z[i]))
            # 逐基因记录（只存 raw；across-cell 对 per-gene 仿射不变）
            for g in range(len(rg)):
                if np.isfinite(rg[g]):
                    res_gene.append(dict(model=tag, eval=Ytag, gene_idx=g, n_cells=int(ns[g]),
                                         r_across_cell=rg[g]))
            summary.append(dict(model=tag, eval=Ytag,
                                n_pc=N_PC, alpha=chosen.get(tag, {}).get("alpha"),
                                cell_raw_median=float(np.nanmedian(r_raw)),
                                cell_shift_median=float(np.nanmedian(r_shift)),
                                cell_zres_median=float(np.nanmedian(r_z)),
                                shift_vs_z_corr=agree,
                                gene_across_cell_median=float(np.nanmedian(rg)),
                                gene_across_cell_frac_pos=float(np.nanmean(rg > 0)),
                                n_genes_across=int(np.isfinite(rg).sum())))
            log("%-22s eval=%-8s raw=%+.4f shift=%+.4f z=%+.4f across=%+.4f (shift~z r=%.3f)"
                % (tag, Ytag, np.nanmedian(r_raw), np.nanmedian(r_shift),
                   np.nanmedian(r_z), np.nanmedian(rg), agree))

    import csv
    with open(R("results", "per_cell_metrics.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res_cell[0].keys())); w.writeheader(); w.writerows(res_cell)
    with open(R("results", "per_gene_across_cell.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res_gene[0].keys())); w.writeheader(); w.writerows(res_gene)
    with open(R("results", "model_metric_summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)

    # --- 2×2（训练靶标 × 评价参考），raw 与 zres 两口径 ---
    def cell2x2(M, Y, Ytag):
        T = Y[te]
        mu, sd = MU_SD[Ytag]
        r_raw, r_shift, r_z, _ = zres_stats(M, T, mte, mu, sd)
        return (float(np.nanmedian(r_raw)), float(np.nanmedian(r_shift)),
                float(np.nanmedian(r_z)))

    grid = {}
    pairs = {"Ridge": ("ridge_CERES", "ridge_Chronos"), "Frozen": ("frozen_CERES", "frozen_Chronos"),
             "Permuted": ("ridge_perm", "ridge_perm")}
    boot_out = {}
    for name, (a_cer, a_chr) in pairs.items():
        cc_raw, cc_s, cc_z = cell2x2(preds[a_cer], Ycer, "CERES")
        ch_raw, ch_s, ch_z = cell2x2(preds[a_cer], Ychr, "Chronos")
        hc_raw, hc_s, hc_z = cell2x2(preds[a_chr], Ycer, "CERES")
        hh_raw, hh_s, hh_z = cell2x2(preds[a_chr], Ychr, "Chronos")
        grid[name] = dict(CC_raw=cc_raw, CH_raw=ch_raw, HC_raw=hc_raw, HH_raw=hh_raw,
                          CC_shift=cc_s, CH_shift=ch_s, HC_shift=hc_s, HH_shift=hh_s,
                          CC_zres=cc_z, CH_zres=ch_z, HC_zres=hc_z, HH_zres=hh_z,
                          interaction_raw=(cc_raw - ch_raw) - (hc_raw - hh_raw),
                          interaction_shift=(cc_s - ch_s) - (hc_s - hh_s),
                          interaction_zres=(cc_z - ch_z) - (hc_z - hh_z))
        log("%s 2x2 raw: CC=%+.4f CH=%+.4f HC=%+.4f HH=%+.4f I=%+.4f | shift I=%+.4f | zres I=%+.4f"
            % (name, cc_raw, ch_raw, hc_raw, hh_raw, grid[name]["interaction_raw"],
               grid[name]["interaction_shift"], grid[name]["interaction_zres"]))

    # --- bootstrap：per-cell residualized 指标 + across-cell 指标 ---
    def cells_of(tag, ev, key):
        return ([r[key] for r in res_cell if r["model"] == tag and r["eval"] == ev],
                [r["cell"] for r in res_cell if r["model"] == tag and r["eval"] == ev])

    for key in ["r_shift", "r_zres"]:
        for tag, ev, nm in [("ridge_CERES", "CERES", "ridge_CERES"), ("ridge_perm", "CERES", "perm"),
                            ("mean_baseline_CERES", "CERES", "mean_baseline"),
                            ("frozen_CERES", "CERES", "frozen_CERES"),
                            ("frozen_Chronos", "Chronos", "frozen_Chronos"),
                            ("ridge_Chronos", "Chronos", "ridge_Chronos")]:
            v, cl = cells_of(tag, ev, key)
            boot_out["%s_%s" % (key, nm)] = boot_median_ci(v, cl)
    for tag in ["ridge_CERES", "frozen_CERES", "ridge_perm", "mean_baseline_CERES"]:
        v = [r["r_across_cell"] for r in res_gene if r["model"] == tag and r["eval"] == "CERES"]
        boot_out["across_cell_" + tag] = boot_median_ci(v, ["g%d" % i for i in range(len(v))])
    json.dump({"grid": grid, "chosen_alpha": chosen,
               "bootstrap": {k: {"point": v[0], "ci95": [v[1], v[2]],
                                 "excludes_zero": v[3]} for k, v in boot_out.items()}},
              open(R("results", "positive_control_summary.json"), "w"), indent=2)
    log("完成，总耗时 %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
