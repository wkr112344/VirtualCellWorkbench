"""核验：全基因残差化 vs 基因子集残差化，训练×评价交互是否都等于 −0.00047。

口径完全照 `GigaScience_Reproducibility_Package/analysis_code/as_supplied/depmap_background_control.py`：
  · 参考矩阵先按 common_observed_mask 与非缺失值构建观测矩阵 obs (136×17,393)；
  · 四格共用同一掩码 Mt（同一细胞下四个相关系数用同一基因集）；
  · 残差化：预测减自身训练臂的逐基因均值，参考减其训练集逐基因均值；
  · 逐细胞交互 I(i)=[CC−CH]−[HC−HH] → 种子内取中位数 → 五种子取均值。

输出：全基因与各"子集"定义下的残差值对照。
"""
import numpy as np
from pathlib import Path

D = Path('depmap_analysis_ready_20260929')
Yce = np.load(D/'reference_test/Y_CERES_test136.npy').astype(float)
Ych = np.load(D/'reference_test/Y_Chronos_test136.npy').astype(float)
mask = np.load(D/'reference_test/common_observed_mask_test136.npy').astype(bool)
muc = np.load(D/'train_gene_mean/train_genemean_Y_CERES.npy').astype(float)
muh = np.load(D/'train_gene_mean/train_genemean_Y_Chronos.npy').astype(float)

obs = np.isfinite(Yce) & np.isfinite(Ych) & mask
frac = obs.mean(axis=0)


def masked_corr(A, B, M):
    A2 = np.where(M, A, np.nan); B2 = np.where(M, B, np.nan)
    ma = np.nanmean(A2, 1, keepdims=True); mb = np.nanmean(B2, 1, keepdims=True)
    da = np.where(M, A - ma, 0.0); db = np.where(M, B - mb, 0.0)
    num = (da*db).sum(1); den = np.sqrt((da*da).sum(1)*(db*db).sum(1)); n = M.sum(1)
    return np.divide(num, den, out=np.full_like(num, np.nan), where=(den > 0) & (n >= 3))


def interaction_residualized(M):
    res = np.empty((5, 136))
    for s in range(5):
        Pc = np.load(D/f'predictions_test/seed{s}_pred_ceres_arm_test.npy').astype(float)
        Ph = np.load(D/f'predictions_test/seed{s}_pred_chronos_arm_test.npy').astype(float)
        res[s] = (masked_corr(Pc-muc, Yce-muc, M) - masked_corr(Pc-muc, Ych-muh, M)
                  - masked_corr(Ph-muh, Yce-muc, M) + masked_corr(Ph-muh, Ych-muh, M))
    per_seed = [float(np.nanmedian(res[s])) for s in range(5)]
    return float(np.mean(per_seed)), per_seed


if __name__ == '__main__':
    print('观测结构：全部 136 细胞都观测到的基因 %d 个；存在缺失的基因 %d 个（占 %.1f%%）'
          % ((frac == 1).sum(), (frac < 1).sum(), 100*(frac < 1).mean()))
    print()
    print('%-44s %-14s %s' % ('口径', '残差交互', '逐种子中位数'))
    m, ps = interaction_residualized(obs)
    print('%-44s %+.8f   %s' % ('全部基因（原始口径，四格共用掩码）', m, np.round(ps, 5).tolist()))
    for th, name in [(1.0, '仅保留 136 细胞全观测的基因'), (0.99, '观测率 ≥99%'), (0.95, '观测率 ≥95%'),
                     (0.90, '观测率 ≥90%'), (0.50, '观测率 ≥50%')]:
        sel = obs & (frac >= th)[None, :]
        m2, ps2 = interaction_residualized(sel)
        print('%-44s %+.8f   %s   (基因数 %d)' % (name, m2, np.round(ps2, 5).tolist(), int((frac >= th).sum())))
    print()
    print('注：17,393 个基因中只有 %d 个（%.1f%%）存在细胞级缺失，'
          '因此任何低于 100%% 的观测率阈值都等价于"全部基因"，数值必然相同。'
          % ((frac < 1).sum(), 100*(frac < 1).mean()))
