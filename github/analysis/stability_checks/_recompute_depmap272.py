"""DepMap 272 细胞共同网格：两参考一致性、冻结预测均值、配对 Δ 及其 bootstrap CI 的重算。"""
import pandas as pd, numpy as np, json
V = 'gigascience_v5/图源数据/'
cer = pd.read_csv(V + 'DepMap_CERES_272x1244.csv', index_col=0)
chr_ = pd.read_csv(V + 'DepMap_Chronos_272x1244.csv', index_col=0)
pred = pd.read_csv('depmap_raw_artifacts_20260929/first_predictor_DeepDEP/DeepDEP_predictor_strict.csv', index_col=0)

def cw(A, B):
    out = {}
    for c in A.index:
        if c not in B.index: continue
        a, b = A.loc[c].to_numpy(float), B.loc[c].to_numpy(float)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 10 or np.std(a[m]) == 0 or np.std(b[m]) == 0: continue
        out[c] = float(np.corrcoef(a[m], b[m])[0, 1])
    return pd.Series(out)

cons, pc, pcx = cw(cer, chr_), cw(pred, cer), cw(pred, chr_)
common = cons.index.intersection(pc.index).intersection(pcx.index)
v = (pc[common] - pcx[common]).to_numpy()
bs = v[np.random.default_rng(74).integers(0, len(v), (5000, len(v)))].mean(1)
arch = json.load(open('gigascience_upload_work/GigaScience_Reproducibility_Package/analysis_code/as_supplied/analysis_summary.json', encoding='utf-8'))['depmap_delta']
print('一致性 mean/median      %.6f / %.6f   (正文 0.9431 / 0.9467)' % (cons[common].mean(), cons[common].median()))
print('冻结预测 vs CERES/Chronos %.6f / %.6f (正文 0.8595 / 0.8226)' % (pc[common].mean(), pcx[common].mean()))
print('配对 delta mean          %.6f           (正文 0.0369)' % v.mean())
print('delta CI (seed=74)       [%.6f, %.6f] (归档 [%.6f, %.6f])' % (np.percentile(bs,2.5), np.percentile(bs,97.5), arch[1], arch[2]))
print('delta>0 细胞 %d/%d  范围 %.6f ~ %.6f' % ((v>0).sum(), len(v), v.min(), v.max()))
