"""within-reference top-k stability（作者提出的关键缺口分析）——最终版。

三个对照：
  V1 行级 bootstrap（全部 2,037 药物）——作者方案：固定参考与药物定义，只重采样组成药物分数的评价行
  V2 行级 bootstrap（仅覆盖 ≥3 行的 1,303 药物子集）——让 bootstrap 真正有力的口径
  V3 药物层面 cluster bootstrap——重采样药物整体（考察"换一批药物"的量级），已修正为按药物身份统计
统一基准（供参考）：V0 随机列表期望重合 k²/N（原稿的弱 null）
输出：within_reference_topk_stability.csv（逐 k 三档 + 跨参考观测 + 随机期望 + 是否低于下界）
"""
import numpy as np
import pandas as pd
from pathlib import Path

SRC = Path('gigascience_figwork2/GigaScience_GTEx_checked/source_data/frozen_prediction_matrices/frozen_perrow_three_metrics.csv')
OUT = Path('gigascience_upload_work/within_reference_topk_stability.csv')
KS = [10, 20, 50, 100, 204, 500, 1000]
B = 200
rng = np.random.default_rng(20260930)
COLS = {'beta': 'pearson_A_beta2020', 'dcic': 'pearson_B_dcic2021'}

d = pd.read_csv(SRC)[['drug', 'cell'] + list(COLS.values())].dropna()
g = d.groupby('drug')
per = {k: {dr: gg.to_numpy(float) for dr, gg in g[c]} for k, c in COLS.items()}
n_row = g.size()
all_drugs = np.sort(d.drug.unique())
sub3 = np.sort(n_row[n_row >= 3].index.to_numpy())


def rank_spearman(a, b):
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def run(drug_list, mode):
    """mode='row' 行级重采样；'drug' 药物整体重采样（按药物身份比较）。"""
    N = len(drug_list)
    obs = {s: np.array([per[s][dr].mean() for dr in drug_list]) for s in COLS}
    res = {}
    for s in COLS:
        reps, rho = [], []
        for b in range(B):
            if mode == 'row':
                v = np.empty(N)
                for i, dr in enumerate(drug_list):
                    a = per[s][dr]; n = a.size
                    v[i] = a[rng.integers(0, n, n)].mean() if n > 1 else a[0]
            else:
                pick = rng.integers(0, N, N)          # 抽药物身份
                cnt = np.bincount(pick, minlength=N)
                v = np.array([per[s][drug_list[i]].mean() if cnt[i] else 0.0 for i in range(N)])
            reps.append(v)
            if b < 20:
                rho.append(v)
        per_k = {}
        for k in KS:
            tops = [set(np.argsort(-v, kind='stable')[:k]) for v in reps]
            ov = [len(tops[i] & tops[j]) for i in range(B) for j in range(i + 1, B)]
            per_k[k] = (float(np.mean(ov)), float(np.percentile(ov, 2.5)), float(np.percentile(ov, 97.5)))
        res[s] = {'obs': obs[s], 'per_k': per_k,
                  'rho': [rank_spearman(rho[i], rho[j]) for i in range(len(rho)) for j in range(i + 1, len(rho))]}
    return res


def cross_overlap(res):
    o1, o2 = res['beta']['obs'], res['dcic']['obs']
    return {k: len(set(np.argsort(-o1, kind='stable')[:k]) & set(np.argsort(-o2, kind='stable')[:k])) for k in KS}


print('=' * 92)
print('V1  行级 bootstrap，全部药物 N=%d' % len(all_drugs))
r1 = run(all_drugs, 'row'); c1 = cross_overlap(r1)
print('%-6s %-22s %-22s %-12s %-10s' % ('k', 'within-β [2.5%,97.5%]', 'within-dcic', 'cross-ref', 'random'))
for k in KS:
    b, c = r1['beta']['per_k'][k], r1['dcic']['per_k'][k]
    print('%-6d %-22s %-22s %-12s %-10.2f' % (k, '%.1f [%.0f, %.0f]' % b, '%.1f [%.0f, %.0f]' % c,
                                              '%d (%.0f%%)' % (c1[k], 100 * c1[k] / k), k * k / len(all_drugs)))
print('replicate 间排名 Spearman：β %.4f，dcic %.4f ｜ 观测跨参考 %.4f' % (
    np.mean(r1['beta']['rho']), np.mean(r1['dcic']['rho']), rank_spearman(r1['beta']['obs'], r1['dcic']['obs'])))

print()
print('=' * 92)
print('V2  行级 bootstrap，仅覆盖 ≥3 行的子集 N=%d' % len(sub3))
r2 = run(sub3, 'row'); c2 = cross_overlap(r2)
print('%-6s %-22s %-22s %-12s %-12s' % ('k', 'within-β [2.5%,97.5%]', 'within-dcic', 'cross-ref', '是否低于下界'))
for k in KS:
    b, c = r2['beta']['per_k'][k], r2['dcic']['per_k'][k]
    flag = '是' if c2[k] < max(b[1], c[1]) else '否（在波动内）'
    print('%-6d %-22s %-22s %-12s %-12s' % (k, '%.1f [%.0f, %.0f]' % b, '%.1f [%.0f, %.0f]' % c,
                                            '%d (%.0f%%)' % (c2[k], 100 * c2[k] / k), flag))
print('replicate 间排名 Spearman：β %.4f，dcic %.4f ｜ 观测跨参考 %.4f' % (
    np.mean(r2['beta']['rho']), np.mean(r2['dcic']['rho']), rank_spearman(r2['beta']['obs'], r2['dcic']['obs'])))

print()
print('=' * 92)
print('V3  药物层面 cluster bootstrap（全药物）')
r3 = run(all_drugs, 'drug'); c3 = cross_overlap(r3)
for k in [10, 20, 50, 100, 204]:
    b, c = r3['beta']['per_k'][k], r3['dcic']['per_k'][k]
    print('k=%-5d within-β %.1f [%.0f, %.0f] | within-dcic %.1f [%.0f, %.0f] | cross-ref %d' % (
        k, b[0], b[1], b[2], c[0], c[1], c[2], c3[k]))

rows = []
for k in KS:
    b1, c1_ = r1['beta']['per_k'][k], r1['dcic']['per_k'][k]
    b2, c2_ = r2['beta']['per_k'][k], r2['dcic']['per_k'][k]
    b3 = r3['beta']['per_k'][k]
    rows.append({
        'k': k,
        'V1_within_beta_mean': b1[0], 'V1_within_beta_lo': b1[1], 'V1_within_beta_hi': b1[2],
        'V1_within_dcic_mean': c1_[0], 'V1_within_dcic_lo': c1_[1], 'V1_within_dcic_hi': c1_[2],
        'V1_cross_reference': c1[k], 'V1_observed_below_lo': bool(c1[k] < max(b1[1], c1_[1])),
        'V2_within_beta_mean': b2[0], 'V2_within_beta_lo': b2[1], 'V2_within_beta_hi': b2[2],
        'V2_within_dcic_mean': c2_[0], 'V2_within_dcic_lo': c2_[1], 'V2_within_dcic_hi': c2_[2],
        'V2_cross_reference': c2[k], 'V2_observed_below_lo': bool(c2[k] < max(b2[1], c2_[1])),
        'V3_within_beta_mean': b3[0], 'V3_cross_reference': c3[k],
        'random_expected_overlap': k * k / len(all_drugs),
        'V2_n_drugs': len(sub3), 'V1_n_drugs': len(all_drugs)})
pd.DataFrame(rows).to_csv(OUT, index=False)
print()
print('输出:', OUT)
