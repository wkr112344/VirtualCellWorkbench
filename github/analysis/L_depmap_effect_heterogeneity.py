"""L_depmap_effect_heterogeneity.py
基于 图源数据/DepMap_percell_delta_eval.csv（272 细胞系 × 27 列），
刻画 DeepDEP fixed-output 评估中 reference-product 效应的异质性：
- delta_pearson / delta_spearman 的分布
- delta 与 product concordance（两参考平均 PCC）的相关
- top-k retention 的 delta 分布
"""
import csv, os
import numpy as np

SRC = 'C:/Users/wkr20/WorkBuddy/Claw/gigascience_v5/图源数据'
RES = 'C:/Users/wkr20/WorkBuddy/Claw/gigascience_v5/results'
os.makedirs(RES, exist_ok=True)

rows = list(csv.DictReader(open(f'{SRC}/DepMap_percell_delta_eval.csv', encoding='utf-8-sig')))
n = len(rows)
def col(name): return np.array([float(r[name]) for r in rows])

dp = col('delta_pearson'); ds = col('delta_spearman')
pc_ceres = col('pearson_vs_ceres'); pc_chronos = col('pearson_vs_chronos')
prod_conc = (pc_ceres + pc_chronos) / 2.0

def desc(x):
    return dict(n=len(x), mean=round(float(np.mean(x)),4), median=round(float(np.median(x)),4),
               std=round(float(np.std(x)),4), min=round(float(np.min(x)),4), max=round(float(np.max(x)),4))

from scipy.stats import pearsonr
corr_dp_prod = pearsonr(dp, prod_conc).statistic
corr_ds_prod = pearsonr(ds, prod_conc).statistic

summary = {
    'n_cell_lines': n,
    'delta_pearson': desc(dp),
    'delta_spearman': desc(ds),
    'corr_delta_pearson_vs_product_concordance': round(float(corr_dp_prod),4),
    'corr_delta_spearman_vs_product_concordance': round(float(corr_ds_prod),4),
    'frac_cells_delta_pearson_gt_0': round(float(np.mean(dp > 0)),4),
    'frac_cells_delta_spearman_gt_0': round(float(np.mean(ds > 0)),4),
}
# per top-k retention delta
for k in [10, 50, 100]:
    rc = col(f'top{k}_retention_vs_ceres'); rk = col(f'top{k}_retention_vs_chronos')
    d = rc - rk
    summary[f'delta_top{k}_retention'] = desc(d)

with open(f'{RES}/L_depmap_effect_modifiers.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['metric', 'value'])
    for k, v in summary.items():
        w.writerow([k, v if not isinstance(v, dict) else str(v)])
print('wrote L_depmap_effect_modifiers.csv')
for k, v in summary.items():
    print('  ', k, '=', v)
