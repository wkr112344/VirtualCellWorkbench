"""全结论数值验算：把正文里的结论数字与归档表逐一比对。

A 部分：头号结论逐条比对（容差 5e-4，区间端点 1e-3）
B 部分：覆盖性检查——正文每个 ≥3 位有效数字的数，是否能在归档文件中找到出处
"""
import json, re, os
import pandas as pd
import numpy as np

PKG = r"C:/Users/wkr20/WorkBuddy/Claw/gigascience_upload_work/GigaScience_Submission_Upload"
SUP = os.path.join(PKG, "supplementary")
import docx
D = docx.Document(os.path.join(PKG, "manuscript/Main_Text_CN_GigaScience_revised.docx"))
TXT = "\n".join(p.text for p in D.paragraphs)

S = {f: pd.read_csv(os.path.join(SUP, f)) for f in os.listdir(SUP) if f.endswith(".csv")}
J = json.load(open(os.path.join(SUP, "analysis_summary.json"), encoding="utf-8"))

def v(x, d=6):
    return round(float(x), d)

checks = []          # (说明, 正文值, 归档值, 容差, 出处)
def C(name, claim, src, tol=5e-4, where=""):
    checks.append((name, claim, v(src, 6), tol, where))

# ---------------- LINCS ----------------
s1 = S["S1_weighting_and_metrics.csv"]
g = lambda **kw: s1.loc[s1[list(kw)].eq(pd.Series(kw)).all(axis=1)].iloc[0]
C("主 Δ（药物等权）", 0.2956, g(metric="pearson", weighting="drug").delta, 5e-4, "S1")
C("主 Δ 区间下", 0.2918, g(metric="pearson", weighting="drug").ci_low, 1e-3, "S1")
C("主 Δ 区间上", 0.2993, g(metric="pearson", weighting="drug").ci_high, 1e-3, "S1")
C("beta2020 均值", 0.3680, g(metric="pearson", weighting="drug").beta_mean, 5e-4, "S1")
C("dcic2021 均值", 0.0724, g(metric="pearson", weighting="drug").dcic_mean, 5e-4, "S1")
C("行等权 Δ", 0.2996, g(metric="pearson", weighting="row").delta, 5e-4, "S1")
C("细胞系等权 Δ", 0.2770, g(metric="pearson", weighting="cell").delta, 5e-4, "S1")
C("细胞系等权 Δ 下", 0.2541, g(metric="pearson", weighting="cell").ci_low, 1e-3, "S1")
C("细胞系等权 Δ 上", 0.3023, g(metric="pearson", weighting="cell").ci_high, 1e-3, "S1")
C("细胞系 Δ 最小值", 0.1554, g(metric="pearson", weighting="cell").min_delta, 1e-3, "S1")
C("细胞系 Δ 最大值", 0.6281, g(metric="pearson", weighting="cell").max_delta, 1e-3, "S1")
C("11,275 评价对", 11275, J["lincs"]["rows"], 0, "analysis_summary")
C("63 细胞系", 63, J["lincs"]["cells"], 0, "analysis_summary")
C("2,037 留出药物", 2037, J["lincs"]["drugs"], 0, "analysis_summary")
C("排序 Spearman", 0.5675, J["lincs"]["rank_rho"], 5e-4, "analysis_summary")
s3 = S["S3_topk_performance_ranking.csv"].set_index("k")
C("k=10 保留率", 0.30, s3.loc[10, "retention"], 5e-4, "S3")
C("k=10 重合数", 3, s3.loc[10, "overlap"], 0, "S3")
C("k=10 Jaccard", 0.1765, s3.loc[10, "jaccard"], 5e-4, "S3")
C("k=10 随机期望", 0.049, s3.loc[10, "chance_expected_overlap"], 1e-3, "S3")
C("k=20 保留率", 0.65, s3.loc[20, "retention"], 5e-4, "S3")
C("k=50 保留率", 0.60, s3.loc[50, "retention"], 5e-4, "S3")
C("k=204 保留率", 0.5980, s3.loc[204, "retention"], 5e-4, "S3")
C("k=204 重合数", 122, s3.loc[204, "overlap"], 0, "S3")
C("k=204 Jaccard", 0.4266, s3.loc[204, "jaccard"], 5e-4, "S3")
C("k=204 随机期望", 20.43, s3.loc[204, "chance_expected_overlap"], 0.02, "S3")
C("k=500 保留率", 0.6040, s3.loc[500, "retention"], 5e-4, "S3")
C("k=1000 保留率", 0.7100, s3.loc[1000, "retention"], 5e-4, "S3")
s2 = S["S2_leave_one_cell_out.csv"]
C("逐一排除下界", 0.2924, s2.drug_weighted_delta.min(), 5e-4, "S2")
C("逐一排除上界", 0.3028, s2.drug_weighted_delta.max(), 5e-4, "S2")
# ---------------- DepMap ----------------
s6 = S["S6_depmap_interaction_bootstrap_summary.csv"].set_index("analysis")
C("原始交互", 0.05632, s6.loc["raw", "estimate"], 5e-5, "S6")
C("原始交互区间下", 0.05589, s6.loc["raw", "ci_low"], 1e-3, "S6")
C("原始交互区间上", 0.05703, s6.loc["raw", "ci_high"], 1e-3, "S6")
C("残差交互", -0.00047, s6.loc["residualized", "estimate"], 5e-5, "S6")
C("残差交互区间下", -0.00471, s6.loc["residualized", "ci_low"], 1e-3, "S6")
C("残差交互区间上", 0.00462, s6.loc["residualized", "ci_high"], 1e-3, "S6")
C("种子残差最小", -0.00413, s6.loc["residualized", "seed_min"], 1e-3, "S6")
C("种子残差最大", 0.00132, s6.loc["residualized", "seed_max"], 1e-3, "S6")
gw = J["gene_wise_median_range"]
C("基因级最小", -0.0054, gw[0], 5e-5, "analysis_summary")
C("基因级最大", 0.0059, gw[1], 5e-5, "analysis_summary")
s7 = S["S7_deepdep_common_grid_summary.csv"]
dd_cell = s7[(s7["axis"].str.contains("cell-wise")) & (s7.metric == "pearson")]
dd_gene = s7[(s7["axis"].str.contains("gene-wise")) & (s7.metric == "pearson")]
C("DeepDEP CERES 细胞级", 0.8690, dd_cell[dd_cell["reference"] == "CERES"]["median"].iloc[0], 5e-4, "S7")
C("DeepDEP Chronos 细胞级", 0.8370, dd_cell[dd_cell["reference"] == "Chronos"]["median"].iloc[0], 5e-4, "S7")
C("DeepDEP CERES 基因级", 0.0422, dd_gene[dd_gene["reference"] == "CERES"]["median"].iloc[0], 5e-4, "S7")
C("DeepDEP Chronos 基因级", 0.1040, dd_gene[dd_gene["reference"] == "Chronos"]["median"].iloc[0], 5e-4, "S7")
s13 = S["S13_depmap_positive_control_model_metrics.csv"]
r13 = s13[(s13["model"] == "ridge_CERES") & (s13["eval"] == "CERES")].iloc[0]
C("Ridge CERES residualized", 0.2164, r13.cell_shift_median, 5e-4, "S13")
C("Ridge CERES z-residualized", 0.1432, r13.cell_zres_median, 5e-4, "S13")
C("Ridge CERES 基因级", 0.1233, r13.gene_across_cell_median, 5e-4, "S13")
C("Ridge CERES 正相关比例", 86.7, r13.gene_across_cell_frac_pos * 100, 0.1, "S13")
r13b = s13[(s13["model"] == "ridge_Chronos") & (s13["eval"] == "Chronos")].iloc[0]
C("Ridge Chronos residualized", 0.2771, r13b.cell_shift_median, 5e-4, "S13")
C("Ridge Chronos z-residualized", 0.1833, r13b.cell_zres_median, 5e-4, "S13")
C("Ridge Chronos 基因级", 0.1621, r13b.gene_across_cell_median, 5e-4, "S13")
C("Ridge Chronos 正相关比例", 92.9, r13b.gene_across_cell_frac_pos * 100, 0.1, "S13")
s15 = S["S15_depmap_positive_control_interaction_bootstrap.csv"]
q = lambda m, mm: s15[(s15["metric"] == m) & (s15["model"] == mm)].iloc[0]
C("Ridge 交互 raw", 0.0620, q("raw", "Ridge").point, 5e-4, "S15")
C("Ridge 交互 residualized", 0.3057, q("shift", "Ridge").point, 5e-4, "S15")
C("冻结交互 raw", 0.0573, q("raw", "Frozen").point, 5e-4, "S15")
C("冻结交互 residualized", -0.0003, q("shift", "Frozen").point, 5e-5, "S15")
C("冻结 residualized 区间下", -0.0172, q("shift", "Frozen").ci95_low, 1e-3, "S15")
C("冻结 residualized 区间上", 0.0167, q("shift", "Frozen").ci95_high, 1e-3, "S15")
C("Ridge 交互 residualized 区间下", 0.2731, q("shift", "Ridge").ci95_low, 1e-3, "S15")
C("Ridge 交互 residualized 区间上", 0.3408, q("shift", "Ridge").ci95_high, 1e-3, "S15")
# ---------------- GTEx ----------------
C("GTEx 交互 Pearson", 0.07096, next(x for x in J["gtex"]["statistics"] if x["statistic"] == "interaction_pearson")["estimate"], 5e-5, "analysis_summary")
C("GTEx 交互区间下", 0.07050, next(x for x in J["gtex"]["statistics"] if x["statistic"] == "interaction_pearson")["ci_low"], 1e-3, "analysis_summary")
C("GTEx 交互区间上", 0.07155, next(x for x in J["gtex"]["statistics"] if x["statistic"] == "interaction_pearson")["ci_high"], 1e-3, "analysis_summary")
C("GTEx 交互 Spearman", 0.06187, next(x for x in J["gtex"]["statistics"] if x["statistic"] == "interaction_spearman")["estimate"], 5e-5, "analysis_summary")
qm = J["gtex"]["pearson_quadrant_medians"]
C("GTEx rRR", 0.9853, qm["rRR"], 5e-4, "analysis_summary")
C("GTEx rRE", 0.9491, qm["rRE"], 5e-4, "analysis_summary")
C("GTEx rER", 0.9488, qm["rER"], 5e-4, "analysis_summary")
C("GTEx rEE", 0.9854, qm["rEE"], 5e-4, "analysis_summary")
C("GTEx 样本数", 3963, J["gtex"]["samples"], 0, "analysis_summary")
C("GTEx 供者数", 189, J["gtex"]["donors"], 0, "analysis_summary")
C("GTEx 组织数", 68, J["gtex"]["tissues"], 0, "analysis_summary")
C("GTEx 组织下界", 0.05972, J["gtex"]["tissue_range"][0], 5e-5, "analysis_summary")
C("GTEx 组织上界", 0.10828, J["gtex"]["tissue_range"][1], 5e-5, "analysis_summary")
C("GTEx 负优势样本", 6, J["gtex"]["sample_adv_E_negative"], 0, "analysis_summary")
C("GTEx 低覆盖组织", 18, J["gtex"]["low_coverage_tissues"], 0, "analysis_summary")

print("=" * 78)
print("A 部分：头号结论逐条比对（正文值 vs 归档值）")
print("=" * 78)
bad = 0
for name, claim, src, tol, where in checks:
    ok = abs(float(claim) - float(src)) <= tol
    if not ok:
        bad += 1
    print('%-26s 正文 %-10s 归档 %-10s %s  [%s]' % (name, claim, src, 'OK' if ok else '✗ 不一致', where))
print('-' * 78)
print('比对 %d 条，不一致 %d 条' % (len(checks), bad))

# ---------------- B 部分：覆盖性 ----------------
print()
print("=" * 78)
print("B 部分：正文所有高精度数字的出处覆盖")
print("=" * 78)
corpus = {}
for f, d in S.items():
    corpus[f] = d.to_csv(index=False)
corpus["analysis_summary.json"] = json.dumps(J, ensure_ascii=False)
nums = set(re.findall(r"\d+\.\d{3,6}", TXT))
missing = []
for x in sorted(nums, key=float):
    hit = [f for f, c in corpus.items() if x in c]
    if not hit:
        missing.append(x)
print('正文高精度数字（≥3 位小数）共 %d 个；在归档表/摘要中能找到原值的 %d 个' %
      (len(nums), len(nums) - len(missing)))
if missing:
    print('未能直接匹配的（多为由归档值换算、区间上限/下限或跨表派生）：')
    print('  ' + ', '.join(missing))
