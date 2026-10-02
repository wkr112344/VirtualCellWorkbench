# 结论 → 脚本 → 输入 → 期望输出

| 论文结论 | 脚本 | 输入 | 期望输出（校验值） |
|---|---|---|---|
| LINCS 主 Δ 0.2956（95% CI 0.2918–0.2993） | `github/analysis/stability_checks/_audit_numbers.py` 或 as-supplied | `zenodo/source_data/frozen_prediction_matrices/frozen_perrow_three_metrics.csv` | S1 表同值 |
| 三种指标跨矩阵同向 | 同上 | 同上 | S1 表（Pearson/Spearman/余弦 × 药物/细胞系） |
| 训练状态反转 Dβ=+0.1975、Ddcic=−0.0810、I=0.2785 | `_patch_A_B.py` | `supplementary/frozen_multimetric_perrow_scores.csv` | 表 3、S20 |
| crossed 2×2 反转保持 2000/2000 | `_patch_A_B.py` | 同上 | S20（direction_probability=1.0） |
| 主 Δ 的 crossed 区间 0.2855–0.3128 | `_crossed_bootstrap_and_disease.py` | frozen per-row 分数 | S18 |
| 排名 Spearman 0.5675 / Kendall τ 0.4110 | `_audit_numbers.py` + `_within_ref_stability.py` | per-drug 排名 | S21 |
| top-10 跨矩阵 3/10；同参考 7.93 / 7.04 | `_within_ref_stability.py`、`_patch_A_B.py` | per-row 分数（B=2000 独立 pair） | S17、S17b、S17c |
| top-k 与随机期望 | `_audit_numbers.py` | 药物数 2,037 | S22 |
| 疾病榜单边界 0.14%–10.85%、带 11–13 | `_crossed_bootstrap_and_disease.py` | `图源数据/TableS_disease_perdrug_scored.csv` | S19 |
| DepMap 交互 0.05632 → −0.00047 | `as_supplied/depmap_background_control.py` | `source_data/depmap_analysis_ready_*` | 表 4、S6 |
| 基因均值基线 ≈0.90–0.94 | 同上 | 同上 | S4 |
| 阳性对照 0.2164 / 0.2771 / 0.1233 / 0.1621 | `positive_control_expression_ridge.py` | `depmap_positive_control`（+ 见 DATA_SOURCES 的 X_TPM 重建） | S13–S15 |
| DepMap 272 细胞一致性 0.9431/0.9467 与预测均值 0.8595/0.8226 | `_recompute_depmap272.py` | `图源数据/DepMap_{CERES,Chronos}_272x1244.csv` + `DeepDEP_predictor_strict.csv` | S16 |
| GTEx 交互 0.07096 / 0.06187（中位数的中位数） | `gtex_reference_sensitivity/scripts/07_bootstrap.py` | `gtex_reference_sensitivity`（逐样本得分） | S8、S10–S12 |
| 图 1–6 与单面板 | `github/analysis/figure_build/*.py` | `zenodo/source_data/图源数据` | `zenodo/figures/` 同名文件 |

> 校验口径：所有关键数字容差 5×10⁻⁴；区间端点 1×10⁻³。逐一比对脚本见 `EXPECTED_OUTPUTS.json` 与 `github/analysis/stability_checks/_audit_numbers.py`。
