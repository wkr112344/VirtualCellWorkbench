# 补充材料清单（GigaScience 投稿清理版 v5）

正文：**《已发布处理后参考产品下的基准稳健性：LINCS L1000 与 DepMap 的可复现分析》**

本投稿包只保留与当前正文主张直接对应的文件。早期 LayerDiag、CPI、癌/非癌分层、遗传侧、旧 baseline/threshold sweep 等探索性材料已从本地投稿包移除，以避免 reviewer 将历史分析误认为当前论文证据链的一部分。冻结 Zenodo version of record 仍负责保存与正文数值对应的完整可复核工件。

## 核心文件

| 编号 | 文件 | 当前用途 |
|---|---|---|
| S13 | `NoteS13_method_and_prior_details.md` | fixed-output/2×2 统计口径、训练暴露边界与 dcic2021 评分对象溯源；S13.7 为逐跳 checksum 链 |
| S14 | `NoteS14_robustness_and_readings.md` | LINCS 对齐/变换检查、排名/shortlist 稳健性、解释边界与 matched 2×2 补充说明 |
| S19 | `TableS19_drug_identity_audit.*` | 药物身份与 ECFP4 重复审计 |
| S20 | `TableS20_绝对topk_筛查排名.csv` | 1,399 药物的跨 reference 性能/疾病签名排名底表 |
| S21 | `TableS21_L1000_provenance_audit.md` | 12 项 L1000 研究的 targeted reporting/provenance audit |
| S21a | `TableS21a_representative_object_tracking_cases.md` | 代表性 score-to-object 追踪案例 |
| S22 | `TableS22_严格双侧未见_逐药.csv` | strict-unseen / fingerprint-clean 逐药结果 |
| S23 | `TableS23_DepMap_fixed_output_replication.*` + `FigS5_DepMap_fixed_output_replication.*` | DepMap CERES/Chronos fixed-output cross-ecosystem stress test |
| S24 | `TableS24_seven_predictor_fixed_output.csv` | 七个已发表预测器固定输出复评的逐模型两臂 PCC、ΔPCC 与 bootstrap 95% CI（臂 A = 各预测器自身的 SDST 参考，臂 B = LINCS dcic2021/CD 产品；均在 978 个 LINCS landmark genes 轴上） |
| S25 | `NoteS25_training_evaluation_2x2.md` + `TableS25_*` | 正文 11,275 对 matched training-product × evaluation-product 2×2 与 E_common=15,990 敏感性 |
| S26 | `NoteS26_gene_caliber_and_seven_model_arms.md` | 978 vs 963 基因口径说明（同一基因轴的两种对齐层次，非维度矛盾）+ 七模型复评两臂定义（含正文 Figure 3 图注修正） |

七个已发表预测器的固定输出复评展示于正文 Figure 3 与补充表 S24；逐模型两臂 PCC、ΔPCC 与 bootstrap 95% CI 见 `TableS24_seven_predictor_fixed_output.csv`，两臂定义与口径说明见 `NoteS26_gene_caliber_and_seven_model_arms.md`。对应 prediction matrices、中间结果和复核工件位于冻结 Zenodo version of record，本地投稿包不重复复制大型矩阵。

## 解释边界

1. 本研究估计的是 **processed reference product replacement sensitivity**，不把 LINCS 的差异唯一归因于某个 preprocessing step。
2. matched 2×2 说明 reference contrast 与 training–evaluation target matching 密切相关；完整对称 2×2 目前只覆盖主预测器。
3. DepMap 是跨数据生态的现象复现，不是对 LINCS `+0.2956` 效应量的直接 replication。
4. strict-unseen/fingerprint-clean 只覆盖当前可追溯的适配训练暴露；更早 base checkpoint 的完整训练组成未保留。

## 数据与代码

- 冻结 version of record：Zenodo DOI `10.5281/zenodo.22720615`
- 开发代码：`https://github.com/wkr112344/G2CP-virtual-cell`
- dcic2021 processed reference 的 SHA256、MD5 与构建链：S13.7
- DepMap fixed-output 结果：S23
- LINCS matched 2×2：S25
