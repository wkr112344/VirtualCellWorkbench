# 补充材料 S14：LINCS 跨参考产品稳健性与对齐检查

本文件仅保留当前 GigaScience 稿件直接使用的跨 reference 检查。旧版中与 LayerDiag、癌/非癌分层、CPI、遗传扰动或其他已退出当前主线的内容已从投稿包移除。

## S14.1 主固定输出比较

在同一 beta-trained frozen prediction、同一 64 个共同细胞系、同一 11,275 个共同 `(cell line, drug)` 对和同一评分协议下，MODZ-based level5beta2020 的逐样本 PCC 为 0.3680，CD-based dcic2021 为 0.0724，`Δ_ref = +0.2956`。该差值是完整 data-product replacement contrast，不解释为某个单一预处理步骤的因果效应。

## S14.2 训练暴露检查

2,037 个受控排名候选中，192 个药物在两臂当前可追溯的适配训练数据中均未出现，其 `ΔPCC = +0.2289`。进一步排除药物 ID 已见或任一训练侧存在完全相同 ECFP4 指纹后，55 个 fingerprint-clean 药物的 `ΔPCC = +0.3118`。这些结果说明主差异不能由当前可追溯适配训练中的直接药物暴露或完全相同指纹单独解释；更早 base checkpoint 的训练组成不完整，因此不作更强外推。逐药底表见 S22。

## S14.3 药物排序与疾病签名短名单

对 2,037 个候选药物的性能排序，跨产品 Spearman = 0.5675，arm A top-10% 在 arm B 的保留率为 59.80%，Jaccard = 0.4266；但 arm A top-10% 中有 12.25% 在 arm B 掉出 top-50%。在 1,399 个共同可评药物上，三条疾病方向性签名的跨产品 Spearman 为 0.498–0.540，top-10% 保留率为 33.6%–42.9%，说明 global rank stability 与 top-k/shortlist stability 不是同一层面的量。排名底表见 S20。

## S14.4 七个已发表预测器的固定输出复评

正文 Figure 3 固定 CIGER、DeepCE、MultiDCP、PertDiT、PRnet、TranSiGen 和 XPert 的既有 prediction matrices，只替换对应 Level-5 reference。七个预测器的 paired `ΔPCC` 方向一致，中位数 0.3993，范围 0.2599–0.6450，各自 bootstrap 95% CI 下界均高于 0。逐模型可复核工件位于正文对应的冻结 Zenodo version of record；当前本地投稿包未重复复制大型中间矩阵。

## S14.5 解释边界

LINCS 的跨产品差异可能同时包含 signature estimation、normalization、replicate aggregation、quality control、dynamic range、measurement repeatability 与 processed-construct definition 等因素。当前设计没有 within-product test-retest 或 technical-replicate noise baseline，因此不进一步拆分这些来源。

## S14.6 数据源层：共同签名的内容差异与常见变换检查

两套 Level-5 产品的 signature-key overlap 为 99.61%，但正文 11,275 个共同可评响应在 978 个 LINCS landmark genes 上的逐行 Pearson 相关中位数仅为 0.4172。为排查键对齐和简单数值变换造成的假象，我们在共享 signature-id 的原始发布矩阵上重新对齐条件名和基因，得到逐行 Pearson 中位数 0.4262、cosine 中位数 0.4224，与主分析同量级。逐行 z-score、正缩放、逐基因中心化/z 化和逐基因标定均未消除主要结构差异。

这些检查只能说明差异不是由所检验的简单对齐/变换问题单独造成；MODZ 与 Characteristic Direction 的生成谱系包含多项差异，本文不将 `Δ_ref` 唯一归因于 normalization、replicate aggregation、quality control 或任何单一算法步骤。

## S14.7 基因面板敏感性

将评分基因面板切换为 978 个 LINCS landmark genes 后，绝对 `ΔPCC` 增大，同时 top-10% 重合率提高。该结果说明绝对性能差异与排名稳定性可以呈现不同方向的变化，因此正文同时报告 score shift 与 rank/top-k 指标。

## S14.8 共同集合敏感性

正文 2×2 使用与主结果完全对齐的 11,275 对；另在 `E_common = 15,990` 的更大共同集合上重复，interaction 方向一致但绝对数值不同。两套结果分别报告，避免把集合变化与 reference replacement 混在同一个效应量中。完整结果见 S25。

## S14.9 其他扰动转录组资源的可对齐性检查

我们评估了 GSE70138、CPJUMP1 等资源是否能够提供与 LINCS 主分析相同任务、同终点、且可在同一评分单元上构造第二套已发布 processed reference product。现有材料未提供可完成同类受控复制的组合，因此本文不把这些资源作为等价 replication；跨数据生态检验转向 DepMap CRISPR gene-effect products（S23）。

## S14.10 matched 2×2 对主固定输出结果的解释

与 11,275 对主 fixed-output comparison 完全对齐的 2×2 结果为：beta-trained 的 reference contrast `+0.2956`，dcic-trained 的 reference contrast `+0.0171`，difference-in-differences interaction `+0.2785 [0.2741, 0.2828]`。两种训练状态下的 reference contrast 明显不同，因此正文将结果表述为 reference-product sensitivity 与 training–evaluation target matching 密切相关，而不是固定不变的“产品效应”。完整四格和更大共同集合敏感性见 S25。

## S14.11 结果解释与报告口径

本稿按 comparative benchmark robustness / reproducibility study 解释结果。LINCS 的整体药物排序保留中等结构，但 top-k 与 disease-signature shortlist 更敏感；DepMap 中 CERES/Chronos 全局高度一致，固定 DeepDEP prediction 的 performance reading 与 top-k recovery 仍有较小但可测量的变化。DepMap 因而提供跨数据生态的现象复现，而不是对 LINCS 具体效应量的直接复制。
