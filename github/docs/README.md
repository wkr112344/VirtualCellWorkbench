# DepMap positive control：expression-PCA → multi-output Ridge

**问题**：在 gene-level background 被 residualize 掉之后，还有没有模型能保留可检测的预测信号？
如果只有"有真实细胞状态信息的模型"能保留，那么"冻结模型 residualize 后归零"就不能被解释成
"residualization 机械地杀死一切信号"。

全部计算在同一条 **Figure 6 划分**上（split.json seed 42：train 636 / val 136 / test 136 细胞），
数据为 DepMap 21Q2：`X_TPM.npy` (908×19177)、`Y_CERES/Y_Chronos.npy` (908×17393)、`common_observed_mask`。

## 对象（4 个 × 2 个评价参考）

| 对象 | 定义 |
|---|---|
| gene-mean baseline | 预测 = 训练集 per-gene 均值（对每个细胞同一向量） |
| **expression-PCA→Ridge（positive control）** | X→训练集标准化→PCA(100)→multi-output Ridge，训练靶标分别 = CERES / Chronos |
| permuted-label Ridge（negative control） | 同上，但训练 Y 的行标签被打乱 |
| frozen G2CP（现模型） | 既有 5-seed 平均预测（`depmap_analysis_ready_20260929`） |

alpha 只用**训练集 5-fold CV** 在 {0.1,1,10,100,1000} 上选（test 全程不参与）；三臂都选中 1000。

## 三个口径（关键定义）

- `raw`：逐细胞跨基因 Pearson(pred, truth)
- `shift`：pred 与 truth **同时减训练集 per-gene 均值**（原设计）。注意：这是**逐基因**常数（c_g），
  不是全局常数，所以 `cov(x−c, y−c) = cov(x,y) − cov(x,c) − cov(y,c) + var(c) ≠ cov(x,y)` ——
  **它会真实改变跨基因相关**，剥掉的正是基因间共享的 essentiality background。
- `z`：在 shift 基础上再除以训练集 per-gene 标准差（额外做尺度归一）
- `across-cell`：每个基因跨 test 细胞的 Pearson（对 per-gene 平移/缩放不变）

## 结果（median，5000× cell bootstrap CI）

| 对象 | raw | shift（残差化） | z（残差化+尺度） | across-cell |
|---|---|---|---|---|
| gene-mean baseline | +0.9338 | **0.0000** | 0.0000 | n/a（恒定预测器） |
| permuted Ridge | +0.9272 | **−0.0001** [−0.0082,+0.0074] | −0.0015 | +0.0015 |
| **expression Ridge (CERES)** | +0.9345 | **+0.2164** [+0.1951,+0.2321] ✓ | **+0.1432** ✓ | **+0.1233** ✓ |
| **expression Ridge (Chronos)** | +0.9434 | **+0.2771** [+0.2568,+0.3063] ✓ | **+0.1833** ✓ | **+0.1621** ✓ |
| frozen G2CP (CERES) | +0.9336 | −0.0012 [−0.0096,+0.0069] ✗ | +0.0013 | +0.0003 |
| frozen G2CP (Chronos) | +0.9389 | −0.0027 ✗ | +0.0040 | +0.0055 |

## 2×2（训练靶标 × 评价参考），interaction = (CC−CH)−(HC−HH)

| 口径 | expression Ridge | frozen G2CP | permuted |
|---|---|---|---|
| raw | +0.0620 [+0.0587,+0.0654] ✓ | +0.0573 [+0.0549,+0.0597] ✓ | 0 |
| **shift（残差化）** | **+0.3057** [+0.2731,+0.3408] ✓ | **−0.0003** [−0.0172,+0.0167] ✗ | 0 |
| z（残差化+尺度） | **+0.2101** [+0.1885,+0.2347] ✓ | +0.0094 ✗ | 0 |

## 结论（审稿人闭环）

- **residualization 不会机械杀掉一切信号**：真预测器残差化后仍 +0.216 / +0.277，CI 排零。
- **流程本身不制造信号**：permuted Ridge ≈ 0，gene-mean baseline 恰好 0。
- **冻结模型的 product-matching interaction 是纯 gene-level background**：raw +0.0573 → 残差化 −0.0003（归零），
  而 expression Ridge 的 interaction 残差化后反而**放大**到 +0.3057。
  → reference-product sensitivity 的两部分被彻底分开：冻结模型那一份来自共享 gene 均值结构，
  Ridge 那一份来自真实细胞状态预测能力。

## 文件

```
scripts/positive_control_expression_ridge.py   主实验（全部指标 + bootstrap）
scripts/bootstrap_interaction.py               interaction 的 5000× cell bootstrap
scripts/make_panelAB_figure.py                 Panel A/B 图
scripts/make_panelC_figure.py                  Panel C（2×2）图
results/model_metric_summary.csv               4 对象 × 2 参考 × 全指标
results/per_cell_metrics.csv                   逐细胞 raw/shift/zres
results/per_gene_across_cell.csv               逐基因 across-cell（17,393 基因）
results/positive_control_summary.json          汇总 + bootstrap
results/interaction_bootstrap.json             interaction CI（三口径 × 三对象）
figures/panel_AB_positive_control.png          Panel A/B
figures/fig_panelC_2x2_heatmap.png             Panel C
```

## 已知边界

- permuted Ridge 的 across-cell 中位 +0.0015，CI [+0.0001,+0.0034] 技术上排零——但比 Ridge 的 +0.1233
  小两个数量级，属数值噪声级别，报数时如实列出。
- across-cell 指标要求该基因在 test 集有 ≥30 个观测细胞；gene-mean baseline 无 across-cell 值（恒定预测器）。
- Ridge 是线性预测器，结论限定为"线性可读的细胞状态信息"。
