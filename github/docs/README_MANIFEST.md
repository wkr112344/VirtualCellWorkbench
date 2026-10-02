# DepMap 原始产物 · 供「共同基因背景对照（残差化 2×2）」实验

生成：2026-09-29（仅复制与校验，未改原始数据）

## 网格与行/列口径（**务必按此读取**）

- 全部 908 × 17393 矩阵的行序 = `index/cell_order.txt`
- 列序 = `index/target_genes.txt`，格式为 `SYMBOL (EntrezID)`，如 `A1BG (1)`
  → 与外部基因名匹配时须去掉 ` (id)` 后缀再比对
- `*_test.npy` 的行序 = `index/test_cell_order.txt`（136 个）
  → 已验证：与 `split.json` 的 `test_cells` **顺序完全一致**
- 切分来自 `index/split.json`（**seed 42**，训练时锁定、不得重切）：
  train 636 / val 136 / test 136

## 目录

| 路径 | 内容 | 形状 |
|---|---|---|
| `reference/Y_CERES.npy` | CERES 参考 | (908, 17393) float32 |
| `reference/Y_Chronos.npy` | Chronos 参考 | (908, 17393) float32 |
| `reference/common_observed_mask.npy` | 共同观测掩码（NaN 处理，勿插补） | (908, 17393) bool |
| `predictions_test/seed{0-4}_pred_ceres_arm_test.npy` | CERES-trained 臂 test 预测 | (136, 17393) float32 |
| `predictions_test/seed{0-4}_pred_chronos_arm_test.npy` | Chronos-trained 臂 test 预测 | (136, 17393) float32 |
| `first_predictor_DeepDEP/DeepDEP_predictor_final.csv` | 第一 predictor 冻结预测 | 273 × 1244 |
| `first_predictor_DeepDEP/DeepDEP_predictor_strict.csv` | 同上 strict 版 | 272 × 1244 |

`index/depmap_index.json` 提供机器可读索引，含
`train_rows_in_cell_order` / `test_rows_in_cell_order` —— 取"仅训练细胞的 gene mean"
直接用 `train_rows_in_cell_order` 索引 Y 即可。

## ⚠ 已知约束（影响"公平比较"能否成立）

第一 predictor（DeepDEP）与第二 predictor（VAE-DeepDEP）**不在同一网格**：

| 项 | final | strict |
|---|---|---|
| shape | 273 × 1244 | 272 × 1244 |
| 基因 ∩ GH 17393 | 1244 | 1244 |
| 细胞 ∩ 全部 908 | 272 | 271 |
| **细胞 ∩ test 136** | **43** | **43** |
| 细胞 ∩ train 636 | 190 | 189 |

**结论**：基因 1244 个全部落在 GH 的 17393 里，但细胞只有 43/136
个落在 test 集。若要把 DeepDEP 与两个训练臂一起做"公平比较"，评测会被迫缩到
**43 个 test 细胞 × 1244 基因**。

这是设计层面的取舍，两种处理：
1. **主分析只用第二 predictor**（136 test 细胞 × 17393 基因，网格完整）；
2. DeepDEP 作为补充对照，在 43×1244 子网格上做，并明确标注该比较受网格限制。

## 原始数据位置（未改动）

- 数据：`C:\Users\wkr20\Desktop\depmap_second_corpus\work\DeepDEP_training_data\aligned_21Q2_908_3omics`
- 预测：`C:\Users\wkr20\Desktop\depmap_second_corpus\GH_second_predictor\predictions`
- 第一 predictor：`C:\Users\wkr20\Desktop\depmap_second_corpus\work\fixed_eval_input`
