# 外部数据来源与重建路径

| 数据 | 版本/标识 | 获取方式 | 本包是否含 | 许可 |
|---|---|---|---|---|
| LINCS L1000 level-3/level-5（beta2020） | GSE92742；MODZ level-5 | GEO / CLUE.io | 否（41 GB gctx） | 公开可用，见 GEO 条款 |
| LINCS dcic2021（Characteristic Direction） | SigCom LINCS / DCIC 2021 | sigcom.lincs.mssm.edu / clue.io | 否 | 公开可用 |
| LINCS 评价网格 11,275 行 | `eval11275_trt_inst_ids.txt`（本包 github 未含，见下） | 由 inst_id 清单确定性重建 | 清单需补 | — |
| DepMap 21Q2 Public v2 依赖矩阵（CERES/Chronos） | 21Q2 Public v2 | depmap.org/portal/download | 派生子集在 zenodo | CC-BY 4.0 |
| DepMap 21Q2 表达矩阵（阳性对照输入） | `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | depmap.org | **否（需重建，见下）** | CC-BY 4.0 |
| DeepDEP 公开预测 | 原论文补充件 | 原论文 | 是（`source_data/.../first_predictor_DeepDEP`） | 按原论文 |
| GTEx v11 RNASeQC v2.4.3 / RSEM v1.3.3 | 2025-08-22 release | gtexportal.org | 派生数据在包内；原始 5.9 GB×2 未含 | GTEx 条款 |
| MSigDB C2:CGP（ACEVEDO_LIVER…/RODRIGUES_THYROID…/CASORELLI_APL…） | v2023.2 | msigdb.org（需注册） | **否（许可限制）** | MSigDB 条款，禁止整库再分发 |

## 两处需要作者补的东西（影响"逐位复现"）

1. **LINCS 11,275×978 参考矩阵**（约 180 MB，float32）：目前只在分析时从 41 GB gctx 现读，未落盘。
   补法二选一：(a) 一次性导出 `Y_beta2020_11275x978.npy` + `Y_dcic2021_11275x978.npy` 放 Zenodo，
   审稿人即可完全跳过 41 GB 下载；(b) 保留 `eval11275_trt_inst_ids.txt` + 重建脚本，审稿人自行下载公开产品重建。
   **建议 (a)+(b) 都做。**
2. **DepMap 阳性对照的 `X_TPM.npy`**（908×19,177，约 69 MB）：工作区内已无此文件。
   补法：从 DepMap 21Q2 `OmicsExpressionProteinCodingGenesTPMLogp1.csv` 按 `cell_order.txt` 抽取 908 个细胞系、
   保留 19,177 个蛋白编码基因，写成 float32 的 `X_TPM.npy` 一并放 Zenodo（并附重建脚本 + 校验和）。
   否则阳性对照只能"重跑到同结果"，无法校验输入矩阵本身。
