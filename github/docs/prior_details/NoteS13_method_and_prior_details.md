# 补充材料 S13：评分对象、统计口径与溯源细节

本文件仅保留与当前 GigaScience 稿件直接对应的内容。分析目标是比较已发布 processed reference products 对 benchmark reading 的影响；不再包含早期 LayerDiag、CPI、遗传侧或其他已退出当前稿件主线的材料。

## S13.1 固定输出参考产品替换对比

主分析定义 `Δ_ref = M(P_fixed, E_A) − M(P_fixed, E_B)`。`P_fixed` 为同一冻结预测，`E_A/E_B` 为两套已发布处理后参考产品，`M` 为固定评分协议。共同可评集合、留出单位、聚合方式和指标均保持不变。

## S13.2 matched training-product × evaluation-product 2×2

正文 2×2 使用与主结果完全一致的 11,275 个共同 `(cell line, drug)` 评测对。两个 prediction states 为 beta-trained 和 dcic2021-trained，分别在 beta 与 dcic2021 reference 上评分。difference-in-differences interaction 定义为 `I = (M_BB − M_BD) − (M_DB − M_DD)`；完整四格结果见 S25。

## S13.3 聚类 bootstrap

主分析两臂、配对差和 2×2 interaction 以 held-out drug 为聚类重采样单位。正文主比较采用 B = 2,000；2×2 同样基于 2,037 个 held-out drugs 进行 compound-cluster bootstrap。

## S13.4 共同可评集合

正文主分析为 64 个共同细胞系与 11,275 个共同 `(cell line, drug)` 对。更大的 `E_common = 15,990` 仅作为 2×2 敏感性分析，避免把不同共同集合上的绝对数值混为一谈。

## S13.5 训练暴露敏感性

strict-unseen 和 fingerprint-clean 分析只覆盖当前可追溯的适配训练暴露。更早 base checkpoint 的完整训练组成未保留，因此这些子集结果不能排除所有历史训练暴露，只用于判断当前可追溯的直接药物暴露或完全相同 ECFP4 指纹是否足以解释主差异。逐药结果见 S22，身份审计见 S19。

## S13.6 评分对象标识原则

本文将 source/accession、processing level、named release/build 与实际进入评分的 processed object 区分记录。正文数值对应的评分对象通过不可变文件、manifest 与 checksum 锁定；12 项 L1000 reporting/provenance audit 见 S21，代表性对象追踪见 S21a。

## S13.7 `dcic2021` processed reference 的完整溯源链

本说明把正文对照臂使用的 dcic2021 评价对象（`data/g2cp_cache_2021/`）从「上游分发」到「实际进入评分的矩阵」逐跳串成一条可复核的链，并给出逐文件校验值。目的：使「被评分的确切文件」不再依赖论文中的文字描述，而可由校验值唯一锁定。

## 溯源链（自顶向下）

1. **上游分发（SigCom / DCIC 2021）**
   - 原始分发文件：`cp_coeff_mat.gctx`（LINCS-sigs-2021 的 CD 系数矩阵；CD = Characteristic Direction，即与本文 β 主源 MODZ 不同的处理族）。
   - 来源 URL：`https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gctx/cd-coefficient/cp_coeff_mat.gctx`（取自本工程 `Downloads/G2CP_push/download_level5.py`）。
   - 本地副本：`Desktop/cp_coeff_mat.gctx`，36,084,518,760 B，MD5 `12ee07d08c0368e77914113b3681eb39`（计算日期 2026-09-27）。
   - 获取日期：未在案（已如实声明）。

2. **`step44a_build_pairs_2021.py`** — 从 gctx 构造 (药物 × 细胞系) 配对均值矩阵
   - 过滤：细胞系取 v7-162 词表（双向别名映射），药物取 v7 32,039 词表。
   - 输出：`rows_gctxorder.npy`（156,931 × 12,327）+ `pairs_meta.npz`（cell_name / cell_idx / drug_id / drug_idx）。

3. **`step44b_map_genes.py`** — v7 HVG（12,328）↔ gctx 列映射
   - 通过 HGNC 规范化（来源 `hgnc_symbols.tsv`、`geneinfo_beta.txt`）把 v7 的 12,328 个 HVG 映射到 gctx 的 12,327 个列。
   - 输出：`gene_map.json`（长度 12,328；gctx 列下标或 −1）。

4. **`step44c_assemble_2021.py`** — 装配最终缓存
   - 输出：`y.npy`（156,931 × 12,328，按 HVG 顺序）+ `meta.npz`（`kind/key/cell/gene_vocab/drug_vocab/cl_names/hvg/col2row`，沿用基座词表索引）+ `drug_fps.npy`（复制自 cache_full，对齐 ckpt 药物词表）。

5. **冻结的评价对象（实际进入评分的矩阵）**
   - `data/g2cp_cache_2021/y.npy` —— 即 dcic2021 臂的「评价对象」。
   - 正文对照臂在该对象上的留出 PCC = 0.0724（受控 11,275 对表）。

## 逐文件校验值（`results/dcic2021_cache_checksums_20260925.json`，计算日期 2026-09-25）

| 文件 | MD5 | 大小 (B) | 备注 |
|---|---|---|---|
| `y.npy` | `be8c09dc312f0ee7ba4c10a1cb0979dd` | 7,738,581,600 | SHA256 `351362fca7e0326b9594fbb92bf9bfae171534d63ff171073922d97767faa9c4`（评价对象本体） |
| `rows_gctxorder.npy` | `88cecb25e51c99e89034ad8f2641e858` | 7,737,953,876 | 步骤 2 中间产物 |
| `pairs_meta.npz` | `6b9fdafac216ab77a008b525d770badb` | 2,690,208 | 配对元数据 |
| `gene_map.json` | `d4fbe9bd879834bc2cdd9cabe605bd48` | 75,177 | 步骤 3 产物 |
| `meta.npz` | `15eb3c8e45a04ab3f0ddb9c0bbb189bc` | 2,169,832 | 词表 / 索引 |
| `drug_fps.npy` | `302f2eec4ad58dcc215c2b21a43dfc2f` | 262,463,616 | 药物指纹 |
| `sigidx_maxDose.npz` | `9828af88cc3907f7a8dc47f860fb03d3` | 3,206,308 | 签名索引 |
| `sigidx_sigALL.npz` | `bf90676af82bbf24c5a0aa4f2f4e8d4f` | 14,531,684 | 签名索引 |
| `sigidx_sigT24.npz` | `a8ef5ff5b755f769e7065ed2d97e5095` | 9,945,060 | 签名索引 |

## 与 β 主源评价对象的对称对照（用于正文「只换文件、其余不变」的断言）

| 跳 | β 主源（level5beta2020） | dcic2021 |
|---|---|---|
| 上游文件 | `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`（GEO GSE92742，Level 5 COMPZ MODZ） | `cp_coeff_mat.gctx`（LINCS-sigs-2021 CD 系数） |
| 构建脚本 | `step58b_beta_build_v2.py`（+ siginfo / cellinfo / geneinfo） | `step44a/b/c_*` |
| 评价对象 | `data/g2cp_cache_beta_v2/y.npy` | `data/g2cp_cache_2021/y.npy` |
| 留出 PCC（受控 11,275 对） | **0.3680** | **0.0724** |
| 唯一被换的层 | — | 处理族（MODZ → CD）+ 发布版本（2020 → 2021） |

> 两臂共用同一基座 `g2cp_full_cpi_v7.pt`、同一词表（cl_names v7-162、drug_vocab v7 32,039）、同一留出划分（后 10% 共 2,037 个药物）与同一 11,275 对评测表；见正文 3.2、4.1.1–4.1.4 与补充材料 S14.6。

