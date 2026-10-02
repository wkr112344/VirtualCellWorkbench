# Supplementary Table S21a. 四个代表性 score-to-object 追踪案例

本表从完整 12 项定向 audit（Table S21）中选择 4 个案例，用于说明 source/accession、processing level、named release、study-specific artifact 与 exact processed object 可以是不同的 recoverability layers。**第 1–2 例是与 Level-5 build ambiguity 最直接的案例；第 3–4 例仅用于展示 provenance layer 的区分，不意味着这些 Level-3 研究本应报告 level5beta2020/dcic2021。**

| # | 研究 | 公开报告到的层级 | 是否直接给出与本文 Level-5 对照同粒度的 build 标识 | 本文如何使用该案例 |
|---|---|---|---|---|
| 1 | DeepCE | 明示 L1000 Level 5；公开 processed train/dev/test 与代码 | 否 | **直接案例**：说明“Level 5”本身不是与 `level5beta2020` / `dcic2021` 同粒度的 build identifier。本文不据此断言原研究不可复现。 |
| 2 | MiTCP | 明示 GSE92742/GSE70138 与 Level 5；预处理数据已归档 | 否 | **直接案例**：即使 broad source 和 Level 5 都清楚，paper-level descriptor 与 exact processed build 仍可处在不同层。 |
| 3 | TranSiGen / XPert | 明示 CMap LINCS Resource 2020、Level 3，并有 processed artifacts | 不适用（其分析为 Level 3） | **层级案例**：说明 named release、processing level 与 study-specific processed artifact 可以分别承担不同的可恢复性功能；不用于指控其缺少 Level-5 build。 |
| 4 | Dr.VAE | 明示 CMap-L1000v1、Level 3 | 不适用（其分析为 Level 3） | **层级案例**：说明 named dataset/version 并不等于本文所说的 exact metric-generating processed object；是否需要更细标识取决于具体 benchmark。 |

**解释边界。** 该 audit 是 targeted context，不估计领域发生率，也不评价被审计论文的总体 reproducibility 质量。本文只用它说明：当一个 benchmark 生态确实存在多个可作为 reference 的 processed products 时，score 的可审计性最好直接绑定到实际评分对象、共同可评集合与 protocol。
