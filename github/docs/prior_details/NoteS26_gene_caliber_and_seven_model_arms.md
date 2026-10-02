# Note S26 — Caliber disambiguation: 978 vs 963 gene axis, and the seven-model re-evaluation arms

> 用途：本说明为**最终英文稿**（正文 Methods / Data availability）与 **Zenodo version-of-record README** 提供可直接粘贴的澄清文字，避免审稿人把 978 与 963、以及"七模型复评两臂"误读为口径矛盾。
> 数据底表：`TableS24_seven_predictor_fixed_output.csv`（逐模型两臂 PCC、ΔPCC 与 bootstrap 95% CI）。
> 数值真值来源：`220_dpb_crossversion.json`（schema `dpb_crossversion/v1`，gene axis = 978，n_rows_paired = 9562，n_drugs = 1422，drug-level PCC，bootstrap B = 2000，seed 20260920）。

---

## A. 基因口径：978 与 963 是同一基因轴的两种对齐层次，不是维度矛盾

**结论（请写入英文稿 + Zenodo README）：**

> **Gene-axis caliber (978 vs 963).** All primary scoring in this study is performed on the **978 LINCS L1000 landmark genes**, identified by Entrez / canonical HGNC symbols and sliced from the 12,328-gene space via the fixed 12,328→978 column index map (`internal_col`, length 978; canonicalization source `hgnc_symbols.tsv`, md5 `faf032b2…`). The figure **963** that appeared in an earlier supplement referred to a *naive* raw-release symbol alignment: 15 of the 978 landmark genes have since undergone HGNC symbol renames (old→new), so a release carrying current symbols fails a literal string match on exactly those 15, yielding 963 = 978 − 15. These 15 genes are **not** genuinely absent — they are recovered by HGNC alias / previous-symbol canonicalization, which restores the full 978. Thus **978 and 963 describe two layers of the same gene-axis alignment (canonical vs naive symbol match), not two different gene sets and not a dimensionality contradiction.** The main text reports 978 (the canonical count used for all primary scoring); the historical 963 reflected a raw-release naive-match intermediate and should not be read as a reduced gene panel.

要点（中文备忘）：
- 978 = LINCS landmarks，**带 HGNC 规范化**（978/978 全命中）。
- 963 = 对原始 release 做**不做规范化的符号直配**时，15 个发生 HGNC 改名的 landmark 没配上 → 978−15=963；规范化后回到 978。
- 二者是**同一基因轴的两种对齐层次**，不是两套基因、更不是维度前后矛盾。
- 注意：算法篇里另一个 "963 列" 是「978 landmarks ∩ 遗传语料实测可测列」的交集，成因不同，勿与本稿 963 混为一谈。

---

## B. 七模型复评的"两臂"定义（与图 3 图注的口径对齐）

**关键事实（避免审稿人误读）：**

- 七模型复评是 **fixed-output replacement contrast**：冻结各已发表预测器（CIGER / DeepCE / MultiDCP / PertDiT / PRnet / TranSiGen / XPert）的既有 prediction matrix，只替换 Level-5 **实测参考产品**重新评分。
- 每个冻结预测在三个参考上分别打分：**(i) 模型自身的 SDST 参考**（第三方加工版本，与 level5beta2020 逐行内容相关 0.3524，是各模型的"训练/原生"参考）、**(ii) LINCS level5beta2020 (MODZ)**、**(iii) LINCS dcic2021 (CD)**。
- 正文 4.1.5 节与图 3 报告的 **"配对 ΔPCC 中位数 0.3993、范围 0.2599–0.6450、各 95% CI 下界 > 0"** 对应的是 **ΔPCC = PCC(自身 SDST 参考) − PCC(dcic2021)**（即 `TableS24` 的 `DeltaPCC_own_minus_dcic2021` 一列；七模型该列 = 0.2889 / 0.4397 / 0.2599 / 0.3993 / 0.3122 / 0.6273 / 0.6450，中位 0.3993）。
- **图 3 图注里"臂 A = MODZ-based level5beta2020"对七模型面板是借用了图 1 的措辞、并不精确**：七模型复评的"臂 A"是各预测器**自身的 SDST 参考输出**（第三方加工，非 LINCS level5beta2020 实测产品）；level5beta2020 只作为 G2CP 主分析（图 1）的臂 A。
- 为完整透明，`TableS24` 同时给出 **PCC on LINCS level5beta2020** 与 **β2020-vs-dcic2021 的 ΔPCC**（= PCC_beta2020 − PCC_dcic2021，量级很小，约 0.01–0.10），因为冻结的外部预测在任一 LINCS 产品上本就得分很低；这进一步说明 0.3993 的幅度来自"自身参考 vs dcic2021"这一对照，而非"两版 LINCS 产品互换"。

**建议写入英文稿（4.1.5 节 + 图 3 图注修正）：**

> **Seven-model re-evaluation arms.** Each published predictor's prediction matrix was held fixed and re-scored against three references: the predictor's own SDST reference (a third-party processed version, nearest to level5beta2020 with row-content correlation 0.3524), the LINCS level5beta2020 (MODZ) product, and the LINCS dcic2021 (CD) product — all on the 978-gene landmark axis. The paired ΔPCC reported here is ΔPCC = PCC(own SDST reference) − PCC(dcic2021); across the seven predictors its median is 0.3993 (range 0.2599–0.6450) and every bootstrap 95% CI lower bound is above 0. (The two LINCS products alone, PCC(level5beta2020) − PCC(dcic2021), give a much smaller ΔPCC because the frozen external predictions score poorly on either LINCS product; see Table S24.) The earlier "arm A = level5beta2020" wording in the figure note applied to the G2CP main analysis (Fig. 1) and is not the arm definition for this seven-model panel.

**建议图 3 图注改为：**

> Figure 3. Cross-product re-evaluation of seven published predictors' fixed outputs. Each prediction matrix is held fixed; only the Level-5 reference product used for scoring is changed — the predictor's own SDST reference versus the LINCS dcic2021 (CD) product, on the 978-gene landmark axis. All seven predictors show a consistent direction of ΔPCC (median 0.3993; range 0.2599–0.6450; every bootstrap 95% CI lower bound > 0), but the magnitude differs across architectures.

---

## C. 与 Zenodo README 的衔接

请在 Zenodo version-of-record (`10.5281/zenodo.22720615`) 的 `README.md` §/Data-availability 段补入 A 节的基因口径说明（978 带 HGNC 规范化、963 为 naive 直配中间量），并明确：
- `results/` 内 `220_dpb_crossversion.json` 与 `dpb_crossversion_boot.npz`（sha256 `3140df42…`）是七模型复评逐模型 ΔPCC 与 bootstrap 向量的真值来源；
- 本补充包的 `TableS24_seven_predictor_fixed_output.csv` 是其逐模型可读底表。
