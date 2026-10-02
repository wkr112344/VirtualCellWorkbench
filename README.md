# 处理后参考产品选择对生物医学模型评测与候选筛选的影响 —— 复现材料

本目录分两部分：
- `github/` —— 推到 GitHub 的**代码与文档**（版本控制、便于提 issue/PR）
- `zenodo/` —— 上传 Zenodo 的**派生数据、图件、稿件与校验**（不可变归档，给 DOI）

## 复现最短路径（审稿人视角）

1. 从 Zenodo 记录下载 `zenodo/source_data/`（约 250 MB）：冻结预测矩阵、DepMap 派生数据、
   GTEx 逐样本得分、图表源数据。
2. 从 GitHub 克隆代码仓库，`pip install -r requirements.txt`。
3. 跑 `run_all.sh`：脚本按顺序重算 表 1–4、图 1–6 与补充表 S1–S22，并与 `EXPECTED_OUTPUTS.json`
   中的关键数字逐条比对（容差 5×10⁻⁴）。
4. 需要从原始公共数据重跑的部分（LINCS level-5、DepMap 21Q2、GTEx v11）见 `DATA_SOURCES.md`；
   其中 LINCS 侧提供了 inst_id 清单用于确定性重建。

## 不能分发的部分（已在 DATA_SOURCES.md 写明）

- LINCS L1000 原始/level-5 产品（GEO GSE92742、SigCom LINCS、DCIC 2021）：公开但体量大，给下载与重建路径
- DepMap 21Q2 表达/依赖矩阵：公开（CC-BY 4.0），给下载与子集构建路径
- GTEx v11 原始表达：公开（dbGaP/GTEx portal），给下载路径
- MSigDB C2:CGP 三个疾病签名：受 MSigDB 许可约束，本包只给集合名与版本，基因清单请自行按许可下载
