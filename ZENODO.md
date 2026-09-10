# Zenodo 流程（5 步，约 2 分钟）

## 前置
- GitHub 账号 + 一个空的 VirtualCellWorkbench 仓库（public）
- Zenodo 账号（用 GitHub 账号登录即可，[zenodo.org](https://zenodo.org)）

## 步骤

### 1. 推 GitHub 仓库（前面已做）
```bash
cd workspace/VirtualCellWorkbench
# 替换 GITHUB_USER 为你的用户名
bash push_to_github.sh
```

### 2. 关联 Zenodo + GitHub
1. 登录 [zenodo.org](https://zenodo.org)（用 GitHub 登录）
2. 右上角头像 → Settings → GitHub
3. 点 "Sync now" 或 "Connect"
4. 在仓库列表里勾选 `VirtualCellWorkbench`
5. 保存

### 3. 创建 GitHub release
1. 打开 `https://github.com/<你的用户名>/VirtualCellWorkbench/releases`
2. 点 "Draft a new release"
3. Tag: `v0.1.0`（前面 push_to_github.sh 已经 push 这个 tag）
4. Title: `VirtualCellWorkbench v0.1.0`
5. Description:
```
VirtualCellWorkbench: A toolkit for chemical-perturbation response prediction at cell-line scale.

Companion code to:
Wei K. (2026) Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual Perturbation Screening at Scale. bioRxiv DOI: 10.64898/2026.08.10.743942

Headline numbers (paper §3.1 protocol):
  v7_raw on official 10% held-out perturbations: PCC = 0.341 [0.334, 0.348]
  Reported in bioRxiv preprint: PCC = 0.305
  Random 10% drug hold-out (E1 baseline): PCC = 0.137
  Leave-3-cell-out (E5 unseen cell lines): PCC = 0.222
  E8 Vemurafenib cross-cell PCC (A375 vs A549): -0.27

Eleven experiments (E0-E11) reproduce all main claims. See docs/EXPERIMENTS.md.
```
6. 点 "Publish release"

### 4. 等 Zenodo 自动打包
- 通常 1-3 分钟
- 打包完成后会在 release 页面下方出现 Zenodo badge

### 5. 拿 DOI
- 打开 release 页面
- 找到 "This DOI" 一行，复制（例如 `10.5281/zenodo.1234567`）
- 把 DOI 填回 `CITATION.cff` 和 `README.md` 的 bibtex 块
- 改完 commit + push

## 验证
DOI 拿到后，访问 `https://doi.org/<你的 DOI>` 应该能跳到 Zenodo 归档页面，看到 `VirtualCellWorkbench-v0.1.0.zip` 可以下载。

## 论文里写
```bibtex
@software{wei2026virtualcellworkbench,
  author = {Wei, Kairui},
  title  = {VirtualCellWorkbench: A toolkit for chemical-perturbation response prediction at cell-line scale},
  year   = {2026},
  version = {0.1.0},
  url    = {https://github.com/<你的用户名>/VirtualCellWorkbench},
  doi    = {10.5281/zenodo.<编号>}
}
```
