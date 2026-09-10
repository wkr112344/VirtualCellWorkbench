#! /bin/bash
# push_to_github.sh
# 一键推 VirtualCellWorkbench 到 GitHub
# 用法：
#   1. 先在 https://github.com/new 创建一个空仓库：VirtualCellWorkbench（public, no README/no .gitignore/no license）
#   2. 把 <你的 GitHub 用户名> 替换成实际用户名
#   3. bash push_to_github.sh

set -e
GITHUB_USER="REPLACE_WITH_YOUR_USERNAME"
REPO="VirtualCellWorkbench"

cd "$(dirname "$0")"

echo "→ 添加 remote"
git remote add origin "https://github.com/${GITHUB_USER}/${REPO}.git" 2>/dev/null || \
    git remote set-url origin "https://github.com/${GITHUB_USER}/${REPO}.git"

echo "→ 推送到 main 分支"
git push -u origin main

echo "→ 创建 v0.1.0 tag"
git tag -a v0.1.0 -m "VirtualCellWorkbench v0.1.0: initial release

Companion code to Wei 2026 (bioRxiv DOI 10.64898/2026.08.10.743942).

E0-E11 reproduction suite: PCC=0.341 (v7_raw) on official 10% drug hold-out,
PCC=0.222 on leave-3-cell-out (E5), Vemurafenib cross-cell PCC=-0.27 (E8)."

git push origin v0.1.0

echo ""
echo "✓ 推送完成。下一步：去 https://zenodo.org 设置 GitHub 集成，详见 ZENODO.md"
