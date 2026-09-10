#!/usr/bin/env bash
# push_to_github.sh - 一键推 VirtualCellWorkbench 到 GitHub
#
# 用法（首次推送）：
#   1. 在 https://github.com/new 创建一个空仓库（public）
#      Repository name: VirtualCellWorkbench
#      不要勾选 Add a README / Add .gitignore / Choose a license
#   2. 编辑本脚本，把下面 GITHUB_USER 改成你的用户名
#      （或者直接 GITHUB_USER=你的用户名 bash push_to_github.sh）
#   3. bash push_to_github.sh
#      首次 push 会弹窗让你输 GitHub 用户名/密码（或 Personal Access Token）
#
# 二次推送（已有 remote）：
#   bash push_to_github.sh

set -e
cd "$(dirname "$0")"

# ---------- 配置 ----------
GITHUB_USER="${GITHUB_USER:-REPLACE_WITH_YOUR_USERNAME}"
REPO="VirtualCellWorkbench"
TAG="v0.1.0"

if [ "$GITHUB_USER" = "REPLACE_WITH_YOUR_USERNAME" ]; then
    echo "ERROR: 还没填 GitHub 用户名"
    echo "  方法 A：编辑本脚本，把 GITHUB_USER 改成你的用户名"
    echo "  方法 B：临时传参：  GITHUB_USER=你的用户名 bash push_to_github.sh"
    exit 1
fi

REMOTE_URL="https://github.com/${GITHUB_USER}/${REPO}.git"

echo "==> 目标：${REMOTE_URL}"
echo "==> 当前 HEAD：$(git rev-parse --short HEAD)"
echo "==> 当前 tag："
git tag --list | sed 's/^/    /' || true
echo

# ---------- 加 remote ----------
if git remote get-url origin >/dev/null 2>&1; then
    echo "==> remote origin 已存在，更新 URL"
    git remote set-url origin "${REMOTE_URL}"
else
    echo "==> 添加 remote origin"
    git remote add origin "${REMOTE_URL}"
fi

# ---------- 推送 main ----------
echo
echo "==> 推送 main 分支"
git push -u origin main

# ---------- 推送 tag ----------
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo
    echo "==> 推送 tag ${TAG}"
    git push origin "$TAG"
else
    echo
    echo "==> 本地没有 ${TAG} tag，跳过 tag 推送"
    echo "    如需打 tag：git tag -a ${TAG} -m '...' && git push origin ${TAG}"
fi

echo
echo "==> 完成。Next：去 https://zenodo.org 设置 GitHub 集成，详见 ZENODO.md"