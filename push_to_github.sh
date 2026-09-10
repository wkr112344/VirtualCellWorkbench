#!/usr/bin/env bash
# push_to_github.sh - 一键推 VirtualCellWorkbench 到 GitHub
#
# 硬护栏（每条都会跑，绝对不会触碰 G2CP-virtual-cell）：
#   1. remote URL 必须包含 "VirtualCellWorkbench"
#   2. remote URL 不能包含 "G2CP"
#   3. push 前记 G2CP-virtual-cell HEAD SHA，push 后必须一字不差
#   4. 只能 push main + tag v0.1.0，不接受 --force
#
# 用法：
#   1. 创 PAT（fine-grained, Public Repos, Contents: R/W）
#   2. export GH_TOKEN=github_pat_xxxxx
#   3. GITHUB_USER=wkr112344 bash push_to_github.sh

set -euo pipefail
cd "$(dirname "$0")"

# ---------- 配置 ----------
GITHUB_USER="${GITHUB_USER:-REPLACE_WITH_YOUR_USERNAME}"
NEW_REPO="VirtualCellWorkbench"        # 新仓（要建、要推的）
PROTECTED_REPO="G2CP-virtual-cell"    # 预印本里引用的仓——绝对不能碰
TAG="v0.1.0"

if [ "$GITHUB_USER" = "REPLACE_WITH_YOUR_USERNAME" ]; then
    echo "ERROR: 还没填 GitHub 用户名"
    echo "  GITHUB_USER=wkr112344 bash push_to_github.sh"
    exit 1
fi
if [ -z "${GH_TOKEN:-}" ]; then
    echo "ERROR: 还没设 GH_TOKEN 环境变量"
    echo "  export GH_TOKEN=github_pat_xxxxx"
    exit 1
fi

# ---------- 硬护栏 1：URL 检查 ----------
NEW_URL="https://github.com/${GITHUB_USER}/${NEW_REPO}.git"
if [[ "$NEW_URL" == *"G2CP"* ]]; then
    echo "FATAL: URL 含 'G2CP'，护栏触发，退出"
    echo "  解析出来的 URL = $NEW_URL"
    exit 1
fi
if [[ "$NEW_URL" != *"VirtualCellWorkbench"* ]]; then
    echo "FATAL: URL 不含 'VirtualCellWorkbench'，护栏触发，退出"
    exit 1
fi
echo "[护栏 1 ✓] URL = $NEW_URL"

# ---------- 硬护栏 2：本地仓库状态 ----------
if [ ! -d .git ]; then
    echo "FATAL: 当前目录不是 git 仓库"
    exit 1
fi
CURRENT_BRANCH=$(git branch --show-current)
CURRENT_HEAD=$(git rev-parse HEAD)
echo "[护栏 2 ✓] branch=$CURRENT_BRANCH  HEAD=${CURRENT_HEAD:0:8}"

if [ "$CURRENT_BRANCH" != "main" ]; then
    echo "FATAL: 当前分支不是 main，是 '$CURRENT_BRANCH'"
    exit 1
fi

# ---------- 硬护栏 3：记 G2CP-virtual-cell push 前 SHA ----------
PROXY="${HTTPS_PROXY:-http://127.0.0.1:51081}"
BEFORE_SHA=$(curl -sS --max-time 15 -x "$PROXY" -H "User-Agent: vcw" \
    "https://api.github.com/repos/${GITHUB_USER}/${PROTECTED_REPO}/commits/master" \
    | python -c "import sys,json; print(json.load(sys.stdin)['sha'][:12])")
echo "[护栏 3 ✓] G2CP-virtual-cell master HEAD (before) = $BEFORE_SHA"

# ---------- 创建新仓（如不存在） ----------
echo
echo "==> 检查 ${NEW_REPO} 是否存在"
EXIST_CODE=$(curl -sS --max-time 15 -x "$PROXY" -H "User-Agent: vcw" \
    -o /dev/null -w "%{http_code}" \
    "https://api.github.com/repos/${GITHUB_USER}/${NEW_REPO}")
echo "    HTTP $EXIST_CODE"

if [ "$EXIST_CODE" = "404" ]; then
    echo "==> 创建新仓 ${NEW_REPO}（空仓，无 README/license/.gitignore）"
    curl -sS --max-time 30 -x "$PROXY" \
        -H "Authorization: token ${GH_TOKEN}" \
        -H "User-Agent: vcw" \
        -H "Accept: application/vnd.github+3" \
        -X POST "https://api.github.com/user/repos" \
        -d "{\"name\":\"${NEW_REPO}\",\"description\":\"Companion code to Wei 2026 (bioRxiv DOI 10.64898/2026.08.10.743942) — VirtualCellWorkbench\",\"private\":false,\"auto_init\":false}" \
        -o /dev/null -w "    create HTTP %{http_code}\n"
else
    echo "==> 仓已存在，跳过创建"
fi

# ---------- 推送（用临时 remote，避免污染全局 git config）----------
echo
echo "==> 推送 main 分支到 ${NEW_REPO}（token 嵌入 URL，不污染 git config）"
AUTH_URL="https://x-access-token:${GH_TOKEN}@github.com/${GITHUB_USER}/${NEW_REPO}.git"

# 硬护栏 4：再次 URL 检查
if [[ "$AUTH_URL" == *"G2CP"* ]] || [[ "$AUTH_URL" != *"VirtualCellWorkbench"* ]]; then
    echo "FATAL: push URL 护栏触发"
    echo "  AUTH_URL = ${AUTH_URL}"
    exit 1
fi

HTTPS_PROXY="$PROXY" git push "$AUTH_URL" "main:main"

echo
echo "==> 推送 tag ${TAG}"
HTTPS_PROXY="$PROXY" git push "$AUTH_URL" "$TAG"

# ---------- 硬护栏 5：G2CP-virtual-cell push 后 SHA 必须未变 ----------
echo
echo "==> 验证 G2CP-virtual-cell 未被触碰"
AFTER_SHA=$(curl -sS --max-time 15 -x "$PROXY" -H "User-Agent: vcw" \
    "https://api.github.com/repos/${GITHUB_USER}/${PROTECTED_REPO}/commits/master" \
    | python -c "import sys,json; print(json.load(sys.stdin)['sha'][:12])")
echo "    G2CP-virtual-cell master HEAD (after) = $AFTER_SHA"

if [ "$BEFORE_SHA" = "$AFTER_SHA" ]; then
    echo "[护栏 5 ✓] G2CP-virtual-cell 未变（SHA 一致）"
else
    echo "[护栏 5 ✗ FATAL] G2CP-virtual-cell HEAD 变了！before=$BEFORE_SHA after=$AFTER_SHA"
    exit 1
fi

# ---------- 验证新仓收到正确内容 ----------
echo
echo "==> 验证新仓 VirtualCellWorkbench HEAD"
NEW_HEAD=$(curl -sS --max-time 15 -x "$PROXY" -H "User-Agent: vcw" \
    "https://api.github.com/repos/${GITHUB_USER}/${NEW_REPO}/commits/main" \
    | python -c "import sys,json; print(json.load(sys.stdin)['sha'][:8])")
echo "    VirtualCellWorkbench main HEAD = $NEW_HEAD (期望 ${CURRENT_HEAD:0:8})"

if [ "$NEW_HEAD" = "${CURRENT_HEAD:0:8}" ]; then
    echo "[护栏 6 ✓] 新仓 HEAD 与本地一致"
else
    echo "[护栏 6 ✗ WARN] 新仓 HEAD 与本地不一致！本地=${CURRENT_HEAD:0:8} 远端=$NEW_HEAD"
    exit 1
fi

echo
echo "==================================="
echo "✅ 推送完成，5 道护栏全过。"
echo "   - G2CP-virtual-cell: 未触碰（SHA 一致）"
echo "   - VirtualCellWorkbench: 新建 + 推送 main + tag $TAG"
echo "   - 预印本 (DOI 10.64898/2026.08.10.743942): 未受影响"
echo
echo "下一步：去 https://zenodo.org 设置 GitHub 集成，详见 ZENODO.md"
echo "==================================="