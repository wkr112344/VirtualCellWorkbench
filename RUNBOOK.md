# 推送 VirtualCellWorkbench 到 GitHub — 手动 runbook

## 当前状态（agent 已经做完的部分）

✅ LieBao 代理（127.0.0.1:51081）已通过 git 配置
✅ GCM 2.9.0 已配置为 github.com 的 credential helper
✅ Remote origin 已加：`https://github.com/wkr112344/VirtualCellWorkbench.git`
✅ 本地仓库有 3 个 commit 在 main 分支（HEAD = `02e816b`）
✅ 5 道硬护栏写在 push_to_github.sh 里

⚠️ **剩余一步必须你亲手做**：跑 `git push`，因为 GCM 需要弹浏览器让你点 OAuth 授权。

---

## 你做 3 步（在 git bash 或 Windows Terminal 里）

### 1. 进仓库目录

```bash
cd "C:/Users/wkr20/WorkBuddy/2026-09-06-00-53-37/workspace/VirtualCellWorkbench"
```

### 2. 推 main 分支（弹 OAuth 窗口）

```bash
git push -u origin main
```

**会发生什么**：
- GCM 检测到没缓存的 github.com 凭证
- 弹 Edge/Chrome 窗口，URL = `https://github.com/login/device` 或 OAuth 同意页
- 你点 **Authorize git-credential-manager**
- 回到终端，git push 自动继续
- 看到 `Writing objects: 100%` 即成功

### 3. 推 tag v0.1.0

```bash
git push origin v0.1.0
```

(tag 已经在本地，commit `02e816b` 之前打的，详见 `git tag -l`)

---

## 万一浏览器没弹

### 方案 A：手动拿 device code

```bash
# 在终端跑，会打印一段 8 位 device code + URL
git credential-manager github login
# 然后手动开 https://github.com/login/device 输入 code
```

### 方案 B：fallback 到 PAT

```bash
# 1. https://github.com/settings/tokens?type=beta → Generate
#    Name: VCW-publish, 7 days, Public Repos (R/W), Contents: R/W
# 2. 拿到 token 后：
export GH_TOKEN=github_pat_xxxxx
git push https://x-access-token:${GH_TOKEN}@github.com/wkr112344/VirtualCellWorkbench.git main
git push https://x-access-token:${GH_TOKEN}@github.com/wkr112344/VirtualCellWorkbench.git v0.1.0
```

---

## 推送后验证（你可以看，agent 也会自动验）

```bash
# 新仓是否建好
https://github.com/wkr112344/VirtualCellWorkbench

# G2CP-virtual-cell 必须**纹丝不动**
# 推送前 HEAD: 3e20583862ef  ("Renumber citations by first-appearance order...")
# 推送后 HEAD: 3e20583862ef  （应该一字不差）
```

如果 G2CP HEAD 变了，立刻来 ping agent，那是 push 脚本出 bug 了。

---

## 推送完之后

### 上 GitHub 创建空仓（**如果你打算用 GCM OAuth 走完上面流程就跳过这步**）

如果 GCM 流程失败要 fallback 到"自己点创建"：
1. 浏览器开 https://github.com/new
2. Repository name: **VirtualCellWorkbench**
3. Description: `Companion code to Wei 2026 (bioRxiv DOI 10.64898/2026.08.10.743942)`
4. **Public** / **不要**勾 Add a README / **不要**勾 Add .gitignore / **不要**勾 Choose a license
5. Create repository

### Zenodo（推送完之后独立做）

按 `ZENODO.md` 5 步：
1. https://zenodo.org → Login with GitHub
2. Settings → GitHub → 勾 VirtualCellWorkbench → Save
3. https://github.com/wkr112344/VirtualCellWorkbench/releases → Draft new release
4. Choose tag: v0.1.0（选已 push 的那个）→ Publish
5. 等 1-2 分钟 → Zenodo 自动存档 → 给你 DOI badge URL
6. 把 DOI 填回 CITATION.cff 和 README.md 的 badge

### 最后通知 agent

```
我推完了，DOI 是 10.5281/zenodo.xxxxx
```

agent 会把 DOI 填进：
- `CITATION.cff`（`identifiers:` 段）
- `README.md`（Zenodo badge）
- 然后 commit 一个新版本

---

## 关键安全提醒

⚠️ **预印本**（DOI 10.64898/2026.08.10.743942）一旦发了就不可修改，与 GitHub 仓无自动同步。

⚠️ **G2CP-virtual-cell**（预印本里引用的仓）最后 push 是 2026-08-09，HEAD = `3e20583862ef`，**任何推送后这个 SHA 必须不变**。

⚠️ **token 不进对话**：你只用 GCM OAuth（不传 PAT）或用 `export GH_TOKEN=...`（terminal 不入库）。