# -*- coding: utf-8 -*-
"""
R_big_independent_cohort.py —— 以「独立 cohort」为核心的大实验（双生态统一）

把 **cohort 独立性**本身作为核心设计变量，在两个数据生态上统一执行，跨正交切分轴检验：

  LINCS（参考产品对比：level5beta2020 vs dcic2021，冻结 β-trained 预测）
    轴 1：药物-based 互不相交 cohort（K=2,4,5）—— 已有
    轴 2：细胞系-based 互不相交 cohort（K=2,3,4）—— 本脚本新增（正交轴）
    三口径并行：Pearson / Spearman / cosine
    聚合单元仍为「药物」（与主文口径一致）：先把 pane 内 Δ 按药物取均值，再按药物 bootstrap

  DepMap（训练靶标对比：CERES vs Chronos，VAE-DeepDEP 5 seed）
    轴：细胞系-based 互不相交 cohort（K=2,4）—— 已有
    per-seed cellwise median → mean over seeds；hierarchical bootstrap (seed→cell)

产出（统一交付）：
  results/big_experiment_independent_cohorts.csv   逐 cohort（双生态）底表
  results/big_experiment_independent_cohorts.json  汇总 + 复制度量
  figures/big_experiment_independent_cohorts.png   双生态 forest 图
"""
import csv
import io
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = r"C:\Users\wkr20\WorkBuddy\Claw\gigascience_v5\results"
FIG = r"C:\Users\wkr20\WorkBuddy\Claw\gigascience_v5\results\for_figures\figures"
N_BOOT = 2000
BOOT_SEED = 12345
_t0 = time.time()


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def ci_of(per_drug, rng, n=N_BOOT):
    b = np.array([np.nanmean(per_drug[rng.integers(0, len(per_drug), len(per_drug))])
                  for _ in range(n)])
    return float(np.nanmean(per_drug)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def relabel(lines):
    """原 CSV 首行列名 -> 三口径 Δ 列名"""
    return None


# ══════════════════════════  LINCS 部分 ══════════════════════════
log("载入 LINCS 逐行三口径底表 ...")
per = []
with open(os.path.join(RES, "frozen_perrow_three_metrics.csv"), encoding="utf-8") as f:
    rd = csv.DictReader(f)
    for r in rd:
        per.append(r)
log("  rows=%d" % len(per))
cells = np.array([r["cell"] for r in per], dtype=object)
drugs = np.array([r["drug"] for r in per], dtype=object)
METRICS = ["pearson", "spearman", "cosine"]
DELTA = {}
for m in METRICS:
    a = np.array([float(r["%s_A_beta2020" % m]) for r in per], dtype=np.float64)
    b = np.array([float(r["%s_B_dcic2021" % m]) for r in per], dtype=np.float64)
    DELTA[m] = a - b

METRIC_KEY = {"pearson": "delta_pearson", "spearman": "delta_spearman", "cosine": "delta_cosine"}


def lincs_cohort(mask, metric, tag_extra=""):
    """cohort 内：先按药物聚合 Δ，再按药物 bootstrap（与主文聚合单元一致）"""
    dg = np.unique(drugs[mask])
    dpos = {d: j for j, d in enumerate(dg)}
    arr = np.empty(len(dg), dtype=np.float64)
    for d in dg:
        arr[dpos[d]] = np.nanmean(DELTA[metric][mask & (drugs == d)])
    mm, lo, hi = ci_of(arr, np.random.default_rng(BOOT_SEED))
    return mm, lo, hi, len(dg), int(mask.sum())


rows_out = []
lincs_summary = {}

# ---- 轴 2：细胞系-based 互不相交 cohort（正交轴，本脚本新增）----
uniq_cells = np.array(sorted(set(cells.tolist())), dtype=object)
log("LINCS 细胞系 %d 个，做 cell-based 独立 cohort ..." % len(uniq_cells))
lincs_summary["axis_cell"] = {}
for K in [2, 3, 4]:
    rngk = np.random.RandomState(31337 + K)
    perm = rngk.permutation(len(uniq_cells))
    belong = np.empty(len(uniq_cells), dtype=np.int64)
    for bi, ci_ in enumerate(perm):
        belong[ci_] = bi % K
    c2g = {c: belong[i] for i, c in enumerate(uniq_cells)}
    row_g = np.array([c2g[c] for c in cells.tolist()], dtype=np.int64)
    rec = {"K": K, "per_cohort": []}
    for c in range(K):
        mask = row_g == c
        ent = {"cohort": int(c)}
        for m in METRICS:
            mm, lo, hi, nd, npair = lincs_cohort(mask, m)
            ent[m] = {"delta": round(mm, 6), "ci95": [round(lo, 6), round(hi, 6)],
                      "excludes_zero": bool(lo * hi > 0)}
            rows_out.append(dict(ecosystem="LINCS", split_axis="cell_line", K=K, cohort=c,
                                 metric=m, n_pairs=int(npair), n_drugs=int(nd),
                                 estimate=round(mm, 6), ci95_low=round(lo, 6), ci95_high=round(hi, 6),
                                 excludes_zero=bool(lo * hi > 0),
                                 estimand="reference contrast Δ=PCC(β2020)−PCC(dcic2021)"))
        rec["per_cohort"].append(ent)
    est = np.array([e["pearson"]["delta"] for e in rec["per_cohort"]], dtype=np.float64)
    rec["consistency"] = {
        "all_same_sign": bool(np.all(est > 0) or np.all(est < 0)),
        "all_ci_exclude_zero": bool(all(e["pearson"]["excludes_zero"] for e in rec["per_cohort"])),
        "estimates_pearson": [round(float(x), 6) for x in est],
        "median_pearson": round(float(np.median(est)), 6),
        "range_pearson": [round(float(est.min()), 6), round(float(est.max()), 6)]}
    lincs_summary["axis_cell"]["K%d" % K] = rec
    log("  K=%d Δ(pearson)=%s 全部同号=%s 全部CI不跨0=%s"
        % (K, np.round(est, 4).tolist(), rec["consistency"]["all_same_sign"],
           rec["consistency"]["all_ci_exclude_zero"]))

# ---- 轴 1：药物-based（引用已算结果，纳入统一底表）----
log("并入 LINCS 药物-based cohort 结果（已有）...")
with open(os.path.join(RES, "frozen_cohort_replication.csv"), encoding="utf-8") as f:
    for r in csv.DictReader(f):
        K = int(r["K"])
        for m in METRICS:
            if r["n_drugs"] == "0":
                continue
            rows_out.append(dict(ecosystem="LINCS", split_axis="drug", K=K, cohort=int(r["cohort"]),
                                 metric=m, n_pairs=int(r["n_pairs"]), n_drugs=int(r["n_drugs"]),
                                 estimate=float(r[METRIC_KEY[m]]), ci95_low=float(r["lo_pearson"] if m == "pearson" else "0"),
                                 ci95_high=float(r["hi_pearson"] if m == "pearson" else "0"),
                                 excludes_zero=(r["pearson_excludes_zero"] == "True") if m == "pearson" else None,
                                 estimand="reference contrast Δ=PCC(β2020)−PCC(dcic2021)"))
lincs_summary["note_drug_axis"] = ("药物轴 CI 逐行写入已有 summary json；"
                                   "本统一底表对 drug 轴仅携带 pearson 的 CI（spearman/cosine CI 见 "
                                   "frozen_multimetric_summary.json cohort_replication）")

# ══════════════════════════  DepMap 部分 ══════════════════════════
log("并入 DepMap 细胞系-based cohort 结果 ...")
with open(os.path.join(RES, "depmap_cohort_replication.csv"), encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows_out.append(dict(ecosystem="DepMap", split_axis="cell_line", K=int(r["K"]),
                             cohort=int(r["cohort"]), metric="pearson",
                             n_pairs=int(r["n_cells"]), n_drugs="",
                             estimate=float(r["interaction"]),
                             ci95_low=float(r["ci95_low"]), ci95_high=float(r["ci95_high"]),
                             excludes_zero=(r["excludes_zero"] == "True"),
                             estimand="2×2 interaction I=(CC−CH)−(HC−HH)"))

with open(os.path.join(RES, "frozen_multimetric_summary.json"), encoding="utf-8") as f:
    lj = json.load(f)
with open(os.path.join(RES, "depmap_cohort_summary.json"), encoding="utf-8") as f:
    dj = json.load(f)

os.makedirs(FIG, exist_ok=True)
with open(os.path.join(RES, "big_experiment_independent_cohorts.csv"), "w", newline="", encoding="utf-8") as f:
    keys = ["ecosystem", "split_axis", "K", "cohort", "metric", "n_pairs", "n_drugs",
            "estimate", "ci95_low", "ci95_high", "excludes_zero", "estimand"]
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows_out)

# ══════════════════════════  复制度量汇总 ══════════════════════════
def counts(sub):
    if not sub:
        return dict(n=0, same_sign=None, all_excl=None)
    e = np.array([r["estimate"] for r in sub], dtype=np.float64)
    excl = [r["excludes_zero"] for r in sub]
    return dict(n=len(sub), same_sign=bool(np.all(e > 0) or np.all(e < 0)),
                all_excl=bool(all(x is True for x in excl)),
                n_excl=int(sum(1 for x in excl if x is True)),
                median=round(float(np.median(e)), 6),
                range=[round(float(e.min()), 6), round(float(e.max()), 6)])


one_per_split = {}
for r in rows_out:
    key = (r["ecosystem"], r["split_axis"], r["K"], r["cohort"])
    if r["metric"] == "pearson" and (key not in one_per_split):
        one_per_split[key] = r
pear = list(one_per_split.values())

big = {
    "schema": "big_experiment_independent_cohort/v1",
    "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    "question": "把 cohort 独立性作为核心设计变量：同一（成对）参考/靶标对比，能否在每个互不相交的独立 cohort 内复现？",
    "reference": {
        "LINCS_full_sample_pearson": lj["metrics"]["beta_trained"]["pearson"]["delta_A_minus_B"],
        "DepMap_full_sample_interaction": dj["full_sample"]["interaction"]},
    "repro_gates": {
        "LINCS_reproduced_maintext": bool(lj["repro_check"]["pass"]),
        "DepMap_reproduced_published": bool(dj["full_sample"]["repro_gate_pass"])},
    "by_view": {
        "LINCS_cell_axis": counts([r for r in pear if r["ecosystem"] == "LINCS" and r["split_axis"] == "cell_line"]),
        "DepMap_cell_axis": counts([r for r in pear if r["ecosystem"] == "DepMap"]),
        "LINCS_drug_axis": counts([r for r in pear if r["ecosystem"] == "LINCS" and r["split_axis"] == "drug"]),
        "ALL": counts(pear)},
    "lincs_cell_axis_detail": lincs_summary["axis_cell"],
    "elapsed_s": round(time.time() - _t0, 1)}
json.dump(big, io.open(os.path.join(RES, "big_experiment_independent_cohorts.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

log("=" * 96)
log("复制度量（pearson 口径，每个 split 取一个代表）")
for k, v in big["by_view"].items():
    log("  %-18s n=%-3d 同号=%-5s 全部CI不跨0=%-5s 中位=%.6f 范围=%s"
        % (k, v["n"], v["same_sign"], v["all_excl"], v.get("median", float("nan")), v.get("range")))
log("  参考：LINCS 全样本 Δ=%.6f  |  DepMap 全样本 I=%.6f"
    % (big["reference"]["LINCS_full_sample_pearson"][0], big["reference"]["DepMap_full_sample_interaction"]))
log("=" * 96)

# ══════════════════════════  图：双生态 forest ══════════════════════════
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3, "axes.axisbelow": True})
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
for ax, eco, title in [(axes[0], "LINCS",
                        "LINCS: reference contrast Δ (β2020 − dcic2021)\nfrozen β-trained prediction"),
                       (axes[1], "DepMap",
                        "DepMap: 2×2 interaction (CERES vs Chronos)\nVAE-DeepDEP, 5 seeds")]:
    sub = [r for r in pear if r["ecosystem"] == eco]
    ys = np.arange(len(sub))
    lab = []
    for r in sub:
        lab.append("%s K=%d c%d (n=%s)" % ({"cell_line": "cell", "drug": "drug"}[r["split_axis"]],
                                           r["K"], r["cohort"], r["n_pairs"]))
    ax.errorbar([r["estimate"] for r in sub], ys,
                xerr=[[max(0.0, r["estimate"] - r["ci95_low"]) for r in sub],
                      [max(0.0, r["ci95_high"] - r["estimate"]) for r in sub]],
                fmt="o", color="#2c7fb8", capsize=3, ecolor="gray")
    ref = (big["reference"]["LINCS_full_sample_pearson"][0] if eco == "LINCS"
           else big["reference"]["DepMap_full_sample_interaction"])
    ax.axvline(ref, color="#1b9e77", lw=1.2, ls="--", label="full-sample estimate")
    ec = "black" if eco == "LINCS" else "black"
    if eco == "DepMap":
        ax.axvline(0, color="r", lw=0.9, ls=":")
    else:
        ax.axvline(0, color="r", lw=0.9, ls=":")
    ax.set_yticks(ys); ax.set_yticklabels(lab, fontsize=7)
    ax.set_xlabel("estimate (95% CI)"); ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
fig.tight_layout()
p = os.path.join(FIG, "big_experiment_independent_cohorts.png")
fig.savefig(p, dpi=150, bbox_inches="tight")
plt.close(fig)
log("saved figure -> %s" % p)
log("saved -> %s | %.0fs" % (os.path.join(RES, "big_experiment_independent_cohorts.json"), time.time() - _t0))
