# -*- coding: utf-8 -*-
"""
S_drugsplit_metrics_ci.py —— 补出「药物轴(drug-split) cohort」的 Spearman / cosine CI

背景：P_frozen_multimetric_cohorts.py 跑的时候，药物轴 cohort 对**三口径都算过** bootstrap CI，
但写 frozen_cohort_replication.csv 时只落了 pearson 的 lo/hi；big_experiment 合并时
spearman/cosine 的 CI 被填成 0/None。本脚本把这些已算出的 CI 全部补齐，不重训、不重抽样。

产出：
  results/frozen_cohort_replication.csv            宽表，补 lo/hi/显著性（三口径）
  results/frozen_drugsplit_cohort_three_metrics.csv 长表（K × cohort × metric）
  results/big_experiment_independent_cohorts.csv   药物轴行改为真实 CI
  results/big_experiment_independent_cohorts.json 加 drug_axis_three_metrics 一致性
"""
import csv
import io
import json
import os

import numpy as np

RES = r"C:\Users\wkr20\WorkBuddy\Claw\gigascience_v5\results"
METRICS = ["pearson", "spearman", "cosine"]
METRIC_KEY = {"pearson": "delta_pearson", "spearman": "delta_spearman", "cosine": "delta_cosine"}


def log(m):
    print(m, flush=True)


with open(os.path.join(RES, "frozen_multimetric_summary.json"), encoding="utf-8") as f:
    J = json.load(f)
CR = J["cohort_replication"]

# ───────────────────────── 1. 宽表：补三口径 lo/hi ─────────────────────────
wide = []
for Kk, rec in CR.items():
    K = int(Kk[1:])
    for e in rec["per_cohort"]:
        row = {"K": K, "cohort": int(e["cohort"]), "n_pairs": e["n_pairs"], "n_drugs": e["n_drugs"]}
        for m in METRICS:
            d = e[m]
            row["delta_%s" % m] = round(d["delta"], 6)
            row["lo_%s" % m] = round(d["ci95"][0], 6)
            row["hi_%s" % m] = round(d["ci95"][1], 6)
            row["%s_excludes_zero" % m] = bool(d["excludes_zero"])
        wide.append(row)
with open(os.path.join(RES, "frozen_cohort_replication.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(wide[0].keys()))
    w.writeheader(); w.writerows(wide)
log("[written] frozen_cohort_replication.csv  (宽表, 三口径含 CI)  rows=%d" % len(wide))

# ───────────────────────── 2. 长表 ─────────────────────────
longr = []
for r in wide:
    for m in METRICS:
        longr.append(dict(split_axis="drug", K=r["K"], cohort=r["cohort"],
                          n_pairs=r["n_pairs"], n_drugs=r["n_drugs"], metric=m,
                          delta=r["delta_%s" % m], ci95_low=r["lo_%s" % m], ci95_high=r["hi_%s" % m],
                          excludes_zero=r["%s_excludes_zero" % m],
                          estimand="reference contrast Δ = score(β2020) − score(dcic2021)"))
with open(os.path.join(RES, "frozen_drugsplit_cohort_three_metrics.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(longr[0].keys()))
    w.writeheader(); w.writerows(longr)
log("[written] frozen_drugsplit_cohort_three_metrics.csv (长表)  rows=%d" % len(longr))

# ───────────────────── 3. 药物轴逐口径一致性 ─────────────────────
drug_axis_metrics = {}
for Kk, rec in CR.items():
    K = int(Kk[1:])
    per = {}
    for m in METRICS:
        est = np.array([e[m]["delta"] for e in rec["per_cohort"]], dtype=np.float64)
        lo = np.array([e[m]["ci95"][0] for e in rec["per_cohort"]], dtype=np.float64)
        hi = np.array([e[m]["ci95"][1] for e in rec["per_cohort"]], dtype=np.float64)
        se = (hi - lo) / (2 * 1.959964)
        wts = 1.0 / np.maximum(se ** 2, 1e-12)
        pooled = float((wts * est).sum() / wts.sum())
        Q = float((wts * (est - pooled) ** 2).sum())
        dfree = max(len(est) - 1, 1)
        I2 = float(max(0.0, (Q - dfree) / max(Q, 1e-12)) * 100.0)
        per[m] = {"n_cohorts": int(len(est)),
                  "estimates": [round(float(x), 6) for x in est],
                  "median": round(float(np.median(est)), 6),
                  "range": [round(float(est.min()), 6), round(float(est.max()), 6)],
                  "all_same_sign": bool(np.all(est > 0) or np.all(est < 0)),
                  "all_ci_exclude_zero": bool(all(e[m]["excludes_zero"] for e in rec["per_cohort"])),
                  "n_ci_exclude_zero": int(sum(e[m]["excludes_zero"] for e in rec["per_cohort"])),
                  "pooled": round(pooled, 6), "Q": round(Q, 4), "I2_pct": round(I2, 2)}
    drug_axis_metrics["K%d" % K] = per

# ───────────────────── 4. 修正 big_experiment CSV 的药物轴行 ─────────────────────
p_big = os.path.join(RES, "big_experiment_independent_cohorts.csv")
with open(p_big, encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
fixed = 0
for r in rows:
    if r["ecosystem"] == "LINCS" and r["split_axis"] == "drug":
        K = int(r["K"]); c = int(r["cohort"]); m = r["metric"]
        src = CR["K%d" % K]["per_cohort"][c][m]
        r["estimate"] = round(src["delta"], 6)
        r["ci95_low"] = round(src["ci95"][0], 6)
        r["ci95_high"] = round(src["ci95"][1], 6)
        r["excludes_zero"] = str(bool(src["excludes_zero"]))
        fixed += 1
with open(p_big, "w", newline="", encoding="utf-8") as f:
    keys = ["ecosystem", "split_axis", "K", "cohort", "metric", "n_pairs", "n_drugs",
            "estimate", "ci95_low", "ci95_high", "excludes_zero", "estimand"]
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
log("[fixed] big_experiment_independent_cohorts.csv  药物轴行 %d 行改为真实 CI" % fixed)

# ───────────────────── 5. 更新 big JSON ─────────────────────
with open(os.path.join(RES, "big_experiment_independent_cohorts.json"), encoding="utf-8") as f:
    B = json.load(f)
B["drug_axis_three_metrics"] = drug_axis_metrics
B.pop("lincs_note_drug_axis_only_pearson_ci", None)
B["notes"] = ("药物轴(2037 留出药切互不相交 cohort)的 Pearson/Spearman/cosine 三口径 CI 均已补出；"
              "见 frozen_drugsplit_cohort_three_metrics.csv 与 frozen_cohort_replication.csv。")
# 重算 by_view 的 ALL（pearson 口径不受影响，但 excludes_zero 现在更完整）
pear = [r for r in rows if r["metric"] == "pearson"]
seen, one = set(), []
for r in pear:
    k = (r["ecosystem"], r["split_axis"], r["K"], r["cohort"])
    if k not in seen:
        seen.add(k); one.append(r)


def counts(sub):
    if not sub:
        return dict(n=0)
    e = np.array([float(r["estimate"]) for r in sub], dtype=np.float64)
    ex = [r["excludes_zero"] == "True" for r in sub]
    return dict(n=len(sub), same_sign=bool(np.all(e > 0) or np.all(e < 0)),
                all_excl=bool(all(ex)), n_excl=int(sum(ex)),
                median=round(float(np.median(e)), 6),
                range=[round(float(e.min()), 6), round(float(e.max()), 6)])


B["by_view"] = {
    "LINCS_drug_axis": counts([r for r in one if r["ecosystem"] == "LINCS" and r["split_axis"] == "drug"]),
    "LINCS_cell_axis": counts([r for r in one if r["ecosystem"] == "LINCS" and r["split_axis"] == "cell_line"]),
    "DepMap_cell_axis": counts([r for r in one if r["ecosystem"] == "DepMap"]),
    "ALL": counts(one)}
json.dump(B, io.open(os.path.join(RES, "big_experiment_independent_cohorts.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
log("[written] big_experiment_independent_cohorts.json + drug_axis_three_metrics")

log("\n药物轴三口径汇总（每 cohort 独立估计 + 药物聚类 bootstrap 2000, seed 12345）")
for Kk, per in drug_axis_metrics.items():
    log("  %s:" % Kk)
    for m in METRICS:
        d = per[m]
        log("    %-9s median=%+.4f range=[%+.4f,%+.4f] 同号=%s CI全不跨0=%s(%d/%d) pooled=%+.4f I²=%.1f%%"
            % (m, d["median"], d["range"][0], d["range"][1], d["all_same_sign"],
               d["all_ci_exclude_zero"], d["n_ci_exclude_zero"], d["n_cohorts"], d["pooled"], d["I2_pct"]))
