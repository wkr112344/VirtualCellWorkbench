# 07 donor bootstrap（预注册：重采样 test donor 5000 次，对逐样本 interaction 取中位数）
# 输出 results/bootstrap.json（CI、p 概率、显著性判定）
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, log, p

cfg = load_config()
nb = cfg["analysis"]["n_bootstrap"]

import csv
rows = list(csv.DictReader(open(p("results", "per_sample_scores.csv"), encoding="utf-8")))
log("test 样本 %d" % len(rows))

for metric in ["interaction_pearson", "interaction_spearman"]:
    by_donor = {}
    for r in rows:
        by_donor.setdefault(r["donor"], []).append(float(r[metric]))
    donors = sorted(by_donor)
    med_of = {d: float(np.median(by_donor[d])) for d in donors}
    point = float(np.median([med_of[d] for d in donors]))      # 点估计：donor 中位数的中位数
    rng = np.random.default_rng(cfg["analysis"]["seed"] + (0 if metric.endswith("pearson") else 1))
    boot = np.empty(nb)
    vals = np.array([med_of[d] for d in donors])
    for b in range(nb):
        pick = rng.integers(0, len(donors), len(donors))
        boot[b] = np.median(vals[pick])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    res = {"metric": metric, "point_median": point,
           "ci95": [float(lo), float(hi)], "excludes_zero": bool(lo * hi > 0),
           "frac_gt_zero": float((boot > 0).mean()),
           "n_donors": len(donors), "n_bootstrap": nb}
    json.dump(res, open(p("results", "bootstrap_%s.json" % metric), "w"), indent=2)
    log("%-24s point=%+.6f CI[%+.6f,%+.6f] 排零=%s frac>0=%.3f"
        % (metric, point, lo, hi, res["excludes_zero"], res["frac_gt_zero"]))

# 汇总
out = {"seed": cfg["analysis"]["seed"], "n_bootstrap": nb,
       "pearson": json.load(open(p("results", "bootstrap_interaction_pearson.json"))),
       "spearman": json.load(open(p("results", "bootstrap_interaction_spearman.json")))}
json.dump(out, open(p("results", "bootstrap.json"), "w"), indent=2)
log("PASS -> results/bootstrap.json")
