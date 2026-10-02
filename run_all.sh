#!/usr/bin/env bash
# 复现入口：按顺序重算全部分析并校验关键数字
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/6] DepMap 共同基因背景对照（表 4、S4–S6）"
python github/analysis/depmap_background_control.py --data zenodo/source_data

echo "[2/6] DepMap 阳性对照（S13–S15）"
python github/analysis/positive_control_expression_ridge.py --data zenodo/source_data/depmap_positive_control

echo "[3/6] GTEx 参考匹配（S8、S10–S12）"
for s in 06_score_2x2 07_bootstrap; do python github/analysis/${s}.py; done

echo "[4/6] LINCS 排序稳定性 / crossed 自助（S17–S22）"
python github/analysis/stability_checks/_within_ref_stability.py
python github/analysis/stability_checks/_crossed_bootstrap_and_disease.py
python github/analysis/stability_checks/_patch_A_B.py

echo "[5/6] DepMap 272 细胞一致性（S16）"
python github/analysis/stability_checks/_recompute_depmap272.py

echo "[6/6] 校验关键数字"
python github/analysis/stability_checks/_audit_numbers.py --expected EXPECTED_OUTPUTS.json

echo "全部完成：若末步无 FAIL，则论文关键数字全部复现。"
