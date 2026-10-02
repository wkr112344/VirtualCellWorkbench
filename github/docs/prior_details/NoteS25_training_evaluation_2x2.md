# Supplementary Note S25. Training-product × evaluation-product matched 2×2

## S25.1 正文口径：与 0.3680 / 0.0724 对齐的 11,275 对

为分析 training target 与 evaluation reference 的匹配关系，我们在正文主分析相同的 11,275 个共同 `(cell line, drug)` 对、同一 978 landmark-gene 面板和同一评分协议下构建 2×2。两个 prediction states 分别为 beta-trained 和 dcic-trained，并分别在 beta 与 dcic2021 reference 上评分。药物聚类 bootstrap 以 2,037 个 held-out drugs 为重采样单位。

| Training product | beta-eval | dcic2021-eval | beta − dcic contrast |
|---|---:|---:|---:|
| beta-trained | 0.3680 [0.3627, 0.3733] | 0.0724 [0.0683, 0.0765] | +0.2956 [0.2917, 0.2996] |
| dcic-trained | 0.1705 [0.1652, 0.1757] | 0.1534 [0.1485, 0.1581] | +0.0171 [0.0140, 0.0204] |

Difference-in-differences interaction：

`I = (M_BB − M_BD) − (M_DB − M_DD) = +0.2785 [0.2741, 0.2828]`.

beta-trained prediction 在两套 reference 上的差值明显大于 dcic-trained prediction。模型相对顺序也随 evaluation product 改变：beta-eval 上 beta-trained > dcic-trained，而 dcic-eval 上 dcic-trained > beta-trained。该结果表明，参考产品变化对 benchmark reading 的影响与 training target 密切相关。

这一 2×2 估计的是训练来源与评价参考共同作用下的分数差异。cross-reference shift 仍可能同时包含 processing-lineage、measurement noise、dynamic range 和 processed construct reliability 等因素；由于缺少独立的 within-product test-retest baseline，本分析不进一步分解这些来源。

## S25.2 更大共同集合敏感性：E_common = 15,990

该口径覆盖 112 个细胞系，用于检验正文结构是否依赖 11,275 对主集合。

| Training product | beta-eval | dcic2021-eval | beta − dcic contrast |
|---|---:|---:|---:|
| beta-trained | 0.3438 | 0.1124 | +0.2314 |
| dcic-trained | 0.1792 | 0.1390 | +0.0402 |

对应 interaction = +0.1912。更大共同集合与正文口径得到相同的交互方向，但绝对数值不同，因此两套四格结果分别报告。
