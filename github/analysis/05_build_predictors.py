# 05 tissue-mean predictor（每 product 各建一个，只用 train 样本）
# 定义（预注册）：在 log2(TPM+1) 空间内，对每个 train 样本按 tissue_label 分组求基因均值；
# test 样本的预测 = 其组织对应列。train 中缺失组织 -> 全 train 全局均值兜底。
# 输出 data_processed/tissue_mean_{rnaseqc,rsem}.npz（tissues × genes float32）
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, log, p

cfg = load_config()
lab = cfg["analysis"]["tissue_label"]
grid = json.load(open(p("metadata", "common_grid.json")))
samples, NG = grid["samples"], grid["n_genes"]
splits = json.load(open(p("metadata", "splits.json")))
tmap = json.load(open(p("metadata", "sample_tissue_map.json")))

tr_s = splits["train_samples"]
tr_idx = np.array([samples.index(s) for s in tr_s], dtype=np.int64)
tissues = sorted({tmap[s][lab] for s in tr_s})
ti = {t: k for k, t in enumerate(tissues)}
lab_of_sample = np.array([ti[tmap[s][lab]] for s in tr_s], dtype=np.int64)
NT = len(tissues)

for prod, fn in [("rnaseqc", "rnaseqc.npy"), ("rsem", "rsem_gene.npy")]:
    M = np.load(p("data_processed", fn), mmap_mode="r")
    t0 = time.time()
    means = np.zeros((NT, NG), dtype=np.float64)
    counts = np.zeros(NT, dtype=np.int64)
    CH = 512  # 样本块
    for lo in range(0, len(tr_s), CH):
        blk = np.log2(M[tr_idx[lo:lo+CH]].astype(np.float64) + 1.0)   # (c, NG)
        for k in range(NT):
            m = lab_of_sample[lo:lo+CH] == k
            if m.any():
                means[k] += blk[m].sum(0)
                counts[k] += int(m.sum())
        log("  %s: %d/%d train 样本, %.0fs" % (prod, min(lo+CH, len(tr_s)), len(tr_s), time.time()-t0))
    means /= counts[:, None]
    # 全局均值兜底（train 中缺该组织时用）
    glob = means.sum(0) / counts.sum()
    np.savez_compressed(p("data_processed", "tissue_mean_%s.npz" % prod),
                        tissues=np.array(tissues), means=means.astype(np.float32),
                        fallback=glob.astype(np.float32), counts=counts,
                        tissue_label=lab)
    log("%s: %d 组织 × %d 基因, %.0fs -> tissue_mean_%s.npz" % (prod, NT, NG, time.time()-t0, prod))
