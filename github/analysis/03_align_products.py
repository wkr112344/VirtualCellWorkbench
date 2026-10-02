# 03 对齐两套 product（快速版）
#   pass1: 裸 gzip 流 + 有限 split 只取基因名列，确定共同基因空间（快）
#   pass2: pandas C 引擎分块读 GCT（header=2），按块转置写入 memmap
#   RSEM:  npy(memmap) 按行序分块映射到共同基因列
# 输出 memmap (n_samples × n_genes) float32：data_processed/{rnaseqc,rsem_gene}.npy
import os, sys, json, time, zlib
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, log, p, iter_numeric_csv_chunks

cfg = load_config()
hdr = json.load(open(p("metadata", "headers.json")))
samples = hdr["gct_sample_columns"]
S = len(samples)
EXP_GCT_ROWS = cfg["expected"]["gct_rows"]


def strip_v(g):
    return g.split(".")[0] if cfg["analysis"]["strip_gene_version"] else g


# --- RSEM gene 矩阵（memmap） ---
log("打开 RSEM gene 矩阵")
meta = json.load(open(p("data_processed", "RSEM_gene_tpm.meta.json")))
R = np.load(p("data_processed", "RSEM_gene_tpm.npy"), mmap_mode="r")
rsem_genes_raw = meta["gene_ids"]
rsem_samples = meta["samples"]
assert set(rsem_samples) == set(samples), "RSEM 与 GCT 样本集合不一致"
B_samplecols = np.array([rsem_samples.index(s) for s in samples], dtype=np.int64)
rsem_row = {strip_v(g): i for i, g in enumerate(rsem_genes_raw)}
log("RSEM: %d 基因 × %d 样本 (memmap %.2f GB)" % (R.shape[0], R.shape[1], R.shape[0]*R.shape[1]*4/1e9))

# --- pass1 已由 02 完成：直接复用其基因列表（= GCT∩RSEM，74,628） ---
common = list(rsem_row.keys())   # 02 的收录顺序（RSEM 文件序），A/B 两 memmap 列序一致即可
NG = len(common)
log("共同基因 %d 个（复用 02 结果，免重扫）" % NG)

A = np.lib.format.open_memmap(p("data_processed", "rnaseqc.npy"), mode="w+",
                              dtype=np.float32, shape=(S, NG))
Bm = np.lib.format.open_memmap(p("data_processed", "rsem_gene.npy"), mode="w+",
                               dtype=np.float32, shape=(S, NG))

# --- pass2: zlib+loadtxt 流式读 GCT，写 A ---
log("pass2: 流式读 GCT -> rnaseqc.npy")
from common import iter_tsv_gz_blocks
gidx = {g: j for j, g in enumerate(common)}
t0 = time.time()
done = 0
for gid, vals in iter_tsv_gz_blocks(cfg["paths"]["gene_tpm_gz"],
                                    skip_lines=3, n_id_cols=2, id_index=0):
    gid = [strip_v(g) for g in gid]
    mask = np.fromiter((g in gidx for g in gid), bool, len(gid))
    if mask.any():
        cols = np.fromiter((gidx[g] for g in gid), np.int64, len(gid))[mask]
        A[:, cols] = vals[mask].T
    done += len(gid)
    log("  GCT %d/%d 行, %.0fs" % (done, EXP_GCT_ROWS, time.time()-t0))
A.flush()

# --- RSEM 分块映射（按 RSEM 行序收集列，块转置写入） ---
log("RSEM -> rsem_gene.npy（分块转置）")
order = [(i, gidx[g]) for g, i in ((g, rsem_row[g]) for g in common)]  # (rsem_row, out_col)
t0 = time.time()
K = 2000
for lo in range(0, len(order), K):
    blk = order[lo:lo+K]
    rrows = np.array([x[0] for x in blk], dtype=np.int64)
    ccols = np.array([x[1] for x in blk], dtype=np.int64)
    Bm[:, ccols] = np.asarray(R[rrows], dtype=np.float32)[:, B_samplecols].T
    log("  %d/%d, %.0fs" % (min(lo+K, len(order)), len(order), time.time()-t0))
Bm.flush()

json.dump({"genes": common, "samples": samples, "n_samples": S, "n_genes": NG},
          open(p("metadata", "common_grid.json"), "w"), indent=1)
log("PASS -> data_processed/{rnaseqc,rsem_gene}.npy (%d × %d) + metadata/common_grid.json" % (S, NG))
