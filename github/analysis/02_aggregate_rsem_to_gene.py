# 02 RSEM transcript -> gene 聚合（最终版 v2：只收 GCT 交集基因，精确预分配）
#   pass1: 扫 GCT 基因名（find 切片，不碰数值列，~2 分钟）
#   pass2: zlib + np.loadtxt 解析 RSEM，np.add.at 只累加 GCT 内的基因
import os, sys, time, json, gzip, zlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, log, p, iter_tsv_gz_blocks

cfg = load_config()
EXP_ROWS = cfg["expected"]["transcripts_rows"]
EXP_COLS = cfg["expected"]["transcripts_sample_cols"]
STRIP = cfg["analysis"]["strip_gene_version"]


def strip_v(g):
    return g.split(".")[0] if STRIP else g


# --- pass1: GCT 基因集 ---
log("pass1: 扫 GCT 基因集")
d = zlib.decompressobj(16 + zlib.MAX_WBITS)
gct_set = set()
t0 = time.time()
with open(cfg["paths"]["gene_tpm_gz"], "rb") as f:
    buf = b""
    while True:
        b = f.read(64 << 20)
        buf += d.decompress(b) if b else d.flush()
        if not b and not buf:
            break
        lines = buf.split(b"\n")
        buf = lines.pop()
        for ln in lines:
            t = ln.find(b"\t")
            if t <= 0:
                continue
            gct_set.add(strip_v(ln[:t].decode("latin1")))
        if not b:
            break
log("GCT 基因 %d 个, %.0fs" % (len(gct_set), time.time()-t0))
MAX_G = len(gct_set)

with gzip.open(cfg["paths"]["transcripts_tpm_gz"], "rt", encoding="utf-8") as f:
    samples = f.readline().rstrip("\r\n").split("\t")[2:]
assert len(samples) == EXP_COLS

# --- pass2: 聚合 ---
gene_pos = {}
data = np.zeros((MAX_G, EXP_COLS), np.float32)
n = 0
t0 = time.time(); nrows = 0; nblk = 0
for gid, vals in iter_tsv_gz_blocks(cfg["paths"]["transcripts_tpm_gz"],
                                    skip_lines=1, n_id_cols=2, id_index=1):
    gid = [strip_v(g) for g in gid]
    mask = np.fromiter((g in gct_set for g in gid), bool, len(gid))
    gid_f = [g for g, m in zip(gid, mask) if m]
    vals_f = vals[mask]
    uniq, inv = np.unique(np.asarray(gid_f), return_inverse=True)
    loc = np.empty(len(uniq), np.int64)
    for i, u in enumerate(uniq):
        j = gene_pos.get(u)
        if j is None:
            j = n; gene_pos[u] = j; n += 1
            assert n <= MAX_G
        loc[i] = j
    out = np.zeros((len(uniq), EXP_COLS), np.float32)
    np.add.at(out, inv, vals_f)
    data[loc] += out
    nrows += len(gid); nblk += 1
    el = time.time()-t0
    log("  块%d: %d/%d 行 (%.0f%%), 收录基因 %d, %.0fs, ETA %.0fs"
        % (nblk, nrows, EXP_ROWS, 100*nrows/EXP_ROWS, n, el, el/nrows*(EXP_ROWS-nrows)))

data = data[:n]
log("完成：%d 转录本 -> %d 基因(GCT∩) × %d 样本, %.0fs" % (nrows, n, EXP_COLS, time.time()-t0))
assert nrows == EXP_ROWS, "行数 %d != 期望 %d（可能文件截断）" % (nrows, EXP_ROWS)

os.makedirs(p("data_processed"), exist_ok=True)
np.save(p("data_processed", "RSEM_gene_tpm.npy"), data)
json.dump({"gene_ids": list(gene_pos.keys()), "samples": samples},
          open(p("data_processed", "RSEM_gene_tpm.meta.json"), "w"))
log("写出 data_processed/RSEM_gene_tpm.npy (%.2f GB) + meta.json" % (data.nbytes/1e9))
