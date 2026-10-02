# 00 检查两个原始 gz 的表头与维度，核对 config.expected；产出 metadata/headers.json
import os, sys, zlib, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, log, p

cfg = load_config()
exp = cfg["expected"]


def head_lines(gz_path, nbytes=4_000_000, n_lines=12):
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    with open(gz_path, "rb") as f:
        txt = d.decompress(f.read(nbytes)).decode("utf-8", "replace")
    return txt.split("\n")[:n_lines]


log("检查 gene_tpm GCT 表头")
gct_lines = head_lines(cfg["paths"]["gene_tpm_gz"])
log("  line0: %r" % gct_lines[0])
log("  line1: %r" % gct_lines[1])          # "74628\t19788"
gct_rows, gct_cols = (int(x) for x in gct_lines[1].split("\t"))
colnames = gct_lines[2].split("\t")
log("  表头行(line2) 前列: %s ...  共 %d 列" % (colnames[:3], len(colnames)))
assert gct_rows == exp["gct_rows"], "GCT 基因行数 %d != 期望 %d" % (gct_rows, exp["gct_rows"])
assert gct_cols == exp["gct_cols"], "GCT 样本列数 %d != 期望 %d" % (gct_cols, exp["gct_cols"])
assert len(colnames) == gct_cols + 2, "GCT 表头列数 %d != %d + 2(Name/Description)" % (len(colnames), gct_cols)

log("检查 transcripts_tpm 表头")
tr_lines = head_lines(cfg["paths"]["transcripts_tpm_gz"])
tr_cols = tr_lines[0].split("\t")
log("  首行 %d 列: %s ..." % (len(tr_cols), tr_cols[:5]))
assert tr_cols[0] == "transcript_id" and tr_cols[1] == "gene_id", \
    "transcripts 首两列不是 transcript_id/gene_id: %s" % tr_cols[:2]
n_tr_samples = len(tr_cols) - 2
assert n_tr_samples == exp["transcripts_sample_cols"], \
    "transcripts 样本列 %d != 期望 %d" % (n_tr_samples, exp["transcripts_sample_cols"])

# 样本顺序与 GCT 是否一致（已验收过一次，这里再快查前 50 个）
gct_samples = colnames[2:]
if gct_samples == tr_cols[2:]:
    sample_order_identical = True
else:
    sample_order_identical = set(gct_samples) == set(tr_cols[2:])
    log("  ⚠ 样本顺序不一致但集合一致，03 会按 GCT 顺序重排")
assert sample_order_identical, "两个文件的样本集合不一致，中止"

out = {"gct_rows": gct_rows, "gct_cols": gct_cols,
       "gct_sample_columns": gct_samples,
       "transcripts_rows_expected": exp["transcripts_rows"],
       "sample_order_identical_to_gct": bool(gct_samples == tr_cols[2:])}
with open(p("metadata", "headers.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
log("PASS：表头/维度核对全部通过 -> metadata/headers.json")
