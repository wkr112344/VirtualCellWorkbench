# 共享工具：配置加载、donor 解析、日志
import os, time, yaml
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # gtex_reference_sensitivity/

def load_config():
    with open(os.path.join(BASE, "config.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)

def donor_id(sample_id: str) -> str:
    """GTEX-1117F-0226-SM-5GZZ7 -> GTEX-1117F"""
    return "-".join(sample_id.split("-")[:2])

def log(msg: str):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)

def p(*parts) -> str:
    return os.path.join(BASE, *parts)


def iter_numeric_csv_chunks(path, id_cols, block_size=1 << 29):
    """pyarrow CSV 流式解析（多线程解析 + arrow 原生 float32 cast + 零拷贝取列）。
    yield (gid_np_str, vals_float32(nrows×ncols), sample_names|None)。
    id_cols: 前几个 ID 列名；gid 取第 2 个（GCT 用 Name，RSEM 用 gene_id）。"""
    import pyarrow as pa
    import pyarrow.csv as pac
    ro = pac.ReadOptions(block_size=block_size)
    po = pac.ParseOptions(delimiter="\t")
    with pac.open_csv(path, read_options=ro, parse_options=po) as reader:
        names = reader.schema.names
        assert names[:len(id_cols)] == list(id_cols), names[:4]
        sample_names = names[len(id_cols):]
        # 只把数值列 cast 成 float32；ID 列保持 string
        target = pa.schema([(nm, pa.string() if nm in id_cols else pa.float32())
                            for nm in names])
        for batch in reader:
            tbl = pa.Table.from_batches([batch]).cast(target)
            gid = tbl.column(id_cols[1]).to_numpy(zero_copy_only=False)
            cols = [tbl.column(i).combine_chunks() for i in range(len(id_cols), len(names))]
            try:
                arrs = [c.to_numpy(zero_copy_only=True) for c in cols]
            except pa.ArrowInvalid:          # 有 null 时退回拷贝
                arrs = [c.to_numpy(zero_copy_only=False) for c in cols]
            vals = np.stack(arrs, axis=1)
            yield gid, vals, sample_names


def iter_tsv_gz_blocks(path, skip_lines=1, n_id_cols=2, id_index=1, block_comp=64 << 20):
    """zlib 流式解压 + 分块 np.loadtxt 解析（numpy 2.x 的 C 实现，~73MB/s，无逐列开销）。
    yield (ids:list[str], vals:float32(nrows×ncols))。
    skip_lines: 跳过文件头行数；n_id_cols: 前 n 个 ID 列；id_index: 用第几个 ID 列做 ids。"""
    import io
    import zlib
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    skip = skip_lines
    f = open(path, "rb")
    buf = ""
    while True:
        b = f.read(block_comp)
        buf += d.decompress(b).decode("utf-8", "replace") if b else d.flush().decode("utf-8", "replace")
        lines = buf.split("\n")
        buf = lines.pop()                      # 块尾半行留到下一轮
        proc = []
        for l in lines:
            if skip > 0:
                skip -= 1
                continue
            if not l:
                continue
            proc.append(l.split("\t", n_id_cols))
        if proc:
            ids = [pp[id_index] for pp in proc]
            rest = "\n".join(pp[-1] for pp in proc)
            vals = np.loadtxt(io.StringIO(rest), dtype=np.float32, delimiter="\t")
            if vals.ndim == 1:
                vals = vals.reshape(1, -1)
            yield ids, vals
        if not b:
            break
    f.close()
