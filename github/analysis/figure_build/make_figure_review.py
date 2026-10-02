"""Rebuild the author figure review (Chinese) for the final 7-figure set."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

R = Path(__file__).resolve().parents[1]
F, A, P = R / "figures", R / "author_notes", R / "preview"
CJK = "Microsoft YaHei"

FILES = ["Figure1_design", "Figure2_fixed_output", "Figure3_training_evaluation",
         "Figure4_ranking_consequence", "Figure5_lincs_heterogeneity",
         "Figure6_depmap_mechanism_control", "Figure7_gtex_reference_matching"]

CAPS = [
    ("图 1　研究设计",
     "三个数据体系分别比较药物扰动转录组、癌症依赖和正常组织表达定量的参考产品。"
     "每个体系内固定共同评价对象与基因空间。DepMap 冻结输出与受控训练是不同分支；"
     "GTEx 按供者划分训练与测试。三个分支支持互补问题，不能合并为同一任务的独立重复。"),
    ("图 2　LINCS 固定输出敏感性：分布 → 散点 → 示例 → 多模型",
     "A：2,037 个药物的 reference-vs-reference 密度散点（hexbin + 边缘分布；对角线 = 无差异；"
     "2031/2037 个药物同方向偏移）。"
     "B：三种口径的逐药 Δ 对齐山脊分布（Pearson/Spearman/cosine，各标中位数与正差值比例）。"
     "C：示例药物小倍数图——最敏感 3 个与相对稳定 3 个的逐细胞两参考分数散点。"
     "D：七个预测来源各自评价集合内的差值（own − dcic2021）及配对 95% CI。"),
    ("图 3　训练 × 参考交互与药物级反转",
     "A：四个状态的药物级分布（2×2；中位数 0.353 / 0.054 / 0.155 / 0.134）。"
     "B：药物级反转象限——x = beta-trained 在 β 参考下的优势，y = 在 dcic 参考下的优势；"
     "反转比例见图内标注。C：逐药交互分布 [(ββ − βD) − (Dβ − DD)]，"
     "竖线 = 总体交互 0.2785 [0.2741, 0.2828]。"
     "D：beta top-50 内排名跌幅最大的 5 个药物的得分坍缩（→ 后为 dcic2021 下的排名）。"),
    ("图 4　排序后果与生物学优先级",
     "A：top-k 保留率与重合候选数（实线 = 观察，虚线 = 随机 k/N；k = 10 时仅 3/10 重合）。"
     "B：两参考各自前 10 候选的并集 rank flow：蓝带 = 共享的 3 个药物，灰/橙带 = 仅单一名单的各 7 个。"
     "C：shortlist 成员矩阵——Global 与三个疾病上下文（Liver/Thyroid ATC/APL）各自 top-10 的药物排名；"
     "颜色编码排名（log10 色标，灰 = 未入名单；HALLMARK 名单不在包内故未列）。"
     "D：10 个上下文的排序 Spearman 与 top-10 重合数汇总。E：一致性 vs 名单稳定性散点。"),
    ("图 5　LINCS 异质性：景观 → 排序 → 结构 → 分布 → 示例",
     "A：63 个细胞系 × [beta2020, dcic2021, Δ, 响应数] 有序景观（按 Δ 降序，色标 97 分位截断）。"
     "B：排序 caterpillar——每个细胞系的平均 Δ（点大小 ∝ 响应数，虚线 = 等权均值 0.277）。"
     "C：药物 × 细胞 Δ 矩阵（转置：行 = 63 个细胞系、列 = 36 个覆盖 ≥25 细胞的药物，列按层次聚类排序，"
     "每 6 列标一个药号；完整清单见 TableS_drug_coverage_matrix_drugs.csv）。"
     "D：单元级分布（per-drug 与 per-cell 对齐山脊，中位 0.294 / 0.240）。"
     "E：示例细胞分布（最敏感 2 个 + 最稳定 2 个，n ≥ 20）。"),
    ("图 6　DepMap 机制 + 阳性对照：raw 背景 → 残差化 → 信号分辨",
     "A：Frozen 模型 raw 逐细胞分布（四状态，136 细胞）。"
     "B：模型 − 训练集基因均值基线的配对差（≈0，说明 raw 分数由稳定基因背景驱动）。"
     "C：Frozen 残差化逐细胞分布（与 A 同版式，坍缩到 ≈0）。"
     "D：17,393 个基因的逐基因跨细胞 r 分布（Expr Ridge / Frozen / Permuted）。"
     "E：阳性对照逐细胞残差分布（Expr Ridge 保留 0.216/0.277，Frozen 与 Permuted ≈0）。"
     "F：交互 raw → residualized 的 bootstrap 分布（Frozen 用 10,000 次抽样，"
     "Ridge/Permuted 用逐细胞重抽样；+0.0563→−0.0005 与 +0.0620→+0.3057）。"),
    ("图 7　GTEx 跨体系泛化：分布 → 供者 → 汇总 → 不确定性 → 组织结构",
     "A：3,963 个测试样本的逐样本跨基因 Pearson 四状态山脊分布（RR/RE/ER/EE）。"
     "B：同 A 的 Spearman 口径。"
     "C：189 个供者的同源参考优势 raincloud（median/IQR/正方向比例）。"
     "D：紧凑 2×2 汇总（交互 = 0.07096）。"
     "E：5,000 次 donor bootstrap 的交互分布（Pearson 0.07096 [0.07050, 0.07164]；"
     "Spearman 0.06186 [0.06132, 0.06224]）。"
     "F：68 组织 × 5 列汇总矩阵（按器官系统分块着色、组内按 Pearson 交互排序；"
     "列：Pearson/Spearman 交互、供者数、两参考的 matched advantage）。"),
]


CHECKLIST = [
    ("主图 8 张压缩为 7 张（证据重排）", "每张主图回答一个完整科学问题；纯 summary/CI panel 移出主图。"),
    ("观察级数据成为主视觉", "2,037 drugs / 63 cells / 136 test cells / 17,393 genes / 3,963 samples / 189 donors / 68 tissues 直接展示。"),
    ("dcic-trained 逐药分数已找回", "由归档冻结预测矩阵 + 原参考 cache 就地重算（11,275 行；逐药均值与 C_perdrug_stability 最大差 6.5e-07，"
     "与 TableS25 聚合 0.1705/0.1534 吻合）——图 3 的分布/象限/交互分布因此可用真实数据。"),
    ("Fig 6 合并机制与阳性对照", "Raw/Residualized 双列 + Frozen/Expr Ridge/Permuted 三模型统一视觉语法，六 panel 闭环。"),
    ("Fig 7 回到分布优先", "四状态 ridgeline（Pearson/Spearman 双列）+ raincloud + 紧凑 2×2 + bootstrap violin + 68 组织×5 列结构矩阵。"),
    ("Figure 1 保留原版", "sha256 与来件逐位一致，未重绘。"),
    ("正文需同步", "图 2–7 全部图注（文末附新图注）+ 图编号引用（原图 6/7/8 → 新 5/6/7）。"),
]


def wrap_cjk(text, width=50):
    out, line, n = [], "", 0
    for ch in text:
        w = 2 if ord(ch) > 0x2E80 else 1
        if ch == "\n":
            out.append(line); line, n = "", 0
            continue
        if n + w > width:
            out.append(line); line, n = "", 0
        line += ch; n += w
    if line:
        out.append(line)
    return "\n".join(out)


def disp(t):
    return t.replace("10⁻³", "1e-3").replace("⁻", "^-")


def build_docx():
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = CJK
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), CJK)
    h = doc.add_paragraph(); r = h.add_run("图表审阅版（主图证据重排后，7 张）")
    r.font.size = Pt(16); r.font.bold = True
    sub = doc.add_paragraph(); rs = sub.add_run(
        "按作者指示完成主图证据重排：8 张压缩为 7 张，每张主图按"
        "「overview → structure → consequence/exemplar」组织，纯 summary/CI 移出主图。"
        "图 2–7 的新图注见文末，正文图注与图编号引用需同步。")
    rs.font.size = Pt(9.5)
    for name, (title, cap) in zip(FILES, CAPS):
        p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(14)
        rr = p.add_run(title); rr.font.bold = True; rr.font.size = Pt(11.5)
        img = doc.add_paragraph(); img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        img.add_run().add_picture(str(F / (name + ".png")), width=Cm(16.5))
        cp = doc.add_paragraph(); cr = cp.add_run(disp(cap)); cr.font.size = Pt(9.5)
    doc.add_page_break()
    p = doc.add_paragraph(); r = p.add_run("本轮证据重排清单")
    r.font.bold = True; r.font.size = Pt(13)
    for k, v in CHECKLIST:
        para = doc.add_paragraph(style="List Bullet")
        a = para.add_run("√ " + k + "："); a.font.bold = True; a.font.size = Pt(10)
        b = para.add_run(v); b.font.size = Pt(10)
    doc.add_paragraph()
    p = doc.add_paragraph(); r = p.add_run("图 2–7 新图注（可直接替换正文）")
    r.font.bold = True; r.font.size = Pt(13)
    for title, cap in CAPS[1:]:
        para = doc.add_paragraph(); rr = para.add_run(disp(cap)); rr.font.size = Pt(10)
        para.paragraph_format.space_after = Pt(8)
    doc.save(A / "Figure_Review_CN.docx")
    print("wrote", A / "Figure_Review_CN.docx")


def build_pdf():
    matplotlib.rcParams["font.family"] = CJK
    matplotlib.rcParams["axes.unicode_minus"] = False
    out = P / "Figure_Review_CN.pdf"
    with PdfPages(out) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.07, 0.94, "图表审阅版（主图证据重排后，7 张）", fontsize=19, weight="bold")
        fig.text(0.07, 0.90, "Figure 1–7 · GigaScience 投稿用", fontsize=12, color="#444444")
        body = ("按作者指示完成主图证据重排，8 张压缩为 7 张：\n\n"
                "· 图 1：study design（原样）。\n"
                "· 图 2：LINCS 固定输出敏感性——分布 / 散点 / 示例 / 多模型。\n"
                "· 图 3：获胜反转与排序后果——2×2 / 交叉 / rank flow / top-k。\n"
                "· 图 4：LINCS 异质性与结构——景观 / 细胞效应 / 上下文热图 / ECDF。\n"
                "· 图 5：DepMap 机制（原图 6）。\n"
                "· 图 6：DepMap 阳性对照（原图 8）。\n"
                "· 图 7：GTEx 泛化（原样）。\n\n"
                "文件名同步：Figure4_heterogeneity / Figure5_depmap_mechanism /\n"
                "Figure6_positive_control；旧文件退役。\n"
                "正文图注与图编号引用需同步（文末附新图注）。\n"
                "本 PDF 与 author_notes/Figure_Review_CN.docx 内容一致。")
        fig.text(0.07, 0.84, wrap_cjk(body, 46), fontsize=10.5, va="top", linespacing=1.7)
        pdf.savefig(fig); plt.close(fig)

        for name, (title, cap) in zip(FILES, CAPS):
            im = Image.open(F / (name + ".png"))
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.06, 0.965, disp(title), fontsize=12.5, weight="bold", va="top")
            ax = fig.add_axes([0.06, 0.26, 0.88, 0.68])
            ax.imshow(im); ax.axis("off")
            fig.text(0.06, 0.225, wrap_cjk(disp(cap), 56), fontsize=9.5, va="top", linespacing=1.7)
            pdf.savefig(fig); plt.close(fig)

        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, "本轮证据重排清单", fontsize=15, weight="bold", va="top")
        y = 0.925
        for k, v in CHECKLIST:
            fig.text(0.06, y, "√ " + k, fontsize=10.5, weight="bold", va="top")
            y -= 0.028
            txt = wrap_cjk(v, 60)
            fig.text(0.08, y, txt, fontsize=9, va="top", linespacing=1.55, color="#333333")
            y -= 0.028 * (txt.count("\n") + 1) + 0.018
        pdf.savefig(fig); plt.close(fig)

        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, "图 2–7 新图注（可直接替换正文）", fontsize=15, weight="bold", va="top")
        yy = 0.91
        for title, cap in CAPS[1:]:
            txt = wrap_cjk(disp(cap), 62)
            fig.text(0.06, yy, txt, fontsize=9.5, va="top", linespacing=1.75)
            yy -= 0.042 * (txt.count("\n") + 2)
        pdf.savefig(fig); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    build_docx()
    build_pdf()
