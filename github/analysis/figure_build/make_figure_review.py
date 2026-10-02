"""Rebuild the author figure review for the final 7-figure set."""
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
    ("Figure 1  Study design",
     "Three data systems each compare the reference products of drug-perturbation transcriptomes, cancer dependencies and normal-tissue expression quantification. "
     "Within each system the common evaluation objects and gene space are held fixed. DepMap frozen output and controlled training are different branches; "
     "GTEx splits training and test by donor. The three branches address complementary questions and cannot be merged into independent repeats of one task."),
    ("Figure 2  LINCS fixed-output sensitivity: distribution -> scatter -> exemplars -> multi-model",
     "A: reference-vs-reference density scatter for 2,037 drugs (hexbin + marginal distributions; the diagonal = no difference; "
     "2031/2037 drugs shift in the same direction). "
     "B: per-drug Δ ridgelines aligned across three calibers (Pearson/Spearman/cosine, each annotated with its median and positive-difference fraction). "
     "C: exemplar-drug small multiples -- per-cell two-reference score scatters for the 3 most sensitive and 3 relatively stable drugs. "
     "D: within-evaluation-set differences for the seven prediction sources (own - dcic2021) with paired 95% CIs."),
    ("Figure 3  Training x reference interaction and drug-level reversal",
     "A: drug-level distributions for the four states (2x2; medians 0.353 / 0.054 / 0.155 / 0.134). "
     "B: drug-level reversal quadrant -- x = the beta-trained advantage under the beta reference, y = under the dcic reference; "
     "the reversal fraction is annotated in the panel. C: per-drug interaction distribution [(bb - bD) - (Db - DD)], "
     "vertical line = overall interaction 0.2785 [0.2741, 0.2828]. "
     "D: score collapse of the 5 drugs with the largest rank drop within the beta top-50 (the number after the arrow is the rank under dcic2021)."),
    ("Figure 4  Ranking consequences and biological priority",
     "A: top-k retention and overlapping candidates (solid = observed, dashed = random k/N; only 3/10 overlap at k = 10). "
     "B: rank flow of the union of the two references' top-10 candidates: blue band = the 3 shared drugs, grey/orange bands = the 7 unique to each list. "
     "C: shortlist membership matrix -- drug ranks within the top-10 of Global and three disease contexts (Liver/Thyroid ATC/APL); "
     "colour encodes rank (log10 scale, grey = not in list; HALLMARK lists are not in the package and are not shown). "
     "D: summary of rank Spearman and top-10 overlap across the 10 contexts. E: agreement vs list-stability scatter."),
    ("Figure 5  LINCS heterogeneity: landscape -> ranking -> structure -> distribution -> exemplars",
     "A: ordered landscape of 63 cell lines x [beta2020, dcic2021, Δ, response count] (sorted by descending Δ, colour scale clipped at the 97th percentile). "
     "B: ranking caterpillar -- the mean Δ per cell line (point size proportional to response count, dashed line = equal-weight mean 0.277). "
     "C: drug x cell Δ matrix (transposed: rows = 63 cell lines, columns = the 36 drugs covering >=25 cells, columns ordered by hierarchical clustering, "
     "every 6th column labelled with a drug index; the full list is in TableS_drug_coverage_matrix_drugs.csv). "
     "D: unit-level distributions (per-drug and per-cell aligned ridgelines, medians 0.294 / 0.240). "
     "E: exemplar-cell distributions (the 2 most sensitive + the 2 most stable, n >= 20)."),
    ("Figure 6  DepMap mechanism + positive control: raw background -> residualization -> signal resolution",
     "A: Frozen model raw per-cell distributions (four states, 136 cells). "
     "B: paired difference between the model and the training-set gene-mean baseline (~0, showing the raw score is driven by a stable gene background). "
     "C: Frozen residualized per-cell distributions (same layout as A, collapsing to ~0). "
     "D: per-gene across-cell r distributions for 17,393 genes (Expr Ridge / Frozen / Permuted). "
     "E: positive-control per-cell residual distributions (Expr Ridge retains 0.216/0.277, Frozen and Permuted ~0). "
     "F: bootstrap distributions of the interaction from raw to residualized (Frozen with 10,000 draws, "
     "Ridge/Permuted with per-cell resampling; +0.0563 -> -0.0005 and +0.0620 -> +0.3057)."),
    ("Figure 7  GTEx cross-ecosystem generalization: distribution -> donor -> summary -> uncertainty -> tissue structure",
     "A: four-state ridgeline distributions of per-sample across-gene Pearson for 3,963 test samples (RR/RE/ER/EE). "
     "B: the same as A in the Spearman caliber. "
     "C: same-origin reference advantage raincloud over 189 donors (median/IQR/positive-direction fraction). "
     "D: compact 2x2 summary (interaction = 0.07096). "
     "E: interaction distribution from 5,000 donor bootstraps (Pearson 0.07096 [0.07050, 0.07164]; "
     "Spearman 0.06186 [0.06132, 0.06224]). "
     "F: 68-tissue x 5-column summary matrix (blocks coloured by organ system, ordered within a group by the Pearson interaction; "
     "columns: Pearson/Spearman interaction, donor count, the two references' matched advantage)."),
]


CHECKLIST = [
    ("Main figures compressed from 8 to 7 (evidence re-ordering)", "Each main figure answers one complete scientific question; purely summary/CI panels are moved out of the main figures."),
    ("Observation-level data becomes the main visual", "2,037 drugs / 63 cells / 136 test cells / 17,393 genes / 3,963 samples / 189 donors / 68 tissues are shown directly."),
    ("dcic-trained per-drug scores recovered", "recomputed in place from the archived frozen prediction matrices + the original reference cache (11,275 rows; the max difference of per-drug means vs C_perdrug_stability is 6.5e-07, "
     "matching the TableS25 aggregates 0.1705/0.1534) -- so Figure 3's distributions/quadrant/interaction distribution use real data."),
    ("Fig 6 merges mechanism and positive control", "Raw/Residualized two columns + a unified visual grammar for the three models Frozen/Expr Ridge/Permuted, closing the loop over six panels."),
    ("Fig 7 returns to distribution first", "Four-state ridgeline (Pearson/Spearman two columns) + raincloud + compact 2x2 + bootstrap violin + the 68-tissue x 5-column structure matrix."),
    ("Figure 1 kept as the original", "sha256 is bit-identical to the supplied file; not redrawn."),
    ("Main text must be synced", "All captions of Figures 2-7 (new captions appended at the end) + figure-number references (original Figures 6/7/8 -> new 5/6/7)."),
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
    h = doc.add_paragraph(); r = h.add_run("Figure review (after main-figure evidence re-ordering, 7 figures)")
    r.font.size = Pt(16); r.font.bold = True
    sub = doc.add_paragraph(); rs = sub.add_run(
        "Completed the main-figure evidence re-ordering as instructed: 8 figures compressed to 7, each organized as "
        "\"overview -> structure -> consequence/exemplar\", with purely summary/CI panels moved out of the main figures. "
        "New captions for Figures 2-7 are at the end; the main-text captions and figure-number references must be synced.")
    rs.font.size = Pt(9.5)
    for name, (title, cap) in zip(FILES, CAPS):
        p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(14)
        rr = p.add_run(title); rr.font.bold = True; rr.font.size = Pt(11.5)
        img = doc.add_paragraph(); img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        img.add_run().add_picture(str(F / (name + ".png")), width=Cm(16.5))
        cp = doc.add_paragraph(); cr = cp.add_run(disp(cap)); cr.font.size = Pt(9.5)
    doc.add_page_break()
    p = doc.add_paragraph(); r = p.add_run("Re-ordering checklist for this round")
    r.font.bold = True; r.font.size = Pt(13)
    for k, v in CHECKLIST:
        para = doc.add_paragraph(style="List Bullet")
        a = para.add_run("√ " + k + "："); a.font.bold = True; a.font.size = Pt(10)
        b = para.add_run(v); b.font.size = Pt(10)
    doc.add_paragraph()
    p = doc.add_paragraph(); r = p.add_run("New captions for Figures 2-7 (ready to replace the main text)")
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
        fig.text(0.07, 0.94, "Figure review (after main-figure evidence re-ordering, 7 figures)", fontsize=19, weight="bold")
        fig.text(0.07, 0.90, "Figure 1-7 - for GigaScience submission", fontsize=12, color="#444444")
        body = ("Completed the main-figure evidence re-ordering as instructed, 8 figures compressed to 7:\n\n"
                "- Figure 1: study design (unchanged).\n"
                "- Figure 2: LINCS fixed-output sensitivity -- distribution / scatter / exemplars / multi-model.\n"
                "- Figure 3: winner reversal and ranking consequences -- 2x2 / crossing / rank flow / top-k.\n"
                "- Figure 4: LINCS heterogeneity and structure -- landscape / cell effects / context heatmap / ECDF.\n"
                "- Figure 5: DepMap mechanism (original Figure 6).\n"
                "- Figure 6: DepMap positive control (original Figure 8).\n"
                "- Figure 7: GTEx generalization (unchanged).\n\n"
                "File names synced: Figure4_heterogeneity / Figure5_depmap_mechanism /\n"
                "Figure6_positive_control; the old files are retired.\n"
                "The main-text captions and figure-number references must be synced (new captions appended at the end).\n"
                "This PDF matches author_notes/Figure_Review_CN.docx.")
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
        fig.text(0.06, 0.965, "Re-ordering checklist for this round", fontsize=15, weight="bold", va="top")
        y = 0.925
        for k, v in CHECKLIST:
            fig.text(0.06, y, "√ " + k, fontsize=10.5, weight="bold", va="top")
            y -= 0.028
            txt = wrap_cjk(v, 60)
            fig.text(0.08, y, txt, fontsize=9, va="top", linespacing=1.55, color="#333333")
            y -= 0.028 * (txt.count("\n") + 1) + 0.018
        pdf.savefig(fig); plt.close(fig)

        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, "New captions for Figures 2-7 (ready to replace the main text)", fontsize=15, weight="bold", va="top")
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
