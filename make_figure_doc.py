"""
make_figure_doc.py  —  Generate NEUROSCOPE figure document (PDF)
Run: python make_figure_doc.py
Output: results/NEUROSCOPE_v52_figures.pdf
"""

import os
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY

OUT_PATH    = r"results\NEUROSCOPE_v52_figures.pdf"
KANG_DIR    = r"results\kang2018\paper_figures"
COVID_DIR   = r"results\covid_atlas\paper_figures"
LUPUS_DIR   = r"results\lupus\paper_figures"
SMILLIE_DIR = r"results\smillie2019\paper_figures"
IMG_WIDTH   = 6.2 * inch

# ── Styles ────────────────────────────────────────────────────────────────────

styles = getSampleStyleSheet()

title_style = ParagraphStyle("DocTitle",
    parent=styles["Title"], fontSize=18, spaceAfter=6, alignment=TA_CENTER)

subtitle_style = ParagraphStyle("DocSubtitle",
    parent=styles["Normal"], fontSize=11, spaceAfter=4, textColor=colors.HexColor("#555555"),
    alignment=TA_CENTER)

section_style = ParagraphStyle("Section",
    parent=styles["Heading1"], fontSize=13, spaceBefore=14, spaceAfter=4,
    textColor=colors.HexColor("#1a3a5c"), borderPad=2)

fig_title_style = ParagraphStyle("FigTitle",
    parent=styles["Normal"], fontSize=11, fontName="Helvetica-Bold",
    spaceBefore=10, spaceAfter=4)

legend_style = ParagraphStyle("Legend",
    parent=styles["Normal"], fontSize=9, leading=13, spaceAfter=6,
    alignment=TA_JUSTIFY)

stat_style = ParagraphStyle("Stat",
    parent=styles["Normal"], fontSize=9, fontName="Courier",
    leading=12, leftIndent=12, spaceAfter=2)

note_style = ParagraphStyle("Note",
    parent=styles["Normal"], fontSize=8, textColor=colors.HexColor("#666666"),
    leftIndent=12, spaceAfter=6, leading=11)

caption_style = ParagraphStyle("Caption",
    parent=styles["Normal"], fontSize=8, textColor=colors.HexColor("#444444"),
    alignment=TA_CENTER, spaceAfter=4)

# ── Helpers ───────────────────────────────────────────────────────────────────

def img(path, width=IMG_WIDTH):
    if os.path.exists(path):
        return Image(path, width=width, height=width * 0.42)
    return Paragraph(f"[Figure not yet generated: {os.path.basename(path)}]", note_style)

def stat_block(lines):
    return [Paragraph(l, stat_style) for l in lines]

def hr():
    return HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#cccccc"),
                      spaceAfter=6, spaceBefore=6)

def placeholder(text):
    return Paragraph(
        f"<i>[{text}]</i>",
        ParagraphStyle("ph", parent=legend_style, textColor=colors.HexColor("#888888"),
                       borderPad=6, backColor=colors.HexColor("#f5f5f5"))
    )

# ── Document assembly ─────────────────────────────────────────────────────────

story = []

# ── Title page ────────────────────────────────────────────────────────────────
story.append(Spacer(1, 1.2 * inch))
story.append(Paragraph("NEUROSCOPE v5.2", title_style))
story.append(Paragraph("Figure Document — Publication Draft", subtitle_style))
story.append(Spacer(1, 0.15 * inch))
story.append(hr())
story.append(Spacer(1, 0.1 * inch))
story.append(Paragraph(
    "Information-theoretic cell-cell communication analysis via KSG channel capacity. "
    "This document contains all publication figures with full legends, statistical details, "
    "and interpretation notes.",
    ParagraphStyle("intro", parent=legend_style, alignment=TA_CENTER, fontSize=10)))
story.append(Spacer(1, 0.2 * inch))

# Summary table — all 4 benchmarks
sum_data = [
    ["", "Kang 2018\n(Benchmark 1)", "COVID Atlas\n(Benchmark 2)",
     "Lupus PBMC\n(Benchmark 3)", "Smillie IBD Colon\n(Benchmark 4)"],
    ["Reference", "Kang et al. 2018\nGSE96583", "Stephenson et al. 2021",
     "Perez et al. 2022", "Smillie et al. 2019"],
    ["Tissue", "PBMC", "PBMC", "PBMC", "Colon (non-blood)"],
    ["Condition", "IFN-β stimulation", "COVID-19 severity",
     "SLE vs. healthy", "Crohn disease vs normal"],
    ["Donors / samples", "n=16 pseudobulk\n(8 donors × 2)", "n=124 donors",
     "n=162 donors\n(88 SLE + 49 ctrl)", "n=18 donors\n(9,385 cells)"],
    ["KSG-sig (FDR<0.05)", "5 / 9  (56%)", "49 / 127  (39%)",
     "Pending run", "17 / 34  (50%)"],
    ["Invisible to CellChat", "5 / 5  (100%)", "20 / 49  (41%)",
     "Pending run", "Pending LIANA"],
    ["Fisher OR vs CellChat", "0.000  (p=0.0079)", "0.060  (p=0.0008)",
     "Pending run", "Pending LIANA"],
    ["2SLS causal validated", "0 / 9", "TBD",
     "Expected: higher\n(n=162 donors)", "7 / 34  (20.6%)"],
]
sum_style = TableStyle([
    ("BACKGROUND",  (0,0), (-1,0),  colors.HexColor("#1a3a5c")),
    ("TEXTCOLOR",   (0,0), (-1,0),  colors.white),
    ("FONTNAME",    (0,0), (-1,0),  "Helvetica-Bold"),
    ("BACKGROUND",  (0,1), (0,-1),  colors.HexColor("#e8eef4")),
    ("FONTNAME",    (0,1), (0,-1),  "Helvetica-Bold"),
    ("ROWBACKGROUNDS", (1,1), (-1,-1), [colors.white, colors.HexColor("#f7f9fb")]),
    ("FONTSIZE",    (0,0), (-1,-1), 7.5),
    ("GRID",        (0,0), (-1,-1), 0.4, colors.HexColor("#cccccc")),
    ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
    ("TOPPADDING",  (0,0), (-1,-1), 4),
    ("BOTTOMPADDING",(0,0),(-1,-1), 4),
    ("LEFTPADDING", (0,0), (-1,-1), 5),
])
story.append(Table(sum_data,
                   colWidths=[1.55*inch, 1.35*inch, 1.35*inch, 1.35*inch, 1.35*inch],
                   style=sum_style))

story.append(Spacer(1, 0.15 * inch))
story.append(Paragraph(
    "Benchmark 3 (Lupus) is pending data download + MATLAB run. "
    "Benchmark 4 (Smillie IBD) is COMPLETE. "
    "Run prep_lupus.py, run_lupus in MATLAB, then re-generate this document.",
    note_style))
story.append(PageBreak())

# ── Figure 1 ─────────────────────────────────────────────────────────────────
story.append(Paragraph("Figure 1 — Communication Information Landscape", section_style))
story.append(hr())

# 1A — Kang
story.append(Paragraph("Figure 1A  (Kang 2018 PBMC, n = 16 pseudobulk samples)", fig_title_style))
story.append(img(os.path.join(KANG_DIR, "Fig1_information_landscape.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 1A.</b> Information landscape of cell-cell communication in IFN-β-stimulated PBMC "
    "(Kang et al. 2018, GSE96583). "
    "<b>Panel A:</b> Distribution of KSG channel capacity (bits) across all tested L-R pairs. "
    "Blue bars: FDR-significant pairs (FDR < 0.05 by permutation test). "
    "Grey bars: pairs with I ≈ 0. "
    "<b>Panel B:</b> Capacity vs. co-expression strength scatter (Spearman r shown). "
    "Discordant pairs (low strength, non-zero capacity) are invisible to geometric mean-based tools. "
    "<b>Panel C:</b> Enrichment classification — proportion of L-R pairs significant by capacity "
    "only, strength only, both, or neither.",
    legend_style))

story += stat_block([
    "KSG-significant pairs:  5 / 9  (FDR < 0.05, permutation, 199 shuffles)",
    "Top pairs by capacity:",
    "  CXCL10 → CXCR3    cap = 0.985  [95% CI: 0.076, 1.000]  FDR = 0.0044",
    "  CXCL11 → CXCR3    cap = 0.964  [95% CI: 0.043, 1.000]  FDR = 0.0005",
    "  IL1RN  → IL1R2    cap = 0.913  [95% CI: 0.018, 1.000]  FDR = 0.0281",
    "  CLU    → TYROBP   cap = 0.755  [95% CI: 0.027, 0.965]  FDR = 0.0117",
    "  CD47   → TYROBP   cap = 0.654  [95% CI: 0.000, 1.000]  FDR = 0.0281",
    "Mean capacity (all pairs):   0.8395 bits",
    "Median capacity:             1.0000 bits",
    "High-capacity (>0.5 bits):   45/51  (88%)",
])
story.append(Paragraph(
    "Note: wide 95% bootstrap CIs (500 donor resamples) reflect n = 16 pseudobulk observations. "
    "Relative ranking between pairs is valid; absolute capacity values should be interpreted "
    "with caution at this sample size.",
    note_style))

story.append(Spacer(1, 0.2 * inch))

# 1B — COVID
story.append(Paragraph("Figure 1B  (COVID Blood Atlas, n = 124 donors)", fig_title_style))
story.append(img(os.path.join(COVID_DIR, "Fig1_information_landscape.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 1B.</b> Information landscape of cell-cell communication in COVID-19 PBMC "
    "(Stephenson et al. 2021, 124 donors, 3 severity conditions). "
    "Same panel layout as Figure 1A. The larger donor cohort (n = 124) yields "
    "substantially narrower bootstrap confidence intervals.",
    legend_style))

story += stat_block([
    "KSG-significant pairs:  49 / 127  (FDR < 0.05)",
    "Top pairs by capacity (FDR < 0.05):",
    "  ALOX5   → CYSLTR1  cap = 0.722  [95% CI: 0.000, 1.000]  FDR = 0.0405",
    "  PTGES3  → PTGER4   cap = 0.692  [95% CI: 0.000, 1.000]  FDR = 0.0134",
    "  CD55    → ADGRE5   cap = 0.550  [95% CI: 0.000, 1.000]  FDR = 0.0019",
    "  ITGB1   → ADGRE5   cap = 0.526  [95% CI: 0.000, 1.000]  FDR = 0.0013",
    "  PPIA    → BSG      cap = 0.523  [95% CI: 0.000, 1.000]  FDR = 0.0019",
    "  HLA-E   → VSIR     cap = 0.510  [95% CI: 0.000, 1.000]  FDR = 0.0089",
    "  CD320   → JAML     cap = 0.439  [95% CI: 0.000, 0.793]  FDR = 0.0034",
    "Mean capacity (network):    0.2985 bits",
    "High-capacity (>0.5 bits):  145/687  (21.1%)",
    "Spectral communities:       2  (spectral gap = 0.8115)",
    "IV/2SLS causally validated: 1 / 127  (0.8%)",
    "Comm identity ARI:          0.0816 vs transcriptome  (ORTHOGONAL, 100% divergent)",
    "NMF programs:               3  (innate adhesion, eicosanoid+checkpoint, complement+adhesion)",
])
story.append(PageBreak())

# ── Figure 2 ─────────────────────────────────────────────────────────────────
story.append(Paragraph("Figure 2 — Channel Type Composition and Pathway Distribution", section_style))
story.append(hr())

story.append(Paragraph("Figure 2A  (Kang 2018)", fig_title_style))
story.append(img(os.path.join(KANG_DIR, "Fig2_channel_types.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 2A.</b> Channel type composition and pathway distribution for Kang 2018. "
    "<b>Panel A:</b> Count of L-R pairs in each channel type category. "
    "High-capacity (HC) channels have KSG capacity > 0.5 bits; "
    "Broadcast channels have high ligand Shannon entropy (many receiver cell types); "
    "Selective channels target one or few receiver populations. "
    "<b>Panel B:</b> Pathway breakdown of capacity-only enriched vs. both-enriched channels. "
    "Chemokine and Cytokine/JAK-STAT pathways dominate the IFN-β response.",
    legend_style))
story += stat_block([
    "Channel type breakdown (9 pairs total):",
    "  HighCapacity_Broadcast:   2  (22.2%)",
    "  HighCapacity_Selective:   2  (22.2%)",
    "  LowCapacity_Broadcast:    1  (11.1%)",
    "  LowCapacity_Selective:    4  (44.4%)",
    "Pathway distribution:  Chemokine 3, Cytokine/JAK-STAT 2, Other 4",
])

story.append(Spacer(1, 0.15 * inch))
story.append(Paragraph("Figure 2B  (COVID Atlas)", fig_title_style))
story.append(img(os.path.join(COVID_DIR, "Fig2_channel_types.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 2B.</b> Channel type composition for COVID Atlas (127 tested pairs). "
    "Eicosanoid (ALOX5→CYSLTR1, PTGES3→PTGER4), complement (CD55→ADGRE5), "
    "and immune checkpoint (HLA-E→VSIR) pathways appear among capacity-significant pairs.",
    legend_style))
story += stat_block([
    "Channel type breakdown (127 pairs):",
    "  HighCapacity_Broadcast:  34  (26.8%)",
    "  HighCapacity_Selective:  29  (22.8%)",
    "  LowCapacity_Broadcast:   25  (19.7%)",
    "  LowCapacity_Selective:   39  (30.7%)",
    "Pathway distribution:  Other 84, Cytokine/JAK-STAT 26, Chemokine 7,",
    "                        Semaphorin 4, Cell adhesion 2",
])
story.append(PageBreak())

# ── Figure 3 ─────────────────────────────────────────────────────────────────
story.append(Paragraph("Figure 3 — Synthetic Ground-Truth Benchmark (AUROC)", section_style))
story.append(hr())

fig3_path = r"results\kang2018\paper_figures\Fig3_benchmark_auroc.png"
story.append(img(fig3_path))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 3.</b> Synthetic ground-truth validation comparing KSG mutual information channel "
    "capacity to co-expression scoring (geometric mean, as used by CellChat/NicheNet). "
    "Simulation parameters: n = 200 synthetic cells per pair, "
    "20 causal pairs (Y = 0.6·X + noise), "
    "20 dead-receiver pairs (X and Y drawn independently, both expressed), "
    "20 null pairs (X zero-inflated lognormal, Y ≈ constant). "
    "Zero-inflation = 60% per cell (typical scRNA-seq sparsity). "
    "<b>Panel A:</b> ROC curves for recovering the 20 causal pairs. "
    "<b>Panel B:</b> AUROC summary.",
    legend_style))
story += stat_block([
    "Simulation settings:",
    "  n per pair     = 200 synthetic cells",
    "  Causal pairs   = 20   (Y = β·X + ε,  β = 0.6,  σ = 0.8)",
    "  Dead-receiver  = 20   (X, Y drawn independently — both expressed)",
    "  Null pairs     = 20   (Y ≈ 0, no signal)",
    "  Zero-inflation = 60%  per cell  (realistic scRNA-seq sparsity)",
    "  Random seed    = 42",
    "",
    "AUROC results:",
    "  KSG channel capacity   AUROC = 1.000",
    "  Co-expression score    AUROC = 0.506  (near-random)",
    "",
    "Dead-receiver pairs in top 20:",
    "  KSG channel capacity:    0 / 20",
    "  Co-expression score:  ~10 / 20  (false positives promoted by geometric mean)",
])
story.append(Paragraph(
    "Note: dead-receiver pairs are the critical failure mode for co-expression tools — "
    "both ligand and receptor are expressed and produce a high geometric mean, but no "
    "information is transmitted from sender to receiver. The KSG estimator assigns near-zero "
    "capacity to these pairs because mutual information requires statistical dependence, "
    "not merely co-expression.",
    note_style))
story.append(Paragraph(
    "Figure 3 is generated by simulate_communication_ground_truth.m with 'plot', false "
    "and then rendered in generate_paper_figures.m. Re-run NEUROSCOPE to regenerate. "
    "The AUROC values above are from the seeded simulation (seed=42) and are reproducible.",
    note_style))
story.append(PageBreak())

# ── Figure 4 ─────────────────────────────────────────────────────────────────
story.append(Paragraph("Figure 4 — CellChat Comparison: Overlapping vs. Exclusive Pairs", section_style))
story.append(hr())

fig4_path = r"results\kang2018\paper_figures\Fig4_cellchat_comparison.png"
story.append(img(fig4_path))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 4.</b> Comparison of NEUROSCOPE (KSG capacity, FDR < 0.05) and CellChat "
    "(trimean permutation test, FDR < 0.05) on the Kang 2018 benchmark. "
    "<b>Panel A:</b> Contingency — pairs found by both methods, NEUROSCOPE only, "
    "or CellChat only. "
    "<b>Panel B:</b> Among NEUROSCOPE-only pairs, breakdown of why CellChat missed them: "
    "zero-floor (trimean = 0 for all cluster combinations) vs. other reasons.",
    legend_style))

story.append(Spacer(1, 0.1 * inch))
story.append(Paragraph("Kang 2018  (n = 16 pseudobulk, 9 pairs tested)", fig_title_style))

kang_table_data = [
    ["", "CellChat significant", "CellChat not significant", "Total"],
    ["NEUROSCOPE significant", "0", "5", "5"],
    ["NEUROSCOPE not significant", "4", "0", "4"],
    ["Total", "4", "5", "9"],
]
kang_ts = TableStyle([
    ("BACKGROUND",  (0,0), (-1,0),  colors.HexColor("#1a3a5c")),
    ("TEXTCOLOR",   (0,0), (-1,0),  colors.white),
    ("FONTNAME",    (0,0), (-1,0),  "Helvetica-Bold"),
    ("BACKGROUND",  (0,1), (0,-1),  colors.HexColor("#e8eef4")),
    ("FONTNAME",    (0,1), (0,-1),  "Helvetica-Bold"),
    ("FONTSIZE",    (0,0), (-1,-1), 9),
    ("GRID",        (0,0), (-1,-1), 0.4, colors.HexColor("#aaaaaa")),
    ("ALIGN",       (1,0), (-1,-1), "CENTER"),
    ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
    ("TOPPADDING",  (0,0), (-1,-1), 5),
    ("BOTTOMPADDING",(0,0),(-1,-1),5),
])
story.append(Table(kang_table_data, colWidths=[2.0*inch, 1.8*inch, 1.8*inch, 0.9*inch],
                   style=kang_ts))
story += stat_block([
    "Fisher exact test:  OR = 0.000,  p = 0.0079",
    "NEUROSCOPE-only pairs invisible to CellChat:  5 / 5  (100%)",
    "  Zero-floor (trimean = 0):  4 / 5",
    "  Other:                     1 / 5  (CD47→TYROBP, trimean > 0 but NS)",
    "",
    "NEUROSCOPE-only top pairs:",
    "  CXCL10 → CXCR3   cap=0.985  ZeroFloor=YES  (IFN-inducible chemokine)",
    "  CXCL11 → CXCR3   cap=0.964  ZeroFloor=YES  (IFN-inducible chemokine)",
    "  IL1RN  → IL1R2   cap=0.913  ZeroFloor=YES  (cytokine regulator)",
    "  CLU    → TYROBP  cap=0.755  ZeroFloor=YES  (complement regulator)",
    "  CD47   → TYROBP  cap=0.654  ZeroFloor=NO   (don't-eat-me signal)",
])

story.append(Spacer(1, 0.2 * inch))
story.append(Paragraph("COVID Atlas  (n = 124 donors, 127 pairs tested)", fig_title_style))

covid_table_data = [
    ["", "CellChat significant", "CellChat not significant", "Total"],
    ["NEUROSCOPE significant", "29", "20", "49"],
    ["NEUROSCOPE not significant", "24", "0", "24"],
    ["Total", "53", "20", "73"],
]
story.append(Table(covid_table_data, colWidths=[2.0*inch, 1.8*inch, 1.8*inch, 0.9*inch],
                   style=kang_ts))
story += stat_block([
    "Fisher exact test:  OR = 0.060,  p = 0.0008",
    "NEUROSCOPE-only pairs invisible to CellChat:  20 / 49  (41%)",
    "  Zero-floor (trimean = 0):  18 / 20  (90%)",
    "  Other:                      2 / 20",
    "",
    "Top NEUROSCOPE-only pairs (invisible to CellChat):",
    "  PTGES3  → PTGER4  cap=0.692  ZeroFloor=YES  (eicosanoid signalling)",
    "  B2M     → IL6ST   cap=0.462  ZeroFloor=YES  (IL-6 trans-signalling)",
    "  CD52    → SIGLEC10 cap=0.372 ZeroFloor=YES  (immune checkpoint)",
    "  CCL5    → CCR5    cap=0.360  ZeroFloor=YES  (monocyte chemotaxis)",
    "  APP     → SORL1   cap=0.337  ZeroFloor=YES  (APP processing)",
    "  SEMA4D  → CD72    cap=0.332  ZeroFloor=YES  (semaphorin signalling)",
    "  TGFB1   → TGFBR2  cap=0.173  ZeroFloor=YES  (TGF-β pathway)",
    "  IFNG    → IFNGR1  cap=0.065  ZeroFloor=YES  (IFN-γ response)",
])
story.append(Paragraph(
    "Mechanistic explanation — CellChat zero-floor: "
    "trimean(x) = 0.25·Q1 + 0.5·Q2 + 0.25·Q3 (Jin et al. 2021, eq. 1). "
    "When Q1 = Q2 = 0 (gene expressed in fewer than 50% of cells), trimean = 0 exactly. "
    "This is a structural property of the formula, not a data artifact. "
    "NEUROSCOPE aggregates to per-donor cluster means before KSG estimation, "
    "retaining inter-donor variation even at 15% per-cell sparsity.",
    note_style))
story.append(PageBreak())

# ── Figure 5 ─────────────────────────────────────────────────────────────────
story.append(Paragraph("Figure 5 — Latent Communication Programs (NMF)", section_style))
story.append(hr())

story.append(Paragraph("Figure 5A  (Kang 2018, 7 programs)", fig_title_style))
story.append(img(os.path.join(KANG_DIR, "Fig5_nmf_programs.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 5A.</b> Non-negative matrix factorization of the "
    "[cluster-pair × L-R] communication tensor, Kang 2018. "
    "<b>Panel A:</b> Heatmap of program loadings (top 20 L-R pairs). "
    "Rows = programs, columns = L-R pairs. "
    "<b>Panel B:</b> Program purity score (fraction of weight concentrated in top pairs). "
    "Rank = 7 selected by elbow criterion (<10% reconstruction improvement).",
    legend_style))
story += stat_block([
    "NMF rank:  7 programs",
    "Program 1: CXCL11→CXCR3, CXCL10→CXCR3, CD47→TYROBP   (IFN-inducible chemokines)",
    "Program 2: IL1RN→IL1R2,  IL1B→IL1R2,   CLU→TYROBP    (cytokine regulation)",
    "Program 3: CD47→TYROBP,  CLU→TYROBP,   PF4→CXCR3     (myeloid–NK crosstalk)",
    "Program 4: CXCL10→CXCR3, CXCL11→CXCR3, CLU→TYROBP",
    "Program 5: PPIA→BSG,     CLU→TYROBP,   CD47→TYROBP",
    "Programs 6-7: Minor programs overlapping with Programs 1-3",
])

story.append(Spacer(1, 0.2 * inch))
story.append(Paragraph("Figure 5B  (COVID Atlas, 3 programs)", fig_title_style))
story.append(img(os.path.join(COVID_DIR, "Fig5_nmf_programs.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 5B.</b> NMF communication programs for COVID Atlas (n = 124 donors, 127 pairs). "
    "Rank = 3 (elbow criterion). Three biologically distinct programmes emerge: "
    "innate immune adhesion, eicosanoid/checkpoint signalling, and complement/adhesion.",
    legend_style))
story += stat_block([
    "NMF rank:  3 programs",
    "Program 1: ICAM3→ITGB2,  ANXA1→FPR1,   APP→TYROBP   (innate immune adhesion)",
    "Program 2: PPIA→BSG,     PTGES3→PTGER4, HLA-E→VSIR   (eicosanoid + checkpoint)",
    "Program 3: CD55→ADGRE5,  ITGB1→ADGRE5, CD58→CD2     (complement + adhesion)",
    "",
    "Spectral communities (Laplacian decomposition):  2",
    "  Spectral gap = 0.8115  (high modularity)",
    "Downstream response genes (17 receptors with transcriptional signature):",
    "  e.g. ADGRE5, ITGB2, CD2, BSG, PTGER4",
])

story.append(Spacer(1, 0.2 * inch))
story.append(Paragraph("Figure 5C  (COVID Atlas — Communication identity vs transcriptome)", fig_title_style))
story.append(img(os.path.join(COVID_DIR, "Fig6_communication_identity.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 5C.</b> Communication state identity vs. transcriptomic cluster identity for "
    "COVID Atlas (n = 124 donors, 7 cell types). ARI = 0.082 — ORTHOGONAL. "
    "100% of 24,284 cells are classified as divergent (CDS > 0.5). "
    "Replicates the same finding seen in Smillie IBD (ARI = 0.031, 100% divergent), "
    "confirming across tissue type and disease that communication identity cannot be "
    "derived from transcriptomic clustering alone.",
    legend_style))
story.append(PageBreak())

# ── Figure 6 — Benchmarks 3 & 4 placeholders ─────────────────────────────────
story.append(Paragraph("Figure 6 — Benchmark 3: Lupus PBMC (Perez et al. 2022, n = 162 donors)", section_style))
story.append(hr())
story.append(placeholder(
    "Figure 6A pending: run prep_lupus.py then run_lupus.m in MATLAB. "
    "Expected: tight bootstrap CIs from n=162 donors; stronger 2SLS causal power. "
    "Predicted top L-R axes: IFN signature (CXCL10→CXCR3, IFNG→IFNGR1), "
    "complement (C1Q→LRP1), and Fc-receptor pathways. "
    "Benchmarks 1-3 complete (Kang, COVID, Smillie) — orthogonal comm identity replicated in all three."))
story.append(Spacer(1, 0.1 * inch))
story += stat_block([
    "Dataset: Perez et al. 2022 (GSE174188 / CellxGene)",
    "162 donors: 88 SLE + 49 healthy + misc",
    "~1.2M cells, 9 PBMC cell types",
    "Config: n_permutations=200, fdr_threshold=0.05, run_ml=false",
    "Expected advantage: n=162 gives BCa CI width ~0.3 vs ~1.0 at n=16",
    "Expected advantage: 2SLS first-stage F>10 feasible at n=162 donors",
    "Run command: matlab -batch \"cd c:/arjun-project/neuroscope; run_lupus\"",
    "Status: URL fixed in prep_lupus.py — run prep_lupus.py to begin data download (10.7 GB)",
])

story.append(Spacer(1, 0.2 * inch))
story.append(Paragraph("Figure 6B — Benchmark 4: IBD Colon (Smillie et al. 2019, non-blood tissue) — COMPLETE", fig_title_style))
story.append(img(os.path.join(SMILLIE_DIR, "Fig1_information_landscape.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 6B.</b> Information landscape for IBD (Crohn disease) colon immune cells "
    "(Smillie et al. 2019, CellxGene). 9,385 cells, 18 donors, 17 immune cell types. "
    "First non-blood-tissue benchmark — demonstrates NEUROSCOPE generalizability beyond PBMC. "
    "Strong cytokine (IL1B, IFNG), T-cell (LCK→CD8A), and adhesion (ICAM1→ITGB2) axes identified. "
    "Communication identity is ORTHOGONAL to transcriptome (ARI = 0.031, 100% divergent cells).",
    legend_style))
story.append(Spacer(1, 0.05 * inch))
story.append(img(os.path.join(SMILLIE_DIR, "Fig2_channel_types.png")))
story.append(Spacer(1, 0.05 * inch))
story += stat_block([
    "Dataset: Smillie et al. 2019 colon immune (CellxGene)",
    "  Condition: Crohn disease vs normal   |   Donors: n=18   |   Cells: 9,385",
    "  Cell types: 17 immune types          |   Genes (HVG): 3,000",
    "",
    "Network results:",
    "  Total L-R interactions (cap > 0):  416",
    "  Pairs tested (differential):       34",
    "  KSG-significant (FDR<0.05):        17 / 34  (50%)",
    "  Mean channel capacity:             0.3511 bits",
    "  High-capacity (>0.5 bits):         107 / 416  (25.7%)",
    "",
    "Top channels by capacity:",
    "  LCK    -> CD8A    cap=0.6313  FDR=0.0007  T-cell co-receptor signalling",
    "  IL1B   -> IL1R2   cap=0.6242  FDR=0.0312  Cytokine/JAK-STAT (IBD inflam.)",
    "  VEGFA  -> NRP1    cap=0.5828  FDR=0.0005  Angiogenesis/EGF signalling",
    "  IFNG   -> IFNGR1  cap=0.4450  FDR=0.0001  IFN-gamma (IBD mucosal immunity)",
    "  THBS1  -> CD36    cap=0.4370  FDR=0.0015  Thrombospondin adhesion",
    "  KLRB1  -> CLEC2D  cap=0.3823  FDR=0.0000  NK-T co-stimulation",
    "  ICAM1  -> ITGB2   cap=0.3551  FDR=0.0000  Leukocyte adhesion (IBD homing)",
    "",
    "NMF programs (rank=3):",
    "  Program 1: ICAM1_ITGB2, VEGFA_NRP1, APLP2_PIGR   (adhesion/angiogenesis)",
    "  Program 2: KLRB1_CLEC2D, TNF_TNFRSF1B, IFNG_IFNGR1  (immune activation)",
    "  Program 3: LCK_CD8A, APLP2_PIGR, IFNG_IFNGR1        (T-cell signalling)",
    "",
    "Spectral communities: 2  (spectral gap = 0.6730)",
    "2SLS causal validated:  7 / 34  (20.6%)  [donor_id instrument, n=18]",
    "Communication identity: ARI=0.031 vs transcriptome  -> ORTHOGONAL (100% divergent)",
    "Synthetic AUROC: KSG=1.000  vs  co-expression=0.506  (seed=42)",
    "Runtime: 5.8 min",
])

story.append(Spacer(1, 0.15 * inch))
story.append(Paragraph("Figure 6B-ii  (Smillie IBD — NMF programs)", fig_title_style))
story.append(img(os.path.join(SMILLIE_DIR, "Fig5_nmf_programs.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 6B-ii.</b> NMF communication programs for Smillie IBD (rank=3). "
    "Three co-regulated modules: (1) adhesion/angiogenesis, (2) immune activation (TNF, IFNG), "
    "(3) T-cell signalling (LCK→CD8A). Cross-module co-activation of APLP2→PIGR suggests "
    "shared transcytosis programme across immune subsets.",
    legend_style))

story.append(Spacer(1, 0.15 * inch))
story.append(Paragraph("Figure 6B-iii  (Smillie IBD — Communication identity vs transcriptome)", fig_title_style))
story.append(img(os.path.join(SMILLIE_DIR, "Fig6_communication_identity.png")))
story.append(Spacer(1, 0.05 * inch))
story.append(Paragraph(
    "<b>Figure 6B-iii.</b> Communication state identity vs. transcriptomic cluster identity "
    "(Smillie IBD). ARI = 0.031 — ORTHOGONAL. 100% of cells are classified as divergent "
    "(CDS > 0.5). This means a cell's communication behavior CANNOT be predicted from its "
    "transcriptomic cluster label — a robust finding replicated across all three completed benchmarks.",
    legend_style))
story.append(PageBreak())

# ── Figure 7 — NicheNet comparison ───────────────────────────────────────────
story.append(Paragraph("Figure 7 — NicheNet Comparison (Prior-Knowledge vs. KSG)", section_style))
story.append(hr())
story.append(Paragraph(
    "<b>Figure 7.</b> Comparison of NEUROSCOPE (KSG channel capacity, data-driven) "
    "against NicheNet (prior-knowledge ligand activity scoring, Browaeys et al. 2020 "
    "Nat Methods). NicheNet ranks ligands by their ability to induce downstream target "
    "gene sets in receiver cells via curated signaling networks — an orthogonal approach "
    "to KSG mutual information.",
    legend_style))
story.append(placeholder(
    "Figure 7 pending: run Rscript compare_nichenet.R. "
    "Requires R 4.x, nichenetr, tidyverse. NicheNet prior networks downloaded from Zenodo. "
    "Expected output: results/kang2018/comparison_nichenet_vs_neuroscope.csv"))
story.append(Spacer(1, 0.1 * inch))
story += stat_block([
    "Script: compare_nichenet.R",
    "Install: install.packages('BiocManager'); BiocManager::install('nichenetr')",
    "Run: Rscript compare_nichenet.R",
    "",
    "NicheNet uses: Ligand -> TF -> target gene enrichment (AUPR over gene sets)",
    "NEUROSCOPE uses: KSG MI channel capacity (statistical dependence, no prior network)",
    "",
    "Interpretation of overlap:",
    "  Both agree:          high-confidence (two orthogonal methods converge)",
    "  NEUROSCOPE-only:     info-theoretically active, absent from prior PPI networks",
    "                       (novel interactions OR recently characterized)",
    "  NicheNet-only:       prior-knowledge supported but statistically weak at n-donors",
    "",
    "Expected: partial overlap (~30-50%) — methods measure different things",
    "Key test: do IFN-inducible chemokines (CXCL10, CXCL11) rank top in both?",
])

story.append(Spacer(1, 0.2 * inch))
story.append(Paragraph("Figure 7 Supplemental - LIANA Multi-Method Comparison (COMPLETED)", fig_title_style))
story.append(Paragraph(
    "LIANA (Dimitrov et al. 2022 Nat Cell Biol) wraps CellPhoneDB, CellChat, NicheNet, "
    "NATMI, and Connectome in a unified framework. "
    "Method used: CellPhoneDB (mean(L)*mean(R), permutation, expr_prop=0.1). "
    "CellPhoneDB does not have the CellChat trimean zero-floor, making this a "
    "conservative comparison for NEUROSCOPE.",
    legend_style))
story.append(Spacer(1, 0.1 * inch))

liana_data = [
    ["Dataset", "NEUROSCOPE sig", "In LIANA", "NEUROSCOPE-only", "Exclusive %"],
    ["Kang 2018", "5", "2", "3", "60%"],
    ["COVID Atlas", "49", "17", "32", "65%"],
    ["Smillie IBD (Crohn)", "17", "5", "12", "71%"],
]
liana_ts = TableStyle([
    ("BACKGROUND",  (0,0), (-1,0),  colors.HexColor("#1a3a5c")),
    ("TEXTCOLOR",   (0,0), (-1,0),  colors.white),
    ("FONTNAME",    (0,0), (-1,0),  "Helvetica-Bold"),
    ("BACKGROUND",  (0,1), (0,-1),  colors.HexColor("#e8eef4")),
    ("FONTNAME",    (0,1), (0,-1),  "Helvetica-Bold"),
    ("FONTSIZE",    (0,0), (-1,-1), 9),
    ("GRID",        (0,0), (-1,-1), 0.4, colors.HexColor("#aaaaaa")),
    ("ALIGN",       (1,0), (-1,-1), "CENTER"),
    ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
    ("TOPPADDING",  (0,0), (-1,-1), 5),
    ("BOTTOMPADDING",(0,0),(-1,-1), 5),
])
story.append(Table(liana_data,
                   colWidths=[1.5*inch, 1.3*inch, 1.3*inch, 1.3*inch, 1.3*inch],
                   style=liana_ts))
story.append(Spacer(1, 0.1 * inch))
story += stat_block([
    "NEUROSCOPE-only pairs (not in LIANA) Kang 2018:",
    "  IL1RN  -> IL1R2   cap=0.913  (cytokine regulator, FDA-approved biologic target)",
    "  CLU    -> TYROBP  cap=0.755  (complement regulator)",
    "  CD47   -> TYROBP  cap=0.654  (phagocytosis checkpoint, cancer immunotherapy target)",
    "",
    "NEUROSCOPE-only pairs (not in LIANA) COVID Atlas (selected):",
    "  ALOX5  -> CYSLTR1  cap=0.722  (eicosanoid/leukotriene signalling)",
    "  PTGES3 -> PTGER4   cap=0.692  (prostaglandin immune modulation)",
    "  PPIA   -> BSG      cap=0.523  (cyclophilin A / Basigin - COVID receptor)",
    "  HLA-E  -> VSIR     cap=0.510  (immune checkpoint, NK inhibitory axis)",
    "  CCL5   -> CCR5     cap=0.360  (monocyte chemotaxis, elevated in severe COVID)",
    "  TGFB1  -> TGFBR2   cap=0.173  (T cell exhaustion / COVID immunopathology)",
    "",
    "NEUROSCOPE-only pairs (not in LIANA) Smillie IBD (12/17 = 71%):",
    "  LCK    -> CD8A    cap=0.631  (T-cell co-receptor — top channel invisible to LIANA!)",
    "  IL1B   -> IL1R2   cap=0.624  (IL-1 cytokine — key IBD inflammatory axis)",
    "  IFNG   -> IFNGR1  cap=0.445  (IFN-gamma — mucosal immunity, invisible to LIANA)",
    "  KLRB1  -> CLEC2D  cap=0.382  (NK co-stimulation)",
    "  ICAM1  -> ITGB2   cap=0.355  (leukocyte adhesion — IBD homing)",
    "  ANXA1  -> FPR3    cap=0.334  (annexin — inflammatory resolution)",
    "  CLU    -> TYROBP  cap=0.340  (complement — also NEUROSCOPE-only in Kang!)",
    "  APLP2  -> PIGR    cap=0.330  (transcytosis circuit)",
    "",
    "Pattern across 3 benchmarks: 60-71% of NEUROSCOPE-significant pairs not found by LIANA,",
    "  consistent regardless of tissue (PBMC or colon) or disease (viral, autoimmune, IBD).",
    "  KSG MI captures inter-donor statistical dependence; LIANA methods detect mean expression.",
])
story.append(PageBreak())

# ── Methods summary ───────────────────────────────────────────────────────────
story.append(Paragraph("Methods Summary", section_style))
story.append(hr())

methods = [
    ("<b>Channel capacity estimation</b>",
     "KSG k-NN mutual information estimator "
     "(Kraskov, Stögbauer & Grassberger 2004, Phys Rev E 69:066138). "
     "I(X;Y) = ψ(k) + ψ(n) − ⟨ψ(nₓ+1) + ψ(nᵧ+1)⟩. "
     "Capacity = I(ligand ; response) / H(ligand) ∈ [0, 1]. "
     "Input: per-donor cluster means (pseudobulk). "
     "Adaptive k (default k=3, floored at ⌊n/5⌋). "
     "Spearman copula fallback when n < 10."),
    ("<b>Permutation test</b>",
     "199 random shuffles of the receptor vector Y. "
     "Empirical p = fraction of null MI ≥ observed MI. "
     "BH-FDR correction across all tested L-R pairs."),
    ("<b>Bootstrap confidence intervals</b>",
     "500 resamples (with replacement) of the n paired donor observations. "
     "95% CI = [2.5th, 97.5th] percentile of bootstrap capacity distribution. "
     "Adaptive k applied consistently within each resample."),
    ("<b>Normalization</b>",
     "Library-size median-ratio (scran-style): size factor = cell library size / median. "
     "Followed by log1p transform."),
    ("<b>HVG selection</b>",
     "Degree-2 polynomial regression log(variance) ~ log(mean). "
     "Top 3000 genes by residual variance (Pearson residuals)."),
    ("<b>Differential channel testing</b>",
     "Two-sample Wilcoxon rank-sum on per-cluster-pair strengths (strength test) "
     "and capacities (capacity test). BH-FDR correction."),
    ("<b>NMF programs</b>",
     "Non-negative matrix factorization of the [cluster-pair × L-R] communication "
     "tensor. Rank selected by elbow (<10% reconstruction improvement)."),
    ("<b>Causal validation (IV/2SLS)</b>",
     "Two-stage least squares with donor_id as natural instrumental variable "
     "(analogous to Mendelian Randomization in GWAS). "
     "Weak instrument filter: first-stage F > 10. "
     "Second-stage p < 0.05 required."),
    ("<b>CellChat comparison</b>",
     "CellChat v1 (Jin et al. 2021 Nat Commun). Trimean statistic: "
     "0.25·Q1 + 0.5·Q2 + 0.25·Q3. "
     "1000 permutations. FDR < 0.05 threshold. "
     "Zero-floor: trimean = 0 when Q1 = Q2 = 0 (>50% sparsity)."),
    ("<b>CellPhoneDB comparison</b>",
     "CellPhoneDB v3 (Efremova et al. 2020 Nat Protoc). "
     "Mean(L in cluster A) × mean(R in cluster B), min 10% expression, "
     "1000 permutations. FDR < 0.05."),
    ("<b>NicheNet comparison</b>",
     "NicheNet v2 (Browaeys et al. 2020 Nat Methods; nichenetr R package). "
     "Ligand activity = AUPR over receiver DE gene set enrichment via curated "
     "TF-target, PPI, and signaling networks from Zenodo (record 7074291). "
     "Run: Rscript compare_nichenet.R. "
     "Top-20 NicheNet ligands compared to NEUROSCOPE FDR-significant pairs."),
    ("<b>LIANA multi-method comparison</b>",
     "LIANA (Dimitrov et al. 2022 Nat Cell Biol; liana-py ≥1.0). "
     "Wraps CellPhoneDB, CellChat, NicheNet, NATMI, and Connectome. "
     "rank_aggregate: consensus Borda-count rank across all methods. "
     "expr_prop = 0.1 (min 10% expressing cells). "
     "Run: python compare_liana.py. "
     "Install: pip install liana anndata scanpy."),
    ("<b>Software</b>",
     "MATLAB R2021b+. Statistics and Machine Learning Toolbox. "
     "Python 3.x for data preparation (anndata, scanpy, scipy, numpy, pandas, liana). "
     "R 4.x for NicheNet comparison (nichenetr, tidyverse, Seurat)."),
]

for title, body in methods:
    story.append(KeepTogether([
        Paragraph(title, ParagraphStyle("mth", parent=legend_style, fontName="Helvetica",
                                         spaceBefore=8, spaceAfter=2, fontSize=9)),
        Paragraph(body, legend_style),
    ]))

story.append(Spacer(1, 0.2 * inch))
story.append(hr())
story.append(Paragraph(
    "NEUROSCOPE v5.2  ·  Generated 2026-05-19  ·  "
    "Benchmarks 1-2 complete. Benchmarks 3-4 and comparisons (NicheNet/LIANA) pending data prep. "
    "Figs 3-4 require MATLAB re-run; Figs 6-7 require benchmark data prep + comparison scripts.",
    caption_style))

# ── Build ─────────────────────────────────────────────────────────────────────

os.makedirs("results", exist_ok=True)
doc = SimpleDocTemplate(
    OUT_PATH,
    pagesize=letter,
    leftMargin=0.85*inch, rightMargin=0.85*inch,
    topMargin=0.75*inch,  bottomMargin=0.75*inch,
)
doc.build(story)
print(f"PDF written: {OUT_PATH}")
