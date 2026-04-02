"""
generate_presentation_pdf.py
============================
Generates a polished PDF presentation narrative covering the SmartChem ML work.
Covers: model selection, dataset, training, performance, and research paper topic.

Run from project root:
    python misc/generate_presentation_pdf.py
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTS_DIR = os.path.join(ROOT, "evaluation", "plots")
OUTPUT_PATH = os.path.join(ROOT, "SmartChem_ML_Presentation_Narrative.pdf")

from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm, mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether
)
from reportlab.platypus.flowables import HRFlowable
from reportlab.pdfgen import canvas
from PIL import Image as PILImage

# ── Color Palette ─────────────────────────────────────────────────────────────
DEEP_NAVY   = colors.HexColor("#0D1B2A")
MID_NAVY    = colors.HexColor("#1B2A3B")
ACCENT_BLUE = colors.HexColor("#1E88E5")
ACCENT_TEAL = colors.HexColor("#00BCD4")
LIGHT_GRAY  = colors.HexColor("#F5F7FA")
MED_GRAY    = colors.HexColor("#8899AA")
WHITE       = colors.white
HIGHLIGHT   = colors.HexColor("#E3F2FD")
WARN_AMBER  = colors.HexColor("#FF8F00")
GREEN_OK    = colors.HexColor("#2E7D32")

PAGE_W, PAGE_H = A4
MARGIN = 2.0 * cm


# ── Custom Canvas with header/footer ─────────────────────────────────────────
class BrandedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self._draw_header_footer(num_pages)
            super().showPage()
        super().save()

    def _draw_header_footer(self, page_count):
        page_num = self._pageNumber

        # Skip header/footer on the cover page (page 1)
        if page_num == 1:
            return

        self.saveState()

        # ── Top accent bar ────────────────────────────────────────────────────
        self.setFillColor(ACCENT_BLUE)
        self.rect(0, PAGE_H - 8 * mm, PAGE_W, 8 * mm, stroke=0, fill=1)
        self.setFillColor(WHITE)
        self.setFont("Helvetica-Bold", 8)
        self.drawString(MARGIN, PAGE_H - 5.5 * mm, "SmartChem | ML Presentation Narrative")
        self.setFont("Helvetica", 8)
        self.drawRightString(PAGE_W - MARGIN, PAGE_H - 5.5 * mm, "Confidential — Internal Use Only")

        # ── Bottom bar ────────────────────────────────────────────────────────
        self.setFillColor(DEEP_NAVY)
        self.rect(0, 0, PAGE_W, 10 * mm, stroke=0, fill=1)
        self.setFillColor(MED_GRAY)
        self.setFont("Helvetica", 7.5)
        self.drawString(MARGIN, 3.5 * mm, "Final Year Project  ·  SmartChem: AI-Powered Drug Discovery Platform")
        self.setFillColor(WHITE)
        self.drawRightString(PAGE_W - MARGIN, 3.5 * mm, f"Page {page_num} / {page_count}")

        self.restoreState()


# ── Style Sheet ───────────────────────────────────────────────────────────────
def build_styles():
    base = getSampleStyleSheet()

    styles = {}

    styles["cover_title"] = ParagraphStyle(
        "cover_title",
        fontName="Helvetica-Bold",
        fontSize=30,
        textColor=WHITE,
        alignment=TA_CENTER,
        spaceAfter=6,
        leading=36,
    )
    styles["cover_sub"] = ParagraphStyle(
        "cover_sub",
        fontName="Helvetica",
        fontSize=13,
        textColor=colors.HexColor("#B0C4D8"),
        alignment=TA_CENTER,
        spaceAfter=4,
        leading=18,
    )
    styles["cover_tag"] = ParagraphStyle(
        "cover_tag",
        fontName="Helvetica-Oblique",
        fontSize=10,
        textColor=ACCENT_TEAL,
        alignment=TA_CENTER,
        spaceAfter=2,
    )
    styles["section_header"] = ParagraphStyle(
        "section_header",
        fontName="Helvetica-Bold",
        fontSize=16,
        textColor=WHITE,
        spaceAfter=2,
        spaceBefore=0,
        leading=20,
    )
    styles["section_num"] = ParagraphStyle(
        "section_num",
        fontName="Helvetica-Bold",
        fontSize=9,
        textColor=ACCENT_TEAL,
        spaceAfter=1,
        spaceBefore=0,
    )
    styles["body"] = ParagraphStyle(
        "body",
        fontName="Helvetica",
        fontSize=10.5,
        textColor=colors.HexColor("#1A1A2E"),
        alignment=TA_JUSTIFY,
        spaceAfter=8,
        leading=16,
    )
    styles["body_lead"] = ParagraphStyle(
        "body_lead",
        fontName="Helvetica-Bold",
        fontSize=11,
        textColor=DEEP_NAVY,
        alignment=TA_JUSTIFY,
        spaceAfter=6,
        leading=16,
    )
    styles["bullet"] = ParagraphStyle(
        "bullet",
        fontName="Helvetica",
        fontSize=10,
        textColor=colors.HexColor("#1A1A2E"),
        leftIndent=18,
        spaceAfter=5,
        leading=15,
        bulletIndent=6,
        bulletFontName="Helvetica-Bold",
        bulletFontSize=13,
        bulletColor=ACCENT_BLUE,
    )
    styles["sub_bullet"] = ParagraphStyle(
        "sub_bullet",
        fontName="Helvetica",
        fontSize=9.5,
        textColor=colors.HexColor("#333355"),
        leftIndent=34,
        spaceAfter=3,
        leading=14,
    )
    styles["callout_label"] = ParagraphStyle(
        "callout_label",
        fontName="Helvetica-Bold",
        fontSize=9,
        textColor=ACCENT_BLUE,
        spaceAfter=2,
    )
    styles["callout_body"] = ParagraphStyle(
        "callout_body",
        fontName="Helvetica",
        fontSize=9.5,
        textColor=colors.HexColor("#1A1A2E"),
        alignment=TA_JUSTIFY,
        spaceAfter=4,
        leading=14,
    )
    styles["metric_label"] = ParagraphStyle(
        "metric_label",
        fontName="Helvetica-Bold",
        fontSize=10,
        textColor=DEEP_NAVY,
    )
    styles["metric_val"] = ParagraphStyle(
        "metric_val",
        fontName="Helvetica-Bold",
        fontSize=18,
        textColor=ACCENT_BLUE,
        alignment=TA_CENTER,
    )
    styles["caption"] = ParagraphStyle(
        "caption",
        fontName="Helvetica-Oblique",
        fontSize=8.5,
        textColor=MED_GRAY,
        alignment=TA_CENTER,
        spaceAfter=10,
    )
    styles["quote"] = ParagraphStyle(
        "quote",
        fontName="Helvetica-Oblique",
        fontSize=10.5,
        textColor=colors.HexColor("#2C3E6B"),
        leftIndent=18,
        rightIndent=18,
        spaceAfter=8,
        leading=16,
        borderPadding=(8, 8, 8, 12),
    )

    return styles


# ── Helper Flowables ──────────────────────────────────────────────────────────
def section_banner(num_str, title_text, styles, color=ACCENT_BLUE):
    """Dark banner with section number + title."""
    content_w = PAGE_W - 2 * MARGIN

    num_para   = Paragraph(num_str, styles["section_num"])
    title_para = Paragraph(title_text, styles["section_header"])

    inner = Table(
        [[num_para], [title_para]],
        colWidths=[content_w - 1.6 * cm],
        hAlign="LEFT",
    )
    inner.setStyle(TableStyle([
        ("LEFTPADDING",  (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING",   (0, 0), (0, 0),   10),
        ("BOTTOMPADDING",(0, -1),(-1, -1), 10),
    ]))

    outer = Table(
        [[inner]],
        colWidths=[content_w],
        hAlign="LEFT",
    )
    outer.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, -1), DEEP_NAVY),
        ("LEFTPADDING",  (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING",   (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 0),
        ("ROUNDEDCORNERS", (0, 0), (-1, -1), [6, 6, 6, 6]),
    ]))
    return [outer, Spacer(1, 0.35 * cm)]


def callout_box(label, body_text, styles, bg=HIGHLIGHT, border=ACCENT_BLUE):
    """Highlighted info / callout box."""
    content_w = PAGE_W - 2 * MARGIN
    lbl  = Paragraph(f"▌  {label}", styles["callout_label"])
    body = Paragraph(body_text, styles["callout_body"])
    tbl  = Table([[lbl], [body]], colWidths=[content_w - 0.4 * cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg),
        ("LEFTPADDING",  (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING",   (0, 0), (0, 0),   10),
        ("BOTTOMPADDING",(0, -1),(-1, -1), 10),
        ("LINEBEFORE", (0, 0), (0, -1), 3, border),
    ]))
    return [tbl, Spacer(1, 0.3 * cm)]


def metric_card(label, value, unit="", bg=HIGHLIGHT):
    """Small metric card for a single KPI."""
    styles = build_styles()
    lbl_p = Paragraph(label, styles["callout_label"])
    val_p = Paragraph(f'<font size="20" color="{ACCENT_BLUE.hexval()}">{value}</font><font size="9" color="#555577"> {unit}</font>', styles["metric_val"])
    tbl = Table([[lbl_p], [val_p]], colWidths=[4.6 * cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg),
        ("LEFTPADDING",  (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING",   (0, 0), (0, 0),   8),
        ("BOTTOMPADDING",(0, -1),(-1, -1), 8),
        ("BOX", (0,0),(-1,-1), 1, ACCENT_BLUE),
        ("ROUNDEDCORNERS",(0,0),(-1,-1),[4,4,4,4]),
        ("ALIGN", (0, 1), (0, 1), "CENTER"),
    ]))
    return tbl


def get_plot(filename, width=None, max_w=None):
    """Return an Image flowable for a plot, scaled to fit."""
    path = os.path.join(PLOTS_DIR, filename)
    if not os.path.exists(path):
        return None
    content_w = PAGE_W - 2 * MARGIN
    target_w  = width or (max_w or content_w)
    try:
        img = PILImage.open(path)
        iw, ih = img.size
        ratio = ih / iw
        rw = min(target_w, content_w)
        rh = rw * ratio
        return Image(path, width=rw, height=rh)
    except Exception:
        return None


def bullet(text, styles, sub=False):
    s = "sub_bullet" if sub else "bullet"
    prefix = "◦" if sub else "•"
    return Paragraph(f"{prefix}  {text}", styles[s])


# ── Cover Page ────────────────────────────────────────────────────────────────
def cover_page(styles):
    content_w = PAGE_W - 2 * MARGIN
    elements  = []

    # Navy background block (simulated via a big dark table)
    cover_tbl = Table(
        [
            [Spacer(1, 1.8 * cm)],
            [Paragraph("SmartChem", styles["cover_title"])],
            [Paragraph("AI-Powered Drug Discovery Platform", styles["cover_sub"])],
            [Spacer(1, 0.4 * cm)],
            [HRFlowable(width=8 * cm, thickness=2, color=ACCENT_TEAL, hAlign="CENTER")],
            [Spacer(1, 0.4 * cm)],
            [Paragraph("ML Component — Presentation Narrative", styles["cover_sub"])],
            [Spacer(1, 0.6 * cm)],
            [Paragraph("Variational Autoencoder · Property Predictor · Latent Space Optimization", styles["cover_tag"])],
            [Spacer(1, 2.0 * cm)],
            [Paragraph("Final Year Project  ·  2025–2026", styles["cover_tag"])],
            [Spacer(1, 0.3 * cm)],
        ],
        colWidths=[content_w],
        hAlign="CENTER",
    )
    cover_tbl.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, -1), DEEP_NAVY),
        ("LEFTPADDING",   (0, 0), (-1, -1), 0),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 0),
        ("TOPPADDING",    (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
        ("ALIGN",         (0, 0), (-1, -1), "CENTER"),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ]))
    elements.append(cover_tbl)

    elements.append(Spacer(1, 0.6 * cm))

    # Subtitle / meta row
    meta_data = [
        ["Dataset", "ZINC 250K"],
        ["Model", "Conv-GRU VAE  +  MLP Predictor"],
        ["Framework", "PyTorch · RDKit · SELFIES"],
        ["Epochs", "40 (VAE)  ·  40 (Predictor)"],
        ["Prepared by", "ML Sub-team"],
    ]
    meta_table = Table(meta_data, colWidths=[4 * cm, content_w - 4 * cm])
    meta_table.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (0, -1), colors.HexColor("#E3F2FD")),
        ("BACKGROUND",    (1, 0), (1, -1), WHITE),
        ("FONTNAME",      (0, 0), (0, -1), "Helvetica-Bold"),
        ("FONTNAME",      (1, 0), (1, -1), "Helvetica"),
        ("FONTSIZE",      (0, 0), (-1, -1), 10),
        ("TEXTCOLOR",     (0, 0), (0, -1), ACCENT_BLUE),
        ("TEXTCOLOR",     (1, 0), (1, -1), DEEP_NAVY),
        ("LEFTPADDING",   (0, 0), (-1, -1), 10),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 10),
        ("TOPPADDING",    (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("GRID",          (0, 0), (-1, -1), 0.5, colors.HexColor("#C8D8E8")),
        ("ROWBACKGROUNDS",(0, 0), (-1, -1), [LIGHT_GRAY, WHITE]),
    ]))
    elements.append(meta_table)
    elements.append(PageBreak())
    return elements


# ── Section 1 — Model & Learning Algorithm ───────────────────────────────────
def section_model(styles):
    el = []
    el += section_banner("SECTION 01", "The Model & Learning Algorithm — Why We Built What We Built", styles)

    el.append(Paragraph(
        "The core machine learning challenge in SmartChem is <b>conditional molecular generation</b>: "
        "given desired drug-like properties (QED, LogP, SAS), produce novel, chemically valid molecules. "
        "This is a fundamentally generative task operating on a discrete, structured chemical language — "
        "a problem space that demands a very particular class of model.",
        styles["body"]
    ))

    # -- Why a VAE? --
    el += section_banner("1.1", "Architecture Choice: Convolutional–GRU Variational Autoencoder (VAE)", styles, ACCENT_BLUE)

    el.append(Paragraph(
        "We chose a <b>Variational Autoencoder (VAE)</b> as our generative backbone. The key insight is that "
        "a VAE learns a <i>continuous, smooth latent space</i> of chemical structures. Once encoded, any point "
        "in this latent space can be decoded back into a molecule — and neighbouring points decode into "
        "structurally similar molecules. This property is essential for the optimization task downstream.",
        styles["body"]
    ))

    decision_data = [
        ["Architecture", "Pros", "Cons", "Verdict"],
        ["VAE\n(our choice)", "Smooth latent space; gradient-optimizable;\nstable under KL annealing", "Reconstruction loss can be high early", "✅ CHOSEN"],
        ["GAN", "Sharp outputs", "Mode collapse; unstable training;\nno latent inference", "❌ Rejected"],
        ["Transformer\n(SELFIES-GPT)", "State-of-the-art language modeling", "No latent space; slow; huge data need", "❌ Rejected"],
        ["LSTM-AE", "Sequence-native", "No variational regularization;\npoor interpolation", "❌ Rejected"],
        ["Flow Model", "Exact likelihood", "Computationally heavy; complex to train", "❌ Rejected"],
    ]
    decision_tbl = Table(
        decision_data,
        colWidths=[3.5 * cm, 5.8 * cm, 4.5 * cm, 2.2 * cm],
        hAlign="LEFT",
    )
    decision_tbl.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, 0), DEEP_NAVY),
        ("TEXTCOLOR",     (0, 0), (-1, 0), WHITE),
        ("FONTNAME",      (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME",      (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE",      (0, 0), (-1, -1), 8.5),
        ("BACKGROUND",    (0, 1), (-1, 1), colors.HexColor("#E8F5E9")),
        ("ROWBACKGROUNDS",(0, 2), (-1, -1), [LIGHT_GRAY, WHITE]),
        ("GRID",          (0, 0), (-1, -1), 0.5, colors.HexColor("#BBCCDD")),
        ("LEFTPADDING",   (0, 0), (-1, -1), 7),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 7),
        ("TOPPADDING",    (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("FONTNAME",      (0, 1), (0, 1), "Helvetica-Bold"),
        ("TEXTCOLOR",     (3, 1), (3, 1), GREEN_OK),
        ("TEXTCOLOR",     (3, 2), (3, -1), colors.red),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ]))
    el.append(decision_tbl)
    el.append(Spacer(1, 0.3 * cm))

    el += callout_box(
        "The Decisive Parameter: Latent Space Smoothness",
        "The single most important selection criterion was whether the model produces a smooth, "
        "differentiable latent space. A smooth latent space is what enables gradient-based optimization "
        "(and future Bayesian optimization) to traverse chemical space from a seed molecule toward a "
        "user-specified property target. GANs and plain autoencoders fail this test — only VAEs and "
        "Normalizing Flows offer this guarantee, and VAEs are significantly cheaper to train.",
        styles,
    )

    el += section_banner("1.2", "Encoder: 1D Convolutional Stack", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The encoder uses three stacked <b>1D convolutional layers</b> (kernel sizes 9, 9, 11; channels 9, 9, 10) "
        "to capture local chemical grammar patterns (e.g. ring closures, branch tokens) from the SELFIES "
        "token sequence. An <b>Adaptive Max Pool</b> collapses the variable-length output to a fixed-size "
        "feature vector, which is then projected to μ (mean) and log σ² (log-variance) via two separate "
        "linear heads — the parameters of the variational posterior q(z|x).",
        styles["body"]
    ))
    el.append(bullet("Embedding Dim: 128  |  Conv Channels: 9 → 9 → 10", styles))
    el.append(bullet("Kernel Sizes: 9, 9, 11 — captures short- and medium-range chemical motifs", styles))
    el.append(bullet("Latent Dim: 128 — balances expressiveness with training stability", styles))

    el += section_banner("1.3", "Decoder: 3-Layer GRU with Temperature Sampling", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The decoder is a <b>3-layer GRU (Gated Recurrent Unit)</b> that autoregressively reconstructs "
        "the SELFIES token sequence from the latent vector z. The GRU was preferred over an LSTM because "
        "it converges faster with fewer parameters while achieving comparable sequence modelling quality "
        "at the token lengths we operate at (max_len=100). Temperature sampling (T=1.5 for random "
        "generation, T=0.8 for targeted decoding) controls the diversity/quality trade-off.",
        styles["body"]
    ))
    el.append(bullet("Hidden Dim: 256  |  Layers: 3  |  Output: Softmax over SELFIES vocab", styles))
    el.append(bullet("Temperature 1.5 → exploration (random generation)", styles))
    el.append(bullet("Temperature 0.8 → exploitation (targeted optimization)", styles))

    el += section_banner("1.4", "Property Predictor: MLP on Latent Space", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "A lightweight <b>Multi-Layer Perceptron (MLP)</b> is attached to the latent space to predict "
        "three drug-likeness properties simultaneously: <b>QED</b> (Quantitative Estimate of Drug-likeness), "
        "<b>LogP</b> (lipophilicity), and <b>SAS</b> (Synthetic Accessibility Score, approximated via "
        "ring + rotatable bond count). The predictor takes the VAE's μ vector as input, making it "
        "gradient-compatible with the encoder — a design that enables end-to-end optimization.",
        styles["body"]
    ))
    el.append(bullet("Architecture: Linear(128→64) → ReLU → Linear(64→32) → ReLU → Linear(32→3)", styles))
    el.append(bullet("Weighted MSE loss: QED weight=10× (most clinically relevant)", styles))
    el.append(bullet("Outputs are directly used to score and rank generated molecules", styles))

    el.append(PageBreak())
    return el


# ── Section 2 — Dataset & Preprocessing ─────────────────────────────────────
def section_dataset(styles):
    el = []
    el += section_banner("SECTION 02", "Dataset Sourcing & Preprocessing Pipeline", styles)

    el.append(Paragraph(
        "A generative model is only as good as the chemical distribution it has learned. "
        "The dataset choice is therefore a foundational decision — it determines the vocabulary, "
        "the structural diversity, and ultimately the validity rate of generated molecules.",
        styles["body"]
    ))

    el += section_banner("2.1", "Source: ZINC 250K", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "<b>ZINC250K</b> is the de-facto benchmark dataset for molecular deep learning. Curated from "
        "the ZINC database — a commercialy available compound library used in virtual screening — "
        "it contains exactly <b>250,000 drug-like molecules</b>, each pre-filtered to satisfy "
        "Lipinski's Rule of Five (Ro5). This makes it ideal for training a model intended to generate "
        "drug candidates rather than arbitrary organic molecules.",
        styles["body"]
    ))

    zinc_data = [
        ["Property", "Value"],
        ["Total molecules", "250,000"],
        ["Source", "ZINC database (Irwin & Shoichet, 2005)"],
        ["Filter applied", "Lipinski Rule of Five (drug-like)"],
        ["Represented as", "SMILES strings + pre-computed SELFIES"],
        ["Training subset used", "150,000 molecules (MVP_LIMIT)"],
        ["Post-filter subset", "~120,000 (length filter: SMILES + SELFIES < 100 chars)"],
    ]
    zt = Table(zinc_data, colWidths=[5.5 * cm, PAGE_W - 2 * MARGIN - 5.5 * cm])
    zt.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, 0), DEEP_NAVY),
        ("TEXTCOLOR",    (0, 0), (-1, 0), WHITE),
        ("FONTNAME",     (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME",     (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTNAME",     (1, 1), (1, -1), "Helvetica"),
        ("FONTSIZE",     (0, 0), (-1, -1), 9.5),
        ("TEXTCOLOR",    (0, 1), (0, -1), ACCENT_BLUE),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[LIGHT_GRAY, WHITE]),
        ("GRID",         (0, 0), (-1, -1), 0.5, colors.HexColor("#BBCCDD")),
        ("LEFTPADDING",  (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING",   (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 7),
    ]))
    el.append(zt)
    el.append(Spacer(1, 0.4 * cm))

    el += section_banner("2.2", "Why SELFIES, Not SMILES?", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "We use <b>SELFIES (SELF-referencing Embedded Strings)</b> as our primary molecular representation, "
        "rather than raw SMILES strings. SMILES are fragile — even a single character error produces an "
        "invalid molecule. SELFIES, by construction, <i>guarantees syntactic validity of every string</i>, "
        "regardless of how they are generated or mutated. This is critical in a deep learning setting: "
        "the decoder's outputs are always valid SELFIES, even before any chemical filtering.",
        styles["body"]
    ))
    el.append(bullet("SMILES validity: syntax-sensitive — 1 error = invalid molecule", styles))
    el.append(bullet("SELFIES validity: 100% syntactically valid by design → higher generation validity rate", styles))
    el.append(bullet("Both vocabularies are built — parallel SMILES tensors are retained for cross-validation", styles))

    el += section_banner("2.3", "Preprocessing Pipeline (misc/preprocess.py)", styles, ACCENT_BLUE)

    pipeline_steps = [
        ("Step 1 — Load Raw CSV", "Read the top 150,000 rows from train_molecules.csv (ZINC 250K download). Using nrows=MVP_LIMIT keeps preprocessing fast and reproducible."),
        ("Step 2 — Length Filter", "Discard molecules where either the SMILES or SELFIES string exceeds 100 characters (max_len=100). This removes rare, large molecules that would significantly inflate training time and vocabulary size while contributing marginally to learning drug-like space."),
        ("Step 3 — Vocabulary Construction", "Scan all sequences to build token dictionaries. SELFIES are split by bracket notation ([C], [Branch1], …) while SMILES are split character-by-character. Special tokens <pad>=0, <sos>=1, <eos>=2 are prepended; all real tokens start at index 3."),
        ("Step 4 — Tensor Encoding", "Convert each molecule string to a fixed-length integer tensor of length 100. Sequences shorter than 100 are right-padded with <pad>; longer sequences are truncated (rare after the length filter)."),
        ("Step 5 — Artifact Saving", "Serialize vocab dictionaries as JSON (selfies_vocab.json, smiles_vocab.json) and tensor sequences as PyTorch .pt files (train_selfies.pt, train_smiles.pt). The backend loads these at inference time to decode generated tokens."),
    ]
    for label, desc in pipeline_steps:
        el += callout_box(label, desc, styles, bg=LIGHT_GRAY, border=ACCENT_TEAL)

    el.append(PageBreak())
    return el


# ── Section 3 — Training ─────────────────────────────────────────────────────
def section_training(styles):
    el = []
    el += section_banner("SECTION 03", "Training: Configuration, Loss Functions & Convergence", styles)

    el.append(Paragraph(
        "Training was performed in two sequential stages: the VAE was trained first to learn the chemical "
        "latent space, and then the Property Predictor was trained on the frozen VAE's encoded representations. "
        "This staged approach prevents the predictor from interfering with the VAE's reconstruction objective "
        "during the critical early epochs.",
        styles["body"]
    ))

    el += section_banner("3.1", "Training Configuration", styles, ACCENT_BLUE)

    config_data = [
        ["Hyperparameter",     "VAE",               "Property Predictor"],
        ["Optimizer",          "Adam",               "Adam"],
        ["Learning Rate",      "5×10⁻⁴",             "1×10⁻³"],
        ["Batch Size",         "64",                 "256"],
        ["Epochs",             "40",                 "40"],
        ["LR Scheduler",       "ReduceLROnPlateau\n(patience=3, factor=0.5)", "—"],
        ["Gradient Clipping",  "L2 norm ≤ 5.0",      "—"],
        ["Device",             "GPU (CUDA) / CPU",   "GPU (CUDA) / CPU"],
        ["Latent Dim",         "128",                 "128 (input)"],
        ["Max Sequence Len",   "100 tokens",          "—"],
    ]
    ct = Table(config_data, colWidths=[5.2 * cm, 4.8 * cm, 5.0 * cm])
    ct.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, 0), DEEP_NAVY),
        ("TEXTCOLOR",    (0, 0), (-1, 0), WHITE),
        ("FONTNAME",     (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME",     (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTNAME",     (1, 1), (-1, -1),"Helvetica"),
        ("FONTSIZE",     (0, 0), (-1, -1), 9),
        ("TEXTCOLOR",    (0, 1), (0, -1), ACCENT_BLUE),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[LIGHT_GRAY, WHITE]),
        ("GRID",         (0, 0), (-1, -1), 0.5, colors.HexColor("#BBCCDD")),
        ("LEFTPADDING",  (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING",   (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 6),
        ("VALIGN",       (0, 0), (-1, -1), "MIDDLE"),
    ]))
    el.append(ct)
    el.append(Spacer(1, 0.4 * cm))

    el += section_banner("3.2", "VAE Loss Function: Reconstruction + KL Divergence", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The VAE is trained with the standard ELBO (Evidence Lower BOund) objective, decomposed into "
        "two terms:",
        styles["body"]
    ))
    el.append(bullet(
        "<b>Reconstruction Loss (BCE/Cross-Entropy):</b> Measures how faithfully the decoder "
        "reconstructs the input SELFIES token sequence. Computed as per-token cross-entropy summed "
        "across the sequence, normalized per sample. This drives the model to learn valid chemical grammar.",
        styles
    ))
    el.append(bullet(
        "<b>KL Divergence:</b> Regularizes the latent space by penalizing posteriors q(z|x) that "
        "deviate from the standard Gaussian prior N(0,I). This ensures the latent space is smooth "
        "and continuous — the core property enabling interpolation and optimization.",
        styles
    ))
    el += callout_box(
        "KL Annealing — Critical Design Decision",
        "A well-known failure mode of VAEs is 'posterior collapse', where the KL term overwhelms the "
        "reconstruction loss early in training, forcing the model to ignore the input. We mitigate this "
        "with KL annealing: the KL weight is linearly ramped from 0 → 1 over the first 4 epochs "
        "(kld_weight = min(1.0, epoch/4.0)). This lets the model first learn to reconstruct, "
        "then gradually regularize its latent space.",
        styles,
        bg=colors.HexColor("#FFF8E1"),
        border=WARN_AMBER,
    )

    el += section_banner("3.3", "VAE Loss Convergence (from vae_training_log.csv)", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The VAE was trained for 40 epochs. Key observations from the training log:",
        styles["body"]
    ))

    # Key metrics extracted from log
    el.append(bullet("Epoch 1 Total Loss: 1.5624  →  Epoch 40 Total Loss: 0.7429  (−52.5% reduction)", styles))
    el.append(bullet("BCE Loss: 1.55 → 0.631  — reconstruction quality improved steadily", styles))
    el.append(bullet("KL Loss: 0.012 → 0.112  — KL rose gradually as annealing took effect, then settled as regularization balanced reconstruction", styles))
    el.append(bullet("No loss spikes observed — gradient clipping (≤5.0) ensured stable training throughout", styles))
    el.append(Spacer(1, 0.25 * cm))

    vae_plot = get_plot("vae_loss.png")
    if vae_plot:
        el.append(vae_plot)
        el.append(Paragraph("Figure 1 — VAE training loss (BCE, KL, Total) over 40 epochs.", styles["caption"]))

    el += section_banner("3.4", "Predictor Loss Convergence (from predictor_log.csv)", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The Property Predictor was trained for 40 epochs on VAE-encoded latent vectors (μ) with "
        "RDKit-computed true labels. The MSE per property:",
        styles["body"]
    ))
    el.append(bullet("MSE QED:  0.0783 → 0.0263  — QED prediction improved by 66%", styles))
    el.append(bullet("MSE LogP: 2.00   → 0.856   — LogP prediction improved by 57%", styles))
    el.append(bullet("MSE SAS:  3.50   → 1.723   — SAS proxy prediction improved by 51%", styles))
    el.append(bullet("QED convergence (lowest absolute MSE) is the most important metric — it is weighted 10× in the combined loss", styles))
    el.append(Spacer(1, 0.25 * cm))

    pred_plot = get_plot("predictor_training_loss.png")
    if pred_plot:
        el.append(pred_plot)
        el.append(Paragraph("Figure 2 — Property Predictor training MSE over 40 epochs (QED, LogP, SAS).", styles["caption"]))

    el.append(PageBreak())
    return el


# ── Section 4 — Current Performance ─────────────────────────────────────────
def section_performance(styles):
    el = []
    el += section_banner("SECTION 04", "Current Model Performance", styles)

    el.append(Paragraph(
        "Performance is evaluated across three dimensions: <b>generation validity</b> "
        "(are the outputs real molecules?), <b>property prediction accuracy</b> "
        "(how well does the predictor approximate RDKit ground truth?), and "
        "<b>optimization effectiveness</b> (can the system navigate the latent space "
        "toward a target?). All numbers below are from the evaluation logs and plots "
        "in evaluation/logs/ and evaluation/plots/.",
        styles["body"]
    ))

    el += section_banner("4.1", "Generation Validity (from validity_stats.json)", styles, ACCENT_BLUE)

    # Metric cards row
    mc1 = metric_card("SELFIES Validity", "94.7%", "of 1,000 generated", bg=colors.HexColor("#E3F2FD"))
    mc2 = metric_card("RDKit Pass Rate", "41.2%", "chemical validity", bg=colors.HexColor("#E8F5E9"))
    mc3 = metric_card("Lipinski Pass", "28.4%", "drug-like (Ro5)", bg=colors.HexColor("#FFF3E0"))
    cards_tbl = Table([[mc1, Spacer(0.3*cm, 1), mc2, Spacer(0.3*cm, 1), mc3]], hAlign="LEFT")
    cards_tbl.setStyle(TableStyle([
        ("LEFTPADDING",  (0,0),(-1,-1), 0),
        ("RIGHTPADDING", (0,0),(-1,-1), 0),
        ("VALIGN",       (0,0),(-1,-1), "TOP"),
    ]))
    el.append(cards_tbl)
    el.append(Spacer(1, 0.35 * cm))

    el += callout_box(
        "Interpreting the Numbers",
        "The 94.7% SELFIES validity is a strong baseline — it means the GRU decoder reliably follows "
        "SELFIES grammar even when decoding from randomly sampled latent points. The 41.2% RDKit pass "
        "rate reflects full chemical valence and ring-closure validity, which is a strict secondary filter. "
        "The 28.4% Lipinski rate shows that roughly 1 in 3.5 valid molecules also satisfies drug-likeness "
        "constraints — a ratio consistent with the literature for ZINC-trained VAEs at this training scale.",
        styles,
    )

    validity_plot = get_plot("validity_stats.png")
    if validity_plot:
        el.append(validity_plot)
        el.append(Paragraph("Figure 3 — Validity funnel: SELFIES validity → RDKit pass → Lipinski pass.", styles["caption"]))

    el += section_banner("4.2", "Property Prediction Accuracy", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The scatter plots (predictor_scatter.png) compare the predictor's output against "
        "RDKit ground truth across 500 test molecules. The closer the points to the diagonal y=x line, "
        "the better. Key takeaways from the final epoch:",
        styles["body"]
    ))
    el.append(bullet("QED: Final MSE = 0.026 — predictions cluster tightly around the diagonal; "
                     "the model reliably distinguishes high-QED (>0.7) from low-QED (<0.3) molecules", styles))
    el.append(bullet("LogP: Final MSE = 0.856 — acceptable for a molecule-level property; "
                     "spread is wider at extreme LogP values (< −1 or > 4)", styles))
    el.append(bullet("SAS: Final MSE = 1.72 — higher variance expected as SAS is only approximated "
                     "by ring+bond count rather than the full Ertl algorithm", styles))
    el.append(Spacer(1, 0.2 * cm))

    scatter_plot = get_plot("predictor_scatter.png")
    if scatter_plot:
        el.append(scatter_plot)
        el.append(Paragraph("Figure 4 — True vs. Predicted scatter plot for QED, LogP, and SAS (500 test samples).", styles["caption"]))

    el += section_banner("4.3", "Latent Space Quality", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The latent space visualization (latent_space.png) shows a 2D PCA projection of 500 molecules "
        "encoded by the VAE. A well-trained VAE should produce a roughly Gaussian blob — molecules should "
        "cluster by structural similarity, with no large empty voids (posterior collapse) or sharp isolated "
        "clusters (no regularization).",
        styles["body"]
    ))

    latent_plot = get_plot("latent_space.png")
    if latent_plot:
        el.append(latent_plot)
        el.append(Paragraph("Figure 5 — PCA projection of the 128-D latent space (coloured by QED score).", styles["caption"]))

    el += section_banner("4.4", "Optimization Convergence", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "The optimization log traces 75 gradient-descent steps in latent space toward targets "
        "QED=0.8 / LogP=2.0 / SAS=2.0. Key observations:",
        styles["body"]
    ))
    el.append(bullet("Starting QED: 0.413  →  Converged QED: 0.82  (+98.5% improvement)", styles))
    el.append(bullet("LogP tracks target closely: starts at 2.003 → stabilizes near 2.22", styles))
    el.append(bullet("L2 Distance from start: 0.001 → ~0.095 — moderate latent travel, "
                     "ensuring generated molecules are structurally related to the seed", styles))
    el.append(bullet("QED plateau observed at step ~50 — convergence achieved well within 75 steps", styles))
    el.append(Spacer(1, 0.2 * cm))

    opt_plot = get_plot("optimization_curve.png")
    if opt_plot:
        el.append(opt_plot)
        el.append(Paragraph("Figure 6 — Optimization curve: QED, LogP, and latent L2-distance over 75 gradient steps.", styles["caption"]))

    recon_plot = get_plot("reconstruction_error.png")
    if recon_plot:
        el.append(Spacer(1, 0.15 * cm))
        el.append(recon_plot)
        el.append(Paragraph("Figure 7 — Reconstruction error distribution across the test set.", styles["caption"]))

    el.append(PageBreak())
    return el


# ── Section 5 — Research Paper ───────────────────────────────────────────────
def section_research(styles):
    el = []
    el += section_banner("SECTION 05", "Research Paper: Gradient vs. Bayesian Optimization in Chemical Latent Space", styles)

    el.append(Paragraph(
        "The SmartChem system is not just an engineering project — it is the experimental platform "
        "for a focused research investigation into <b>how the choice of optimization strategy affects "
        "the quality and diversity of molecules found in the latent space</b>.",
        styles["body"]
    ))

    el += section_banner("5.1", "Research Question", styles, ACCENT_BLUE)
    el += callout_box(
        "Core Research Question",
        "Does replacing gradient-based latent space optimization with Bayesian Optimization "
        "(specifically, Naive Bayes / Gaussian Process surrogate models) produce molecules with "
        "superior drug-likeness metrics, and does it explore a more diverse region of chemical space "
        "compared to greedy gradient descent?",
        styles,
        bg=colors.HexColor("#EDE7F6"),
        border=colors.HexColor("#7B1FA2"),
    )

    el += section_banner("5.2", "Motivation", styles, ACCENT_BLUE)
    el.append(Paragraph(
        "Current latent space optimization methods in molecular VAEs almost universally rely on "
        "<b>gradient-based approaches</b> — specifically backpropagating a property loss signal "
        "through the predictor back into the latent vector z. While effective at hill-climbing "
        "toward a single target, gradient descent has two fundamental limitations in this setting:",
        styles["body"]
    ))
    el.append(bullet("It follows the gradient greedily — it converges to the nearest local optimum "
                     "and cannot escape saddle points or explore multiple modes of the property distribution simultaneously.", styles))
    el.append(bullet("It requires differentiable objectives — any scoring function that is not "
                     "backprop-compatible (e.g. RDKit-computed SAS) must be approximated via the predictor, "
                     "introducing an additional source of error.", styles))
    el.append(Spacer(1, 0.2 * cm))
    el.append(Paragraph(
        "<b>Bayesian Optimization (BO)</b> treats the latent space as a black-box function and builds "
        "a probabilistic surrogate model (Gaussian Process or Naive Bayes classifier) to model the "
        "property landscape. It balances exploitation (move toward high-scoring regions) with "
        "<b>exploration</b> (sample uncertain regions via the acquisition function), making it "
        "inherently better suited for navigating a complex, multi-modal chemical space.",
        styles["body"]
    ))

    el += section_banner("5.3", "Proposed Experimental Design", styles, ACCENT_BLUE)
    comparison_data = [
        ["Aspect",                "Gradient Descent (Baseline)",      "Bayesian Optimization (Proposed)"],
        ["Optimization target",   "Predictor MSE Loss",               "RDKit / Predictor score (black-box)"],
        ["Landscape model",       "None (direct gradient)",           "Gaussian Process / Naive Bayes surrogate"],
        ["Exploration",           "None — greedy descent",            "UCB / EI acquisition function"],
        ["Diversity of results",  "Narrow — converges to one point",  "Broad — explicitly incentivized"],
        ["Steps required",        "75 (from optimization log)",       "Expected: 30–50 (BO is sample-efficient)"],
        ["Evaluation metrics",    "Final QED, LogP, SAS",             "Final QED, LogP, SAS + diversity (Tanimoto)"],
    ]
    comp_tbl = Table(
        comparison_data,
        colWidths=[4.2 * cm, 5.5 * cm, 6.3 * cm],
    )
    comp_tbl.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, 0), DEEP_NAVY),
        ("TEXTCOLOR",    (0, 0), (-1, 0), WHITE),
        ("FONTNAME",     (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME",     (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTNAME",     (1, 1), (-1, -1),"Helvetica"),
        ("FONTSIZE",     (0, 0), (-1, -1), 8.5),
        ("TEXTCOLOR",    (0, 1), (0, -1), ACCENT_BLUE),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[LIGHT_GRAY, WHITE]),
        ("GRID",         (0, 0), (-1, -1), 0.5, colors.HexColor("#BBCCDD")),
        ("LEFTPADDING",  (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING",   (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 6),
        ("VALIGN",       (0, 0), (-1, -1), "MIDDLE"),
    ]))
    el.append(comp_tbl)
    el.append(Spacer(1, 0.4 * cm))

    el += section_banner("5.4", "Expected Contributions", styles, ACCENT_BLUE)
    el.append(bullet(
        "<b>Empirical benchmark</b>: First systematic comparison of gradient-based vs. BO optimization "
        "in a SELFIES-VAE latent space on ZINC 250K molecules.", styles
    ))
    el.append(bullet(
        "<b>Diversity metric</b>: We will introduce Tanimoto similarity clustering as a diversity measure "
        "alongside standard drug-likeness KPIs — a dimension absent from most optimization papers.", styles
    ))
    el.append(bullet(
        "<b>Practical guidance</b>: Results will directly inform whether future SmartChem versions "
        "should replace the gradient optimizer with a BO backend for the targeted generation API.", styles
    ))
    el.append(bullet(
        "<b>Open analysis</b>: All optimization logs (optimization_log.csv) and latent trajectories "
        "will be released as supplementary data, making the comparison reproducible.", styles
    ))

    el += section_banner("5.5", "Tentative Paper Title", styles, ACCENT_BLUE)
    el += callout_box(
        "Working Title",
        '"Gradient Descent vs. Bayesian Optimization in Chemical Latent Space: '
        'A Comparative Study on Drug-Likeness and Structural Diversity Using a '
        'SELFIES Variational Autoencoder"',
        styles,
        bg=colors.HexColor("#E8F5E9"),
        border=GREEN_OK,
    )

    el.append(Spacer(1, 0.5 * cm))
    el += callout_box(
        "Why This Research Matters",
        "The pharmaceutical industry loses billions to failed drug candidates every year. Improving the "
        "efficiency and diversity of in-silico molecular exploration can directly reduce that cost. "
        "By quantifying exactly how much Bayesian optimization outperforms (or fails to outperform) "
        "gradient descent in continuous chemical latent space, this paper provides actionable evidence "
        "for practitioners building the next generation of AI drug discovery pipelines.",
        styles,
        bg=LIGHT_GRAY,
        border=ACCENT_TEAL,
    )

    el.append(PageBreak())
    return el


# ── Closing Summary ───────────────────────────────────────────────────────────
def section_summary(styles):
    el = []
    el += section_banner("SUMMARY", "Narrative at a Glance — Key Takeaways", styles)

    summary_data = [
        ["#", "Question", "Answer"],
        ["1", "What model + algo?",
         "Conv-GRU VAE (latent_dim=128) trained with ELBO + KL annealing; "
         "MLP Predictor on latent vectors for QED/LogP/SAS prediction."],
        ["2", "Why this model?",
         "VAE is the only practical architecture that provides a smooth, "
         "differentiable latent space — the sine qua non for gradient-based "
         "(and future BO) chemical optimization."],
        ["3", "Dataset?",
         "ZINC 250K — 150,000 drug-like molecules preprocessed into SELFIES "
         "tensors; length-filtered to max_len=100; vocab built from scratch."],
        ["4", "How was training measured?",
         "BCE loss + KL divergence (VAE); per-property MSE for QED, LogP, SAS "
         "(Predictor). All metrics logged per epoch to CSV and plotted."],
        ["5", "Current performance?",
         "94.7% SELFIES validity; 41.2% RDKit pass; 28.4% Lipinski; "
         "QED MSE=0.026; optimization converges QED from 0.41→0.82 in 50 steps."],
        ["6", "Research paper?",
         "Gradient Descent vs. Bayesian Optimization in Chemical Latent Space — "
         "a comparative study on drug-likeness and structural diversity."],
    ]
    st = Table(summary_data, colWidths=[0.8*cm, 4.5*cm, PAGE_W - 2*MARGIN - 5.3*cm])
    st.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, 0), DEEP_NAVY),
        ("TEXTCOLOR",    (0, 0), (-1, 0), WHITE),
        ("FONTNAME",     (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME",     (0, 1), (1, -1), "Helvetica-Bold"),
        ("FONTNAME",     (2, 1), (2, -1), "Helvetica"),
        ("FONTSIZE",     (0, 0), (-1, -1), 9),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[LIGHT_GRAY, WHITE]),
        ("GRID",         (0, 0), (-1, -1), 0.5, colors.HexColor("#BBCCDD")),
        ("LEFTPADDING",  (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING",   (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 7),
        ("VALIGN",       (0, 0), (-1, -1), "TOP"),
        ("TEXTCOLOR",    (0, 1), (0, -1), ACCENT_TEAL),
        ("TEXTCOLOR",    (1, 1), (1, -1), ACCENT_BLUE),
    ]))
    el.append(st)
    el.append(Spacer(1, 0.8*cm))

    el.append(HRFlowable(width="100%", thickness=1.5, color=ACCENT_BLUE))
    el.append(Spacer(1, 0.3*cm))
    el.append(Paragraph(
        "SmartChem represents a complete, end-to-end AI-powered drug discovery pipeline — "
        "from raw ZINC data to a deployed web application capable of generating and optimizing "
        "novel drug candidates in real time. The ML component, detailed in this document, "
        "is both the technical heart of the project and the foundation for ongoing academic research.",
        styles["quote"]
    ))

    return el


# ── Main Build ────────────────────────────────────────────────────────────────
def build_pdf():
    doc = SimpleDocTemplate(
        OUTPUT_PATH,
        pagesize=A4,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN + 0.8 * cm,   # extra space for header bar
        bottomMargin=MARGIN + 0.5 * cm, # extra space for footer bar
        title="SmartChem ML Presentation Narrative",
        author="SmartChem ML Sub-team",
        subject="Final Year Project — ML Component Presentation",
    )

    styles = build_styles()
    story  = []

    story += cover_page(styles)
    story += section_model(styles)
    story += section_dataset(styles)
    story += section_training(styles)
    story += section_performance(styles)
    story += section_research(styles)
    story += section_summary(styles)

    doc.build(story, canvasmaker=BrandedCanvas)
    print(f"\n✅ PDF generated successfully:\n   {OUTPUT_PATH}\n")


if __name__ == "__main__":
    build_pdf()
