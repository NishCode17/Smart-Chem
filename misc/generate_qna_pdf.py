"""
generate_qna_pdf.py
===================
Generates a polished PDF for the QnA preparation related to the ML Presentation.

Run from project root:
    python misc/generate_qna_pdf.py
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_PATH = os.path.join(ROOT, "SmartChem_ML_QnA_Prep.pdf")

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.units import cm, mm
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, KeepTogether, Table, TableStyle, HRFlowable, PageBreak
    )
    from reportlab.pdfgen import canvas
except ImportError:
    print("ReportLab is not installed. Run 'pip install reportlab first'.")
    sys.exit(1)

# ── Color Palette ─────────────────────────────────────────────────────────────
DEEP_NAVY   = colors.HexColor("#0D1B2A")
ACCENT_BLUE = colors.HexColor("#1E88E5")
ACCENT_TEAL = colors.HexColor("#00BCD4")
LIGHT_GRAY  = colors.HexColor("#F5F7FA")
WHITE       = colors.white
MED_GRAY    = colors.HexColor("#8899AA")
HIGHLIGHT   = colors.HexColor("#E3F2FD")

PAGE_W, PAGE_H = A4
MARGIN = 2.0 * cm

# ── Custom Canvas ─────────────────────────────────────────────────────────────
class QnACanvas(canvas.Canvas):
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

        if page_num == 1:
            return

        self.saveState()
        self.setFillColor(ACCENT_BLUE)
        self.rect(0, PAGE_H - 8 * mm, PAGE_W, 8 * mm, stroke=0, fill=1)
        self.setFillColor(WHITE)
        self.setFont("Helvetica-Bold", 8)
        self.drawString(MARGIN, PAGE_H - 5.5 * mm, "SmartChem | ML Viva & Presentation Q&A")
        
        self.setFillColor(DEEP_NAVY)
        self.rect(0, 0, PAGE_W, 10 * mm, stroke=0, fill=1)
        self.setFillColor(MED_GRAY)
        self.setFont("Helvetica", 7.5)
        self.drawString(MARGIN, 3.5 * mm, "Final Year Project  ·  SmartChem Internal Review")
        self.setFillColor(WHITE)
        self.drawRightString(PAGE_W - MARGIN, 3.5 * mm, f"Page {page_num} / {page_count}")
        self.restoreState()


# ── Style Sheet ───────────────────────────────────────────────────────────────
def build_styles():
    styles = {}
    
    styles["cover_title"] = ParagraphStyle(
        "cover_title", fontName="Helvetica-Bold", fontSize=26, textColor=WHITE,
        alignment=TA_CENTER, spaceAfter=6, leading=32
    )
    styles["cover_sub"] = ParagraphStyle(
        "cover_sub", fontName="Helvetica", fontSize=14, textColor=colors.HexColor("#B0C4D8"),
        alignment=TA_CENTER, spaceAfter=4, leading=18
    )
    styles["cover_tag"] = ParagraphStyle(
        "cover_tag", fontName="Helvetica-Oblique", fontSize=10, textColor=ACCENT_TEAL,
        alignment=TA_CENTER, spaceAfter=2
    )
    
    styles["q_num"] = ParagraphStyle(
        "q_num", fontName="Helvetica-Bold", fontSize=11, textColor=ACCENT_TEAL,
        spaceAfter=2
    )
    styles["question"] = ParagraphStyle(
        "question", fontName="Helvetica-Bold", fontSize=12, textColor=DEEP_NAVY,
        spaceAfter=6, leading=16
    )
    styles["answer_body"] = ParagraphStyle(
        "answer_body", fontName="Helvetica", fontSize=10.5, textColor=colors.HexColor("#1A1A2E"),
        alignment=TA_JUSTIFY, leading=16, spaceAfter=8
    )
    styles["bullet"] = ParagraphStyle(
        "bullet", fontName="Helvetica", fontSize=10.5, textColor=colors.HexColor("#1A1A2E"),
        leftIndent=15, bulletIndent=5, spaceAfter=6, leading=16, bulletFontName="Helvetica-Bold",
        bulletFontSize=12, bulletColor=ACCENT_BLUE
    )
    styles["section_header"] = ParagraphStyle(
        "section_header", fontName="Helvetica-Bold", fontSize=16, textColor=WHITE,
        alignment=TA_CENTER
    )
    
    return styles

def section_banner(title, styles):
    content_w = PAGE_W - 2 * MARGIN
    p = Paragraph(title, styles["section_header"])
    tbl = Table([[p]], colWidths=[content_w])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), DEEP_NAVY),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("ROUNDEDCORNERS", (0, 0), (-1, -1), [4, 4, 4, 4]),
    ]))
    return [tbl, Spacer(1, 0.5 * cm)]

def qa_block(num, question, answer_paragraphs, styles):
    blocks = []
    
    # Question Card
    content_w = PAGE_W - 2 * MARGIN
    q_num_p = Paragraph(f"QUESTION {num}", styles["q_num"])
    q_para = Paragraph(question, styles["question"])
    
    q_table = Table([[q_num_p], [q_para]], colWidths=[content_w])
    q_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), HIGHLIGHT),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING", (0, 0), (0, 0), 10),
        ("BOTTOMPADDING", (0, 1), (0, 1), 10),
        ("LINEBEFORE", (0, 0), (-1, -1), 4, ACCENT_BLUE),
    ]))
    
    blocks.append(q_table)
    blocks.append(Spacer(1, 0.3 * cm))
    
    # Answer paragraphs
    for p in answer_paragraphs:
        if isinstance(p, str):
            if p.startswith("• "):
                blocks.append(Paragraph(f"• {p[2:]}", styles["bullet"]))
            else:
                blocks.append(Paragraph(p, styles["answer_body"]))
        else:
            blocks.append(p)
            
    blocks.append(Spacer(1, 0.7 * cm))
    return [KeepTogether(blocks)]

# ── Content ───────────────────────────────────────────────────────────────────
def build_story():
    styles = build_styles()
    story = []
    
    # Cover Page
    content_w = PAGE_W - 2 * MARGIN
    cover_tbl = Table(
        [
            [Spacer(1, 4 * cm)],
            [Paragraph("SmartChem", styles["cover_title"])],
            [Paragraph("ML Technical Defense & Q&A Prep", styles["cover_sub"])],
            [Spacer(1, 0.4 * cm)],
            [HRFlowable(width=8 * cm, thickness=2, color=ACCENT_TEAL, hAlign="CENTER")],
            [Spacer(1, 0.4 * cm)],
            [Paragraph("Comprehensive Question Bank for Viva and Presentation", styles["cover_sub"])],
            [Spacer(1, 3.0 * cm)],
            [Paragraph("Final Year Project  ·  2025–2026", styles["cover_tag"])],
            [Spacer(1, 2 * cm)],
        ],
        colWidths=[content_w],
        hAlign="CENTER",
    )
    cover_tbl.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), DEEP_NAVY)]))
    story.append(cover_tbl)
    story.append(PageBreak())

    # ── Section I: Model & Architecture 
    story += section_banner("Part I: Model Selection & Architecture", styles)
    
    story += qa_block(
        1, "Why did you choose a Variational Autoencoder (VAE) instead of a GAN or Transformer for molecular generation?",
        [
            "We chose a VAE primarily because it learns a <b>smooth, continuous, and differentiable latent space</b>. This is the absolute prerequisite for gradient-based property optimization.",
            "• <b>GANs</b> often suffer from mode collapse and do not natively provide an encoder to map discrete molecules into a continuous space. Optimization in GANs requires complex inversions.",
            "• <b>Transformers</b> are state-of-the-art for sequence modeling, but they do not naturally construct a low-dimensional manifold we can traverse. They are autoregressive block-boxes.",
            "• <b>VAEs</b> map molecules to continuous Gaussian distributions. We can encode a seed molecule, follow the gradient of our Property Predictor in this latent space, and decode back to a structurally similar but optimized molecule."
        ], styles
    )

    story += qa_block(
        2, "You mentioned KL Annealing. What is 'posterior collapse', and why does KL annealing fix it?",
        [
            "<b>Posterior collapse</b> is a common failure mode in text/sequence VAEs utilizing powerful autoregressive decoders (like our GRU). The model learns to ignore the latent variable z entirely, relying solely on the GRU's autoregressive power, while the encoder simply outputs the prior N(0,I).",
            "This happens because the KL divergence penalty is too strong early in training. The model 'panics' and forces the posterior to match the prior to reduce KL loss to zero, rather than learning to reconstruct the molecule.",
            "<b>KL Annealing</b> fixes this by introducing a weight (from 0 to 1) applied to the KL term. In our case, we linearly ramp this weight during the first 4 epochs. This allows the model to first minimize reconstruction loss (learning chemical grammar) before we slowly introduce the KL penalty to smooth the latent space."
        ], styles
    )

    story += qa_block(
        3, "Why use an Adaptive Max Pool between your convolution layers and linear heads in the encoder?",
        [
            "We use Adaptive Max Pooling because the sequence lengths can technically vary, and we need a consistent pathway from the variable-length sequence dimension to the fixed-size dense layers representing our latent distributions (μ and log variance).",
            "By pooling across the spatial sequence dimension, we extract the most important learned features (e.g., the presence of a specific functional group or ring motif) regardless of exactly where it appeared in the sequence."
        ], styles
    )

    # ── Section II: Dataset & Processing 
    story += section_banner("Part II: Dataset & Preprocessing", styles)

    story += qa_block(
        4, "Why did you use SELFIES instead of SMILES? Doesn't SELFIES have any drawbacks?",
        [
            "We used SELFIES because it is designed to be <b>100% physically and syntactically robust</b>. If you randomly mutate a SMILES string (e.g., forget a closing parenthesis for a ring), it becomes invalid and RDKit will throw an error. SELFIES strings, by definition, always decode to a chemically valid graph.",
            "This is crucial for our decoder, which samples tokens probabilistically. Using SELFIES drastically improves the baseline validity rate of our generation.",
            "The <b>drawback</b> of SELFIES is that the sequences are typically 20-30% longer than SMILES because of the bracket notation constraints, and the vocabulary is larger, which slightly increases the model parameter count and training time."
        ], styles
    )

    story += qa_block(
        5, "Your Lipinski Rule of Five pass rate is ~28%. If your training data (ZINC250K) was already filtered for Lipinski rules, why is this rate relatively low?",
        [
            "This is expected behavior in continuous latent space generation. While the training data follows Lipinski's rules, the VAE does not have hard-coded chemical constraints—it is learning a continuous probability distribution.",
            "When sampling randomly from the standard normal prior N(0,I) in high dimensions, we are exploring regions of the latent space that interpolate between training points. These interpolated regions occasionally decode into molecules that are syntactically valid (due to SELFIES) and chemically valid (passing RDKit), but fall outside the strict thresholds of the Rule of Five (e.g., slightly too heavy or hydrophobic).",
            "28% is actually standard/strong for an unconstrained VAE at this parameter scale."
        ], styles
    )

    # ── Section III: Training & Evaluation 
    story += section_banner("Part III: Training, Predictor & Optimization", styles)

    story += qa_block(
        6, "Why is QED weighted 10x heavier than SAS and LogP in your predictor loss function?",
        [
            "QED (Quantitative Estimate of Drug-likeness) is a normalized score between 0 and 1, whereas LogP and SAS are not bounded (LogP can span -5 to +5, SAS typically 1-10).",
            "Without the 10x multiplier on QED, the MSE error for QED would be orders of magnitude smaller than the MSE for SAS or LogP.",
            "The optimizer would primarily focus on minimizing the larger error scales (LogP and SAS) and essentially ignore QED. Up-weighting the QED loss ensures gradients flow strongly enough for the predictor to learn to distinguish drug-like molecules from non-drug-like ones."
        ], styles
    )

    story += qa_block(
        7, "How are you calculating SAS (Synthetic Accessibility Score) during training? What are the limits of this?",
        [
            "We are using a proxy for SAS: the sum of the number of rings and number of rotatable bonds. The full Ertl SAS algorithm in RDKit requires deep molecular fragment matching, which is computationally expensive to run on millions of molecules during data generation.",
            "The trade-off is precision. Our proxy captures the fundamental intuition of synthetic complexity (more rings and floppy bonds = harder to make), but lacks the nuanced fragment contributions of true SAS. This is why our MSE for SAS is highest (~1.72) compared to our other predictors."
        ], styles
    )

    story += qa_block(
        8, "Your optimization curve shows convergence in about 50 steps. What exactly are you doing mathematically during these steps?",
        [
            "We are performing <b>Gradient Descent in the Latent Space</b>.",
            "1. We take a seed molecule, pass it through the VAE encoder to get a latent vector z.",
            "2. We freeze the predictor model. We pass z through the predictor to get predicted properties (ŷ).",
            "3. We calculate the Mean Squared Error loss between our predictions (ŷ) and our desired Target values (y_target).",
            "4. We take the gradient of the loss <i>with respect to the input vector z</i>: ∇z(Loss).",
            "5. We update z : z = z - learning_rate * ∇z(Loss).",
            "We repeat this ~50 times until the predicted properties map closely to the target, and then decode the final optimized z back into a SELFIES string."
        ], styles
    )

    story += qa_block(
        9, "What happens if your property predictor outputs a QED of 0.9, but when you decode the molecule and check with RDKit, the actual QED is 0.4?",
        [
            "This is known as <b>Predictor Pathology</b> or an adversarial example. Gradient descent effectively finds 'loopholes' in the predictor's approximation of the RDKit function, traversing into regions of latent space where the predictor is wildly overconfident.",
            "To mitigate this, during optimization, we constrain the L2 distance the latent vector z can travel from the original seed. If z drifts too far into unexplored space, it enters the out-of-distribution domain for the predictor, where accuracy collapses. Our optimization logs show the L2 distance stabilizing nicely around 0.09."
        ], styles
    )

    # ── Section IV: Research Paper
    story += section_banner("Part IV: The Research Paper Component", styles)

    story += qa_block(
        10, "For your research paper, why do you expect Bayesian Optimization (BO) to provide better structural diversity than Gradient Descent?",
        [
            "Gradient descent is purely <b>greedy algorithm (exploitation)</b>. It computes a local gradient and marches straight up the hill to the nearest local property maximum. It will always return roughly the same optimized structure for a given region.",
            "Bayesian Optimization deliberately models <b>uncertainty</b> using surrogate models (like a Gaussian Process). An acquisition function (like Upper Confidence Bound) forces BO to sample regions of the latent space where uncertainty is high.",
            "Because BO explicitly combines exploitation with exploration, it can step out of local minima and discover novel, distinctly different molecular scaffolds that satisfy the target properties, hence increasing structural diversity."
        ], styles
    )

    story += qa_block(
        11, "How does Temperature sampling work in your decoder, and why use 1.5 for random generation but 0.8 for targeted optimization?",
        [
            "Temperature sampling scales the logits (raw outputs) of the decoder before the softmax layer: probabilities = softmax(logits / T).",
            "• At <b>T = 1.0</b>, it's a standard softmax.",
            "• At <b>T = 0.8</b> (Optimization), probabilities become sharper. The model exploits its learned grammar heavily and only samples the vastly most probable next tokens. This improves stability, ensuring our highly-optimized latent vector z doesn't get ruined by a noisy decoding step.",
            "• At <b>T = 1.5</b> (Random Gen), the distribution is flattened. The model is forced to pick less probable tokens, encouraging exploration and novelty, which allows us to discover carbon chains and branching structures we might not otherwise see.",
        ], styles
    )

    # Build PDF
    doc = SimpleDocTemplate(
        OUTPUT_PATH, pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=MARGIN + 0.8 * cm, bottomMargin=MARGIN + 0.5 * cm,
        title="SmartChem ML QnA", author="ML Sub-team"
    )
    doc.build(story, canvasmaker=QnACanvas)
    print(f"\n✅ QnA PDF generated successfully:\n   {OUTPUT_PATH}\n")

if __name__ == "__main__":
    build_story()
