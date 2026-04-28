# SmartChem — AI-Powered De Novo Drug Discovery Platform

> **B.Tech Major Project · Department of Computer Science and Engineering**  
> United College of Engineering & Research, Prayagraj  
> Dr. A.P.J. Abdul Kalam Technical University, Lucknow · May 2026

---

## Overview

SmartChem is an end-to-end generative AI platform for **de novo drug molecule design**. It combines a tripartite Variational Autoencoder (VAE) architecture with a full-stack web application, enabling medicinal chemists to generate, optimise, and evaluate novel drug candidates in a single unified interface.

**Key capabilities:**
- **100% structurally valid** molecule generation via SELFIES grammar
- **Three encoder modalities** — 1D dilated CNN, GINEConv GNN, and a novel Hybrid encoder with learned gated modality fusion
- **Latent-space property optimisation** toward user-specified QED, LogP, and SAS targets
- **Four-stage toxicity filtering** (Lipinski, PAINS, Brenk, NIH)
- **Full-stack web platform** — React + TypeScript frontend, FastAPI backend, MongoDB, 3Dmol.js 3D viewer, Llama-3.3-70B AI assistant

---

## Repository Structure

```
Smart Chem/
│
├── backend/                      # FastAPI server — REST API, auth, job queue
├── frontend/                     # React + TypeScript + Vite web application
│
├── models/                       # PyTorch model definitions
│   ├── cnn_vae.py                #   1D Dilated CNN VAE encoder
│   ├── gnn_vae.py                #   GINEConv GNN VAE encoder
│   └── hybrid_vae.py             #   Gated modality fusion hybrid encoder
│
├── training/                     # All training & evaluation scripts
│   ├── preprocess.py             #   ZINC-250k preprocessing pipeline
│   ├── run_training.py           #   Main unified training entry point
│   ├── run_cnn_training.py       #   CNN VAE training
│   ├── run_gnn_training.py       #   GNN VAE training
│   ├── run_hybrid_training.py    #   Hybrid encoder training
│   ├── run_targeted_eval.py      #   Full targeted generation evaluation
│   ├── train_vae.py              #   VAE training loop
│   ├── train_predictor.py        #   Property predictor training loop
│   ├── losses.py                 #   Loss functions (BCE, KL, free-bits)
│   ├── sascorer.py               #   Synthetic Accessibility Score module
│   ├── fpscores.pkl.gz           #   SA scorer fragment score data
│   └── generate_all_plots.py     #   Regenerates evaluation/plots/ figures
│
├── evaluation/                   # Evaluation artefacts
│   ├── logs/                     #   Structured metrics (CSV / JSON)
│   │   └── raw_run_logs/         #   Raw training run stdout logs
│   ├── plots/                    #   Generated report figures (fig1–fig8)
│   └── verify_latent_space.py    #   Active-unit & latent analysis
│
├── data/                         # Preprocessed dataset tensors (gitignored)
├── checkpoints/                  # Saved model weights (gitignored)
├── images/                       # Report figures referenced by LaTeX
│
├── Docs/                         # All project documentation
│   ├── main.tex                  #   LaTeX project report source
│   ├── Research.tex              #   Research paper (original)
│   ├── Research Paper V2.tex     #   Research paper (revised)
│   ├── Major_Project_Report.pdf  #   Compiled final report
│   ├── Project PPT.pptx          #   Viva presentation
│   ├── Research_IEEE.docx        #   IEEE Word submission
│   ├── Research_Paper.docx       #   Research paper Word format
│   ├── PROJECT_DATA.md           #   Detailed project data reference
│   ├── images/                   #   Copy of all report figures
│   └── project_data/             #   Structured JSON knowledge base
│
├── README.md
├── requirements.txt              # Python dependencies
├── .env                          # Environment variables (gitignored)
└── .gitignore
```

---

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Preprocess data

```bash
python training/preprocess.py
```

### 3. Train models (in order)

```bash
python training/run_cnn_training.py       # Stage 1: CNN VAE
python training/run_gnn_training.py       # Stage 2: GNN VAE
python training/run_hybrid_training.py    # Stage 3: Hybrid encoder (base weights frozen)
```

### 4. Run evaluation

```bash
python training/run_targeted_eval.py      # Generates evaluation/logs/ CSVs
python training/generate_all_plots.py     # Regenerates evaluation/plots/ figures
```

### 5. Launch web platform

```bash
# Backend
cd backend && uvicorn main:app --reload

# Frontend (separate terminal)
cd frontend && npm install && npm run dev
```

---

## Evaluation Results Summary

| Metric | CNN | GNN | **Hybrid** |
|--------|-----|-----|-----------|
| Structural Validity | 100% | 100% | **100%** |
| Novelty | 100% | 100% | **100%** |
| Active Latent Units | 37 / 128 | 70 / 128 | **45 / 128** |
| QED Hit Rate | 6.1% | 7.0% | **7.6%** |
| LogP Hit Rate | 13.0% | 10.3% | **11.1%** |
| All-Filter Pass Rate | 10.0% | 17.6% | **18.2%** |
| Lipinski Compliance | 93.7% | 93.4% | **92.3%** |
| Pairwise Diversity | ~0.895 | ~0.896 | **~0.897** |

Training hardware: NVIDIA GeForce RTX 3050 (4 GB VRAM) · Intel Core i5 (12th Gen) · 16 GB RAM

---

## Architecture

```
SELFIES sequence ──► CNN Encoder ──────────┐
                                           ├──► Gated Fusion ──► z ──► GRU Decoder ──► SELFIES
Molecular graph  ──► GNN Encoder ──────────┘
                                           
z ──► MLP Predictor ──► (QED, LogP, SAS)
z ◄── Adam gradient ascent (targeted generation)
```

**Anti-collapse stack:** Cyclical KL annealing (4 cycles) · Free bits · Word dropout (50%) · Latent injection at every GRU step · Weakened single-layer decoder

---

## Tech Stack

| Layer | Technology |
|-------|------------|
| ML Framework | PyTorch 2.1 + PyTorch Geometric |
| Cheminformatics | RDKit · SELFIES · SAscore |
| Backend | FastAPI · Motor (async MongoDB) · PyJWT |
| Frontend | React 18 · TypeScript · Vite · Tailwind CSS |
| 3D Visualisation | 3Dmol.js (WebGL) |
| AI Assistant | Groq Llama-3.3-70B |
| Database | MongoDB |

---

## Authors

| Name | Roll No. |
|------|----------|
| Nishant Singh Bisht | 2200100100228 |
| Navya Singh | 2200100100219 |
| Nishant Kumar Singh | 2200100100225 |
| Shivam Sahani | 2200100100313 |

**Supervisor:** Mr. Sachin Sonker, Assistant Professor  
**Department:** Computer Science and Engineering  
**Institution:** United College of Engineering & Research, Prayagraj
