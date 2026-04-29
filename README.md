# SmartChem — AI-Powered De Novo Drug Discovery Platform

> **B.Tech Major Project · Department of Computer Science and Engineering**  
> United College of Engineering & Research, Prayagraj  
> Dr. A.P.J. Abdul Kalam Technical University, Lucknow · May 2026

---

## Overview

SmartChem is an end-to-end generative AI platform for **de novo drug molecule design**. It combines a tripartite Variational Autoencoder (VAE) architecture with a full-stack web application, enabling medicinal chemists to generate, optimise, and evaluate novel drug candidates in a single unified interface.

**Key capabilities:**
- **100% structurally valid** molecule generation via SELFIES grammar
- **Three encoder modalities** — 1D dilated CNN, GINEConv GNN, and a Hybrid encoder with learned gated modality fusion
- **Latent-space property optimisation** toward user-specified QED, LogP, and SAS targets
- **Four-stage toxicity filtering** (Lipinski, PAINS, Brenk, NIH)
- **Comprehensive ADMET profiling** — 20+ endpoints (HIA, Caco-2, BBB, CYP450, hERG, AMES, DILI, etc.)
- **Full-stack web platform** — React + TypeScript frontend, FastAPI backend, MongoDB Atlas, 3Dmol.js 3D viewer, Llama-3.3-70B AI assistant

---

## Repository Structure

```
Smart Chem/
│
├── backend/                      # FastAPI server — REST API, auth, job queue
│   ├── main.py                   #   App entry point + all route definitions
│   ├── admet.py                  #   Comprehensive ADMET predictor (20+ endpoints)
│   ├── chem_utils.py             #   RDKit property calculation & filtering
│   ├── ml_executor.py            #   VAE generation & optimization logic
│   ├── optimizer.py              #   Latent space gradient optimization
│   ├── assistant.py              #   Groq LLM chat integration
│   ├── database.py               #   MongoDB Atlas connection
│   ├── models.py                 #   Pydantic data models
│   ├── auth_utils.py             #   JWT auth helpers
│   └── routers/                  #   auth, projects, molecules, jobs
│
├── frontend/                     # React + TypeScript + Vite web application
│
├── models/                       # PyTorch model definitions
│   ├── vae.py                    #   VAE (CNN / GNN / Hybrid encoder)
│   └── predictor.py              #   Property predictor (QED, LogP, SAS)
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
│   │   └── raw_run_logs/         #   Raw training stdout logs
│   └── plots/                    #   Generated report figures (fig1–fig8)
│
├── data/                         # Preprocessed dataset tensors (gitignored)
├── checkpoints/                  # Saved model weights (gitignored)
├── images/                       # Report figures referenced by LaTeX
│
├── Docs/                         # All project documentation
│   ├── main.tex                  #   LaTeX project report source
│   ├── Major_Project_Report.pdf  #   Compiled final report
│   ├── Project PPT.pptx          #   Viva presentation
│   └── ...
│
├── README.md
├── requirements.txt
├── .env                          # Environment variables (gitignored)
└── .gitignore
```

---

## Prerequisites

- **Conda environment:** `smartchem` (Python 3.10+)
- **Node.js** 18+ (for frontend)
- **MongoDB Atlas** cluster with credentials in `.env`

Install Python dependencies (first time only):
```bash
conda activate smartchem
pip install -r requirements.txt
```

---

## Running the Project

> **All commands are run from the project root:**  
> `d:\Final Year Project\Smart Chem\`

### Terminal 1 — Backend (FastAPI)

```bash
conda activate smartchem
python -m uvicorn backend.main:app --reload --port 8000
```

The API will be live at **http://localhost:8000**  
Interactive docs at **http://localhost:8000/docs**

### Terminal 2 — Frontend (React + Vite)

```bash
cd frontend
npm install        # first time only
npm run dev
```

The web app will be live at **http://localhost:8080**

---

> **Important:** Always run the backend with `python -m uvicorn backend.main:app` from the **project root**, not from inside the `backend/` subfolder. Running from inside `backend/` breaks relative import paths (`models/`, `checkpoints/`, `data/`).

---

## Training (optional — checkpoints already included)

```bash
conda activate smartchem

# Stage 1: CNN VAE
python training/run_cnn_training.py

# Stage 2: GNN VAE
python training/run_gnn_training.py

# Stage 3: Hybrid encoder (requires Stage 1 & 2 checkpoints)
python training/run_hybrid_training.py
```

## Evaluation

```bash
python training/run_targeted_eval.py   # generates evaluation/logs/ CSVs
python training/generate_all_plots.py  # regenerates evaluation/plots/ figures
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

Training hardware: NVIDIA GeForce RTX 3050 (4 GB VRAM) · Intel Core i5-12th Gen · 16 GB RAM

---

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/auth/register` | Register new user |
| POST | `/auth/login` | Login (returns JWT) |
| POST | `/generate` | Random molecule generation |
| POST | `/generate/targeted` | Targeted generation (QED/LogP/SAS) |
| POST | `/optimize/lead` | Lead compound optimisation |
| POST | `/utils/analyze` | Full property analysis for a SMILES |
| POST | `/utils/admet` | Comprehensive 20+ endpoint ADMET profile |
| POST | `/utils/3d` | Generate 3D MOL block |
| POST | `/assistant/chat` | LLM chemistry assistant |
| GET | `/projects/` | List user projects |
| GET | `/molecules/` | List saved molecules |

---

## Architecture

```
SELFIES sequence --> CNN Encoder ──────────┐
                                           ├──> Gated Fusion --> z --> GRU Decoder --> SELFIES
Molecular graph  --> GNN Encoder ──────────┘

z --> MLP Predictor --> (QED, LogP, SAS)
z <-- Adam gradient ascent (targeted generation)
```

**Anti-collapse stack:** Cyclical KL annealing (4 cycles) · Free bits · Word dropout (50%) · Latent injection at every GRU step · Weakened single-layer decoder

---

## Tech Stack

| Layer | Technology |
|-------|------------|
| ML Framework | PyTorch 2.1 + PyTorch Geometric |
| Cheminformatics | RDKit · SELFIES · SAscore |
| Backend | FastAPI · Motor (async MongoDB) · PyJWT |
| Frontend | React 18 · TypeScript · Vite |
| 3D Visualisation | 3Dmol.js (WebGL) |
| AI Assistant | Groq Llama-3.3-70B |
| Database | MongoDB Atlas |

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
