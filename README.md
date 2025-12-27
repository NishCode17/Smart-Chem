# SmartChem

## Overview

SmartChem is an application that demonstrates how **generative AI and machine learning** can be used for **molecular generation and optimization** in computational drug discovery. The project explores how ML models can generate novel molecular structures and evaluate them using standard chemical properties.

The goal of SmartChem is to showcase the **end-to-end workflow** of molecule generation, evaluation, and optimization, combining machine learning with cheminformatics tools in a usable application.

---

## Problem Context

In early-stage drug discovery, researchers often need to explore a large chemical space to identify molecules with desirable properties such as drug-likeness and synthetic feasibility. This process is computationally intensive and difficult to do manually.

SmartChem addresses this by:
- generating candidate molecules using a generative ML model,
- evaluating their chemical properties,
- allowing targeted optimization based on desired constraints.

---

## System Architecture (Hybrid)

SmartChem uses a **Hybrid Monolithic Architecture** to support both real-time interaction and heavy background processing:

- **FastAPI Backend (API)**  
  Handles synchronous requests (`/generate`) and manages the job queue. Serves as the entry point for the frontend.

- **ML Worker (Background Processor)**  
  A standalone Python process that consumes long-running optimization tasks (`/jobs`) from MongoDB. It uses a shared ML engine to ensure consistency with the API.
  
- **Shared ML Executor**
  A centralized module (`backend/ml_executor.py`) that contains the core logic for VAE inference, property prediction, and RDKit validation, ensuring zero duplication between the API and Worker.

- **MongoDB**  
  Acts as both the persistence layer and the asynchronous job queue.

---

## Workflow

1. **Synchronous**: User requests random generation -> API runs inference immediately -> Returns results (Best for quick interaction).
2. **Asynchronous**: User requests Lead Optimization -> API creates PENDING job -> Worker claims & processes job -> DB Updated -> User polls for results (Best for heavy computation).

---

## Machine Learning & Chemistry Components

SmartChem integrates machine learning and cheminformatics as follows:

- **Generative Model**: Variational Autoencoder (1D CNN/GRU) mapping SELFIES to latent space.
- **Property Predictor**: MLP predicting QED, LogP, and SAS from latent vectors.
- **Optimization**: Gradient Ascent in latent space to maximize predicted properties.
- **RDKit Integration**: Validates chemical structure and calculates physical properties.

---

## How to Run

### Prerequisites
- Python 3.9+
- MongoDB running locally or via URI
- Node.js (for Frontend)

### 1. Setup Environment
```bash
# Install dependencies
pip install -r requirements.txt
```

### 2. Start the API Server
```bash
uvicorn backend.main:app --reload
```
*The API is now running at `http://localhost:8000`*

### 3. Start the ML Worker (New)
In a new terminal window:
```bash
python -m backend.worker
```
*The worker is now listening for optimization jobs.*

### 4. Start Frontend
```bash
cd frontend
npm run dev
```

---

## Tech Stack

- **Backend**: FastAPI
- **Worker**: Python `asyncio` + MongoDB `find_one_and_update`
- **Database**: MongoDB
- **ML / Chemistry**: PyTorch, RDKit, SELFIES
