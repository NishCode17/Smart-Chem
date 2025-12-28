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

## System Architecture (Event-Driven Async)

SmartChem uses a **Fully Asynchronous Event-Driven Architecture** to ensure high scalability and responsiveness. All heavy Machine Learning computations are decoupled from the user-facing API.

- **FastAPI Gateway (API)**  
  Acts as a lightweight "Receptionist". It authenticates users, validates inputs, and submits tasks to the job queue (`/jobs`). It **never** runs ML inference directly, ensuring the server stays responsive.

- **ML Worker (Background Processor)**  
  A dedicated process that consumes tasks (Targeted Gen, Lead Opt) from MongoDB. It claims jobs, runs the heavy VAE/RDKit math, and updates the database with results.
  
- **Shared ML Executor**
  The "Brain" of the operation. A centralized module loaded by the Worker to perform the actual math.

- **MongoDB (Queue & State)**  
  Acts as the message broker. Jobs transition atomically from `PENDING` $\to$ `PROCESSING` $\to$ `COMPLETED`.

---

## Workflow

1. **Submission**: User sends request (e.g., "Optimize this Lead") -> API creates a Job Ticket -> Returns `job_id` ($< 100ms$).
2. **Queueing**: The Job sits in MongoDB as `PENDING`.
3. **Processing**: The Worker (looping in background) claims the job, locks it, and runs the ML algorithms (`~2s - 20s`).
4. **Completion**: Worker saves results to DB.
5. **Pollling**: Frontend automatically checks job status and displays results once `COMPLETED`.

---

## Machine Learning & Chemistry Components

## Machine Learning & Chemistry Components

SmartChem integrates machine learning and cheminformatics as follows:

### 1. Generative Models
- **Architecture**: Variational Autoencoder (VAE) trained on SELFIES representation.
- **Latent Space**: 128-dimensional continuous space where similar molecules map to similar vectors.
- **Property Predictor**: A feed-forward neural network (MLP) acts as a surrogate model to predict QED, LogP, and SAS directly from latent vectors ($z$).

### 2. Optimization Algorithms
- **Targeted Generation (Gradient Ascent)**:  
  Starts from random noise ($z$) and iteratively updates the vector using gradients from the Property Predictor to maximize the match with user-defined targets (e.g., QED=0.9).
  
- **Lead Optimization (Neighborhood Search)**:  
  Takes an existing molecule, encodes it to ($z_{lead}$), and samples a "cloud" of 200+ perturbations (Gaussian noise). These are decoded and filtered to find structural analogs with improved properties.

### 3. Validation & Cheminformatics
- **RDKit Integration**: Acts as the ground-truth validator.
  - Filters invalid SMILES strings.
  - Rejects molecules with only Carbon atoms (sanity check).
  - Calculates accurate physicochemical properties for the final results.

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
