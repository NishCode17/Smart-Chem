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

## System Architecture (High-Level)

SmartChem is implemented as a web-based application with a backend API that coordinates ML execution:

- **FastAPI Backend**  
  Provides REST APIs for molecule generation and optimization requests.

- **Asynchronous Processing**  
  ML tasks can take time, so requests are handled asynchronously to keep the application responsive.

- **MongoDB**  
  Used to store request metadata, job status, and generated molecular results.

---

## Workflow

1. A user submits a request to generate or optimize molecules.
2. The backend schedules the ML task and immediately returns a request identifier.
3. The ML model performs molecule generation or optimization.
4. Generated molecules are evaluated using chemical property metrics.
5. Results are stored and made available for retrieval.

---

## Machine Learning & Chemistry Components

SmartChem integrates machine learning and cheminformatics as follows:

- **Generative Model**  
  Uses a generative ML model (based on a Variational Autoencoder) to produce novel molecular structures.

- **Molecular Representation**  
  Molecules are processed in a structured representation suitable for ML-based generation.

- **Property Evaluation**  
  Each generated molecule is evaluated using standard metrics:
  - QED (Quantitative Estimation of Drug-likeness)
  - LogP (lipophilicity)
  - SAS (Synthetic Accessibility Score)

- **RDKit Integration**  
  RDKit is used for molecule validation, property calculation, and filtering of invalid structures.

The emphasis of the project is on **applying generative models to a real problem**, rather than on proposing new ML architectures.

---

## Key Learnings

- Applying generative ML models to chemical structure generation
- Integrating ML inference into an application workflow
- Handling long-running ML tasks in an application setting
- Using cheminformatics tools (RDKit) alongside ML models

---

## Tech Stack

- **Backend**: FastAPI
- **Database**: MongoDB
- **ML / Chemistry**: VAE, PyTorch, RDKit
- **APIs**: REST
