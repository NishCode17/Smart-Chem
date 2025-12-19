# SmartChem

## Introduction

SmartChem is a backend-orchestrated machine learning system designed to address the challenge of molecule generation and property optimization. It demonstrates the integration of generative models within a web service architecture, focusing on the engineering of long-running computational tasks.

This project functions as a system for chemical space exploration, decoupling complex ML workloads from synchronous API interaction. It is designed to illustrate backend orchestration patterns rather than serve as a consumer product.

## Architecture

The system is architected around a **FastAPI** service that acts as the central entry point.

*   **Entry Point**: RESTful API endpoints handle request validation and job dispatch.
*   **Asynchronous Execution**: Long-running ML tasks are executed asynchronously to ensure the API remains responsive.
*   **State Management**: **MongoDB** is used to persist job states, intermediate progress, and final molecular results.

## Asynchronous Execution

To manage the computational latency inherent in generative models, the system implements a controlled asynchronous workflow:

1.  **Job Submission**: Client requests for molecule generation are accepted immediately, returning a unique Job ID.
2.  **Background Execution**: The specific ML task (e.g., generation or optimization) creates a background workload that runs independently of the HTTP response cycle.
3.  **Result Persistence**: Upon completion, results are serialized and stored in the database, updating the job status for subsequent retrieval.

## ML Components

The machine learning logic is scoped to specific cheminformatics tasks:

*   **VAE-based Generation**: Utilizes Variational Autoencoders to map molecular structures to and from a latent space.
*   **Property Evaluation**: Computes standard metrics such as QED (Quantitative Estimation of Drug-likeness), LogP, and SAS (Synthetic Accessibility Score).
*   **Cheminformatics Engine**: Integrates **RDKit** for molecular validation, canonicalization, and rule enforcement.

## Explicit Non-Goals

To clarify the scope and intent of this project:

*   **Not a production drug-discovery platform**: The focus is on software architecture, not novel chemical discovery.
*   **Not a distributed training system**: The system uses pre-trained models for inference and optimization.
*   **Not focused on scaling or deployment**: Scaling strategies and cloud deployment pipelines are outside the scope of this implementation.
