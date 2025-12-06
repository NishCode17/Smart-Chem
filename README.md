# Smart Chem 🧪

Smart Chem is a deep learning-based project designed for *de novo* molecule generation. It leverages Variational Autoencoders (VAEs) to learn the continuous latent space of chemical structures, enabling the improved generation and optimization of drug-like molecules.

## 🚀 Key Features

*   **Dual Representation**: Supports both **SMILES** and **SELFIES** (robust molecular representation) for training.
*   **Deep Learning Model**: Implements a Variational Autoencoder (VAE) to encode and decode molecular structures.
*   **Chemical Analysis**: Backend utilities to calculate key drug-discovery metrics:
    *   **QED** (Quantitative Estimation of Drug-likeness)
    *   **LogP** (Octanol-water partition coefficient)
*   **Visualization**: Automatic generation of 2D molecular structure images encoded in Base64 for easy frontend integration.

## 📂 Project Structure

```
Smart Chem/
├── backend/            # Backend logic for chemical property calculation & image generation
│   └── chem_utils.py   # Utilities using RDKit and Selfies
├── checkpoints/        # Saved model weights during training
├── data/               # Dataset storage
│   ├── raw/            # Original CSV datasets (e.g., train_molecules.csv)
│   └── processed/      # Tokenized tensors and vocabulary files (.pt, .json)
├── models/             # PyTorch model definitions (VAE)
├── frontend/           # Frontend application files
├── preprocess.py       # Script to digest raw CSVs into training tensors
└── train.py            # Main training loop for the VAE model
```

## 🛠️ Setup & Usage

### 1. Prerequisites
Ensure you have the following libraries installed:
*   Python 3.8+
*   PyTorch
*   RDKit
*   Selfies
*   Pandas
*   Tqdm

### 2. Data Preprocessing
Before training, the raw molecular data (SMILES/SELFIES) must be tokenized and converted into tensors.
This step:
*   Filters molecules by length.
*   Builds a vocabulary mapping.
*   Saves processed data to `data/processed/`.

```bash
python preprocess.py
```
*Note: Ensure your raw data is placed at `data/raw/train_molecules.csv`.*

### 3. Model Training
Train the VAE model on the processed data. You can configure hyperparameters (Epochs, Batch Size, Learning Rate) directly in `train.py`.

```bash
python train.py
```
*   **Checkpoints**: Model weights are saved to `checkpoints/` after every few epochs.
*   **Loss Function**: Uses a combination of Reconstruction Loss (CrossEntropy) and KL Divergence (with annealing).

## 📊 Data Pipeline

1.  **Input**: Raw CSV with `smiles` and `SELFIES` columns.
2.  **Processing**: Molecules are tokenized, padded to a fixed length, and mapped to integers.
3.  **Training**: The VAE learns to map these integer sequences to a latent space and back.
4.  **Inference**: The trained decoder serves as a generator for new, valid molecules.

---
*Final Year Project*
