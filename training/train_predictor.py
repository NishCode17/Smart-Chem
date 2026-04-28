"""
training/train_predictor.py

Trains the property predictor MLP on latent vectors (z) produced by the VAE.

Expected data
-------------
z : (N, 128)  — latent means from the VAE encoder (train split ONLY)
y : (N, 3)    — [QED, LogP, SAS] molecular properties

Data split
----------
- 80 % of molecules → VAE training + predictor training  (train split)
- 20 % of molecules → held-out evaluation only            (test split)
NEVER mix train and test splits.

Use build_latent_dataset() from train_vae.py to produce (z, y).
"""

import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


# ---------------------------------------------------------------------------
# Main training function
# ---------------------------------------------------------------------------

def train_predictor(predictor, z_dataset, optimizer, device, epochs=12,
                    batch_size=256, save_path="predictor.pt"):
    """
    Train the property predictor on a latent dataset.

    Args:
        predictor  : PropertyPredictor model instance
        z_dataset  : tuple (z_tensor, y_tensor)
                       z_tensor : (N, latent_dim)
                       y_tensor : (N, 3)  — [QED, LogP, SAS]
        optimizer  : torch.optim optimizer
        device     : torch.device
        epochs     : number of training epochs (default 12)
        batch_size : mini-batch size (default 256)
        save_path  : where to save the trained predictor weights

    Returns:
        predictor with trained weights
    """
    z_tensor, y_tensor = z_dataset
    z_tensor = z_tensor.to(device)
    y_tensor = y_tensor.to(device)

    loader = DataLoader(
        TensorDataset(z_tensor, y_tensor),
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
    )

    criterion = nn.MSELoss()

    predictor.to(device)
    predictor.train()

    print(f"\n[train_predictor] {len(z_tensor)} samples  "
          f"| {epochs} epochs on {device}\n")

    for epoch in range(1, epochs + 1):
        epoch_loss = 0.0
        n_batches  = 0

        for z_batch, y_batch in loader:
            optimizer.zero_grad()

            preds = predictor(z_batch)             # (batch, 3)
            loss  = criterion(preds, y_batch)

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            n_batches  += 1

        avg_loss = epoch_loss / n_batches
        print(f"Epoch [{epoch:02d}/{epochs}]  MSE Loss: {avg_loss:.6f}")

    torch.save(predictor.state_dict(), save_path)
    print(f"\n[train_predictor] Predictor saved → {save_path}")
    return predictor


# ---------------------------------------------------------------------------
# Entry-point (CNN VAE → latent dataset → predictor)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import json
    import math
    from torch.utils.data import DataLoader, TensorDataset

    sys.path.insert(0, ROOT)
    from models.vae import VAE
    from models.predictor import PropertyPredictor
    from training.train_vae import train_vae, build_latent_dataset

    PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
    CKPT_DIR      = os.path.join(ROOT, "checkpoints")
    MODE          = "selfies"
    DEVICE        = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    BATCH_SIZE    = 64
    VAE_EPOCHS    = 25
    PRED_EPOCHS   = 12
    LATENT_DIM    = 128
    TRAIN_RATIO   = 0.8

    os.makedirs(CKPT_DIR, exist_ok=True)

    # --- Load data ---
    with open(os.path.join(PROCESSED_DIR, f"{MODE}_vocab.json")) as fh:
        vocab = json.load(fh)
    vocab_size = len(vocab) + 3

    data = torch.load(
        os.path.join(PROCESSED_DIR, f"train_{MODE}.pt"), weights_only=True
    )

    # FIX 2 — 80/20 train/test split (no data leakage)
    n_total = len(data)
    n_train = math.floor(n_total * TRAIN_RATIO)
    train_data = data[:n_train]
    test_data  = data[n_train:]       # held-out — NOT used for VAE or predictor
    print(f"[split] total={n_total}  train={n_train}  test={n_total - n_train}")

    # Load aligned SMILES and apply the SAME split
    smiles_path  = os.path.join(PROCESSED_DIR, "train_smiles.txt")
    train_smiles = None
    # test_smiles kept separate for downstream evaluation
    if os.path.isfile(smiles_path):
        with open(smiles_path) as fh:
            all_smiles = [line.strip() for line in fh if line.strip()]
        train_smiles = all_smiles[:n_train]   # FIX 3 — only train split
        test_smiles  = all_smiles[n_train:]   # kept separate, not used here
        print(f"[split] train_smiles={len(train_smiles)}  "
              f"test_smiles={len(test_smiles)}")

    # Unwrap TensorDataset tuples so train_vae receives plain tensors
    class _UnwrapLoader:
        def __init__(self, dl):
            self._dl = dl
        def __iter__(self):
            for (x,) in self._dl:
                yield x
        def __len__(self):
            return len(self._dl)

    # --- Stage 1: Train VAE on train split only ---
    train_raw_loader = DataLoader(
        TensorDataset(train_data),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=True,
        num_workers=0,
    )
    vae      = VAE(vocab_size=vocab_size).to(DEVICE)
    vae_opt  = optim.Adam(vae.parameters(), lr=1e-3)
    vae_path = os.path.join(CKPT_DIR, "vae_cnn.pt")

    train_vae(vae, _UnwrapLoader(train_raw_loader), vae_opt, DEVICE,
              epochs=VAE_EPOCHS, save_path=vae_path)

    # --- Stage 2: Build latent dataset from train split only ---
    label_loader = DataLoader(
        TensorDataset(train_data),
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=0,
    )

    z_tensor, y_tensor = build_latent_dataset(
        vae,
        _UnwrapLoader(label_loader),
        DEVICE,
        smiles_list=train_smiles,   # FIX 3 — train SMILES only
    )

    # --- Stage 3: Train predictor on train latents only ---
    predictor = PropertyPredictor(latent_dim=LATENT_DIM).to(DEVICE)
    pred_opt  = optim.Adam(predictor.parameters(), lr=1e-3)
    pred_path = os.path.join(CKPT_DIR, "predictor.pt")

    train_predictor(
        predictor,
        (z_tensor, y_tensor),
        pred_opt,
        DEVICE,
        epochs=PRED_EPOCHS,
        save_path=pred_path,
    )

    print("\n[Done] All stages complete.")
    print(f"  VAE checkpoint       → {vae_path}")
    print(f"  Predictor checkpoint → {pred_path}")
    print(f"  Test split size      → {len(test_data)} molecules (held out)")
