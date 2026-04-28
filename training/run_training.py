"""
run_training.py

Full pipeline orchestration for SmartChem molecular VAE system.

Stages
------
1. Train CNN-VAE
2. Train GNN-VAE
3. Train Hybrid-VAE
4. Build latent dataset (from CNN-VAE encoder)
5. Train property predictor

Run from project root:
    python run_training.py
"""

import os
import sys
import json
import math
import torch
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, Subset

# ---------------------------------------------------------------------------
# Project root on sys.path
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
from models.vae            import VAE
from models.gnn_encoder    import GNNEncoder
from models.hybrid_encoder import HybridEncoder
from models.predictor      import PropertyPredictor

from training.train_vae       import train_vae, build_latent_dataset
from training.train_predictor import train_predictor

from torch_geometric.loader import DataLoader as PyGDataLoader
from torch_geometric.data   import Batch as PyGBatch

# ---------------------------------------------------------------------------
# SECTION 1 — CONFIG
# ---------------------------------------------------------------------------
DEVICE           = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE       = 64
VAE_EPOCHS       = 25
PREDICTOR_EPOCHS = 12
LR               = 1e-3
TRAIN_RATIO      = 0.8

NODE_DIM   = 17
EDGE_DIM   = 6
LATENT_DIM = 128

PROCESSED_DIR    = os.path.join(ROOT, "data", "processed")
RAW_DATA_PATH    = os.path.join(ROOT, "data", "raw", "train_molecules.csv")
CKPT_DIR         = os.path.join(ROOT, "checkpoints")
TARGET_MOLECULES = 100_000

os.makedirs(CKPT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# SECTION 2 — DATA  (100 k molecules via preprocess_zinc pipeline)
# ---------------------------------------------------------------------------
from data.preprocess_zinc import build_full_dataset

splits = build_full_dataset(
    input_path  = RAW_DATA_PATH,
    output_dir  = PROCESSED_DIR,
    target_size = TARGET_MOLECULES,
    max_atoms   = 60,
    max_len     = 100,
    train_ratio = TRAIN_RATIO,
)

train_seq_data      = splits["train_seq"]
test_seq_data       = splits["test_seq"]
train_graph_dataset = splits["train_graph"]
train_hybrid_ds     = splits["train_hybrid"]
train_smiles        = splits["train_smiles"]
test_smiles         = splits["test_smiles"]
vocab               = splits["vocab"]

VOCAB_SIZE = len(vocab)   # already includes <pad>/<sos>/<eos>
n_train    = len(train_seq_data)
n_test     = len(test_seq_data)
N          = n_train + n_test

print(f"  Total molecules : {N:,}")
print(f"  Train           : {n_train:,}")
print(f"  Test (held-out) : {n_test:,}")
print(f"  Vocab size      : {VOCAB_SIZE}")
print(f"  Device          : {DEVICE}")

# ---------------------------------------------------------------------------
# DataLoaders
# ---------------------------------------------------------------------------

class _UnwrapLoader:
    """Strips TensorDataset tuple wrapper → plain tensor per batch."""
    def __init__(self, dl):
        self._dl = dl
    def __iter__(self):
        for (x,) in self._dl:
            yield x
    def __len__(self):
        return len(self._dl)

cnn_loader = _UnwrapLoader(
    DataLoader(TensorDataset(train_seq_data), batch_size=BATCH_SIZE,
               shuffle=True, drop_last=True, num_workers=0)
)

gnn_loader = PyGDataLoader(
    train_graph_dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True
)

def hybrid_collate(batch):
    seqs   = torch.stack([item[0] for item in batch])
    graphs = PyGBatch.from_data_list([item[1] for item in batch])
    return seqs, graphs

hybrid_loader = DataLoader(
    train_hybrid_ds, batch_size=BATCH_SIZE, shuffle=True,
    drop_last=True, num_workers=0, collate_fn=hybrid_collate
)

cnn_label_loader = _UnwrapLoader(
    DataLoader(TensorDataset(train_seq_data), batch_size=BATCH_SIZE,
               shuffle=False, drop_last=False, num_workers=0)
)


# ---------------------------------------------------------------------------
# SECTION 3 — MODEL INITIALISATION
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print("  Initialising models")
print(f"{'='*60}")

# --- CNN VAE (default encoder) ---
vae_cnn = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = 100,
).to(DEVICE)

# --- GNN VAE ---
vae_gnn = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = 100,
).to(DEVICE)
vae_gnn.encoder = GNNEncoder(
    node_dim   = NODE_DIM,
    edge_dim   = EDGE_DIM,
    hidden_dim = 128,
    latent_dim = LATENT_DIM,
).to(DEVICE)

# --- Hybrid VAE ---
vae_hybrid = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = 100,
).to(DEVICE)
_gnn_for_hybrid = GNNEncoder(
    node_dim   = NODE_DIM,
    edge_dim   = EDGE_DIM,
    hidden_dim = 128,
    latent_dim = LATENT_DIM,
).to(DEVICE)
vae_hybrid.encoder = HybridEncoder(
    cnn_encoder = vae_hybrid.encoder,   # original CNNEncoder already on device
    gnn_encoder = _gnn_for_hybrid,
    latent_dim  = LATENT_DIM,
).to(DEVICE)

print(f"  CNN VAE    : {sum(p.numel() for p in vae_cnn.parameters()):,} params")
print(f"  GNN VAE    : {sum(p.numel() for p in vae_gnn.parameters()):,} params")
print(f"  Hybrid VAE : {sum(p.numel() for p in vae_hybrid.parameters()):,} params")

# ---------------------------------------------------------------------------
# SECTION 4 — TRAIN ALL THREE VAEs
# ---------------------------------------------------------------------------

# ---- Stage 1: CNN VAE ----
print(f"\n{'='*60}")
print("  STAGE 1 — CNN VAE Training")
print(f"{'='*60}")

opt_cnn = optim.Adam(vae_cnn.parameters(), lr=LR)
train_vae(
    model      = vae_cnn,
    dataloader = cnn_loader,
    optimizer  = opt_cnn,
    device     = DEVICE,
    epochs     = VAE_EPOCHS,
    save_path  = os.path.join(CKPT_DIR, "vae_cnn.pt"),
)
print("  ✅ CNN VAE training complete")

# ---- Stage 2: GNN VAE ----
print(f"\n{'='*60}")
print("  STAGE 2 — GNN VAE Training")
print(f"{'='*60}")

opt_gnn = optim.Adam(vae_gnn.parameters(), lr=LR)
train_vae(
    model      = vae_gnn,
    dataloader = gnn_loader,
    optimizer  = opt_gnn,
    device     = DEVICE,
    epochs     = VAE_EPOCHS,
    save_path  = os.path.join(CKPT_DIR, "vae_gnn.pt"),
)
print("  ✅ GNN VAE training complete")

# ---- Stage 3: Hybrid VAE ----
print(f"\n{'='*60}")
print("  STAGE 3 — Hybrid VAE Training")
print(f"{'='*60}")

opt_hybrid = optim.Adam(vae_hybrid.parameters(), lr=LR)
train_vae(
    model      = vae_hybrid,
    dataloader = hybrid_loader,
    optimizer  = opt_hybrid,
    device     = DEVICE,
    epochs     = VAE_EPOCHS,
    save_path  = os.path.join(CKPT_DIR, "vae_hybrid.pt"),
)
print("  ✅ Hybrid VAE training complete")

# ---------------------------------------------------------------------------
# SECTION 5 — BUILD LATENT DATASET (CNN-VAE encoder, train split only)
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print("  STAGE 4 — Building latent dataset (CNN-VAE)")
print(f"{'='*60}")

z_tensor, y_tensor = build_latent_dataset(
    model       = vae_cnn,
    dataloader  = cnn_label_loader,
    device      = DEVICE,
    smiles_list = train_smiles,     # train split only — no leakage
)

# ---------------------------------------------------------------------------
# SECTION 6 — TRAIN PREDICTOR
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print("  STAGE 5 — Predictor Training")
print(f"{'='*60}")

predictor = PropertyPredictor(latent_dim=LATENT_DIM).to(DEVICE)
opt_pred  = optim.Adam(predictor.parameters(), lr=LR)

train_predictor(
    predictor  = predictor,
    z_dataset  = (z_tensor, y_tensor),
    optimizer  = opt_pred,
    device     = DEVICE,
    epochs     = PREDICTOR_EPOCHS,
    save_path  = os.path.join(CKPT_DIR, "predictor.pt"),
)
print("  ✅ Predictor training complete")

# ---------------------------------------------------------------------------
# SECTION 7 — SAVE ALL MODELS
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print("  Saving model weights")
print(f"{'='*60}")

torch.save(vae_cnn.state_dict(),    os.path.join(CKPT_DIR, "vae_cnn.pt"))
torch.save(vae_gnn.state_dict(),    os.path.join(CKPT_DIR, "vae_gnn.pt"))
torch.save(vae_hybrid.state_dict(), os.path.join(CKPT_DIR, "vae_hybrid.pt"))
torch.save(predictor.state_dict(),  os.path.join(CKPT_DIR, "predictor.pt"))

print(f"  vae_cnn.pt    → {CKPT_DIR}")
print(f"  vae_gnn.pt    → {CKPT_DIR}")
print(f"  vae_hybrid.pt → {CKPT_DIR}")
print(f"  predictor.pt  → {CKPT_DIR}")

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
print(f"\n{'#'*60}")
print("  SmartChem — Training Complete")
print(f"  Device          : {DEVICE}")
print(f"  Train molecules : {n_train}")
print(f"  Test molecules  : {n_test}  (held-out, not used)")
print(f"  Latent samples  : {len(z_tensor)}")
print(f"{'#'*60}\n")
