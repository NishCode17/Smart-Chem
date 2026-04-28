"""
run_cnn_training.py

Standalone CNN-VAE training script for 100k molecules.
Loads cached artefacts from data/processed/ (no preprocessing).
Saves checkpoint to checkpoints/vae_cnn.pt

Key hyperparameter changes vs old 50k config:
  - AdamW (weight_decay=0.01) instead of Adam
  - Batch size 256 (was 64) — fewer gradient steps/epoch, same as 50k×64 dynamics
  - lambda_recon=0.5 (was 0.01) — reconstruction must compete with KL
  - Cyclical KL annealing with 4 cycles (was linear)
  - free_bits=0.5 nats/dim (was 0.1) — per-dim floor prevents selective collapse
  - word_drop=0.5 (was 0.3) — weaker decoder forces z usage
  - 30 epochs (was 25) — more cycles to stabilise

Run from project root:
    python run_cnn_training.py
"""

import os
import sys
import json
import torch
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from models.vae         import VAE
from training.train_vae import train_vae

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
DEVICE       = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE   = 256          # raised from 64 → 256 (fewer steps/epoch ≈ 50k dynamics)
EPOCHS       = 30           # raised from 25 → 30 (more cyclical annealing cycles)
LR           = 3e-4         # lowered from 1e-3 → 3e-4 (AdamW standard for transformers/VAEs)
WEIGHT_DECAY = 0.01         # AdamW weight decay
LATENT_DIM   = 128
MAX_LEN      = 100
TRAIN_RATIO  = 0.8

PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
CKPT_DIR      = os.path.join(ROOT, "checkpoints")
os.makedirs(CKPT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# LOAD CACHED DATA
# ---------------------------------------------------------------------------
pt_path     = os.path.join(PROCESSED_DIR, "train_selfies.pt")
vocab_path  = os.path.join(PROCESSED_DIR, "selfies_vocab.json")

for p in [pt_path, vocab_path]:
    if not os.path.isfile(p):
        raise FileNotFoundError(
            f"Missing cached artefact: {p}\n"
            "Run  python data/preprocess_zinc.py  first to build the dataset."
        )

seq_tensor = torch.load(pt_path, weights_only=True)
with open(vocab_path) as fh:
    vocab = json.load(fh)

VOCAB_SIZE = len(vocab)
N          = len(seq_tensor)
n_train    = int(N * TRAIN_RATIO)

train_seq = seq_tensor[:n_train]

print(f"\n{'='*60}")
print("  CNN VAE — Standalone Training (100k, Anti-Collapse Config)")
print(f"{'='*60}")
print(f"  Total molecules  : {N:,}")
print(f"  Train split      : {n_train:,}")
print(f"  Vocab size       : {VOCAB_SIZE}")
print(f"  Batch size       : {BATCH_SIZE}  (raised from 64 to reduce steps/epoch)")
print(f"  Batches/epoch    : {n_train // BATCH_SIZE:,}")
print(f"  Device           : {DEVICE}")
print(f"{'='*60}\n")

# ---------------------------------------------------------------------------
# DATALOADER
# ---------------------------------------------------------------------------
class _UnwrapLoader:
    def __init__(self, dl): self._dl = dl
    def __iter__(self):
        for (x,) in self._dl: yield x
    def __len__(self): return len(self._dl)

cnn_loader = _UnwrapLoader(
    DataLoader(TensorDataset(train_seq), batch_size=BATCH_SIZE,
               shuffle=True, drop_last=True, num_workers=0,
               pin_memory=(DEVICE.type == "cuda"))
)

# ---------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------
vae_cnn = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = MAX_LEN,
).to(DEVICE)

print(f"  Parameters : {sum(p.numel() for p in vae_cnn.parameters()):,}\n")

# AdamW: decoupled weight decay prevents L2 penalty on Adam's adaptive LR
# (Loshchilov & Hutter 2019) — better regularisation than plain Adam
optimizer = optim.AdamW(vae_cnn.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

# ---------------------------------------------------------------------------
# TRAIN  (anti-collapse configuration)
# ---------------------------------------------------------------------------
train_vae(
    model            = vae_cnn,
    dataloader       = cnn_loader,
    optimizer        = optimizer,
    device           = DEVICE,
    epochs           = EPOCHS,
    word_drop        = 0.5,      # stronger decoder weakening at 100k scale
    n_cycles         = 4,        # cyclical KL: 4 cycles over 30 epochs
    max_beta         = 1.0,      # KL ceiling
    lambda_recon     = 0.5,      # CRITICAL: was 0.01 — reconstruction must dominate
    free_bits        = 0.5,      # CRITICAL: per-dim 0.5 nats (was global 0.1)
    lr_warmup_epochs = 2,        # linear LR warmup
    base_lr          = LR,
    save_path        = os.path.join(CKPT_DIR, "vae_cnn.pt"),
)

print(f"\n{'#'*60}")
print("  CNN VAE Training Complete")
print(f"  Checkpoint → {os.path.join(CKPT_DIR, 'vae_cnn.pt')}")
print(f"{'#'*60}\n")
