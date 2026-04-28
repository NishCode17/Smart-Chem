import os
import sys
import torch
import torch.optim as optim
from torch_geometric.loader import DataLoader

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from models.vae         import VAE
from models.gnn_encoder import GNNEncoder
from training.train_vae import train_vae

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
DEVICE       = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE   = 128          # raised from 64 → 128 for GNN (graph batching heavier)
EPOCHS       = 30           # raised to fit 4 cyclical cycles
LR           = 3e-4         # AdamW standard
WEIGHT_DECAY = 0.01
TRAIN_RATIO  = 0.8

NODE_DIM   = 17
EDGE_DIM   = 6
LATENT_DIM = 128
MAX_LEN    = 100

PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
RAW_DATA_PATH = os.path.join(ROOT, "data", "raw", "train_molecules.csv")
CKPT_DIR      = os.path.join(ROOT, "checkpoints")
os.makedirs(CKPT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# LOAD EXISTING DATASET
# ---------------------------------------------------------------------------
from data.preprocess_zinc import build_full_dataset

print("\n" + "="*60)
print("  Loading Dataset")
print("="*60)

splits = build_full_dataset(
    input_path  = RAW_DATA_PATH,
    output_dir  = PROCESSED_DIR,
    target_size = 100_000,
    max_atoms   = 60,
    max_len     = MAX_LEN,
    train_ratio = TRAIN_RATIO,
)

train_graph_dataset = splits["train_graph"]
train_smiles        = splits["train_smiles"]
vocab               = splits["vocab"]
VOCAB_SIZE          = len(vocab)

print(f"  Vocab size      : {VOCAB_SIZE}")
print(f"  Train graph size: {len(train_graph_dataset):,}")
print(f"  Device          : {DEVICE}")

# ---------------------------------------------------------------------------
# DATALOADER
# ---------------------------------------------------------------------------
train_loader = DataLoader(
    train_graph_dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True
)

# ---------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("  Initializing GNN VAE Model")
print("="*60)

vae_gnn = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = MAX_LEN,
).to(DEVICE)

# Replace the default (CNN) encoder with the robust GNN encoder
vae_gnn.encoder = GNNEncoder(
    node_dim   = NODE_DIM,
    edge_dim   = EDGE_DIM,
    hidden_dim = 128,
    latent_dim = LATENT_DIM,
).to(DEVICE)

print(f"  GNN VAE Parameters: {sum(p.numel() for p in vae_gnn.parameters()):,}")

# ---------------------------------------------------------------------------
# TRAINING
# ---------------------------------------------------------------------------
optimizer = optim.AdamW(vae_gnn.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

print("\n" + "="*60)
print("  Training GNN VAE")
print("="*60)

vae_gnn = train_vae(
    model            = vae_gnn,
    dataloader       = train_loader,
    optimizer        = optimizer,
    device           = DEVICE,
    epochs           = EPOCHS,
    word_drop        = 0.5,
    n_cycles         = 4,
    max_beta         = 1.0,
    lambda_recon     = 0.5,
    free_bits        = 0.5,
    lr_warmup_epochs = 2,
    base_lr          = LR,
    save_path        = os.path.join(CKPT_DIR, "vae_gnn.pt"),
)

# Save one explicitly just to be safe
torch.save(vae_gnn.state_dict(), os.path.join(CKPT_DIR, "vae_gnn.pt"))

print("\nGNN VAE training complete")
