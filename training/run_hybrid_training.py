import os
import sys
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from torch_geometric.data import Batch as PyGBatch

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from models.vae import VAE
from models.cnn_encoder import CNNEncoder
from models.gnn_encoder import GNNEncoder
from models.hybrid_encoder import HybridEncoder
from training.train_vae import train_vae
from data.preprocess_zinc import build_full_dataset

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
DEVICE       = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE   = 128
EPOCHS       = 30
LR           = 3e-4
WEIGHT_DECAY = 0.01
TRAIN_RATIO  = 0.8

NODE_DIM   = 17
EDGE_DIM   = 6
LATENT_DIM = 128
MAX_LEN    = 100

PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
RAW_DATA_PATH = os.path.join(ROOT, "data", "raw", "train_molecules.csv")
CKPT_DIR      = os.path.join(ROOT, "checkpoints")

# ---------------------------------------------------------------------------
# DATA LOAD
# ---------------------------------------------------------------------------
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

train_hybrid_ds = splits["train_hybrid"]
VOCAB_SIZE      = len(splits["vocab"])

print(f"  Vocab size       : {VOCAB_SIZE}")
print(f"  Train hybrid size: {len(train_hybrid_ds):,}")
print(f"  Device           : {DEVICE}")

def hybrid_collate(batch):
    seqs   = torch.stack([item[0] for item in batch])
    graphs = PyGBatch.from_data_list([item[1] for item in batch])
    return seqs, graphs

train_loader = DataLoader(
    train_hybrid_ds, batch_size=BATCH_SIZE, shuffle=True,
    drop_last=True, num_workers=0, collate_fn=hybrid_collate
)

# ---------------------------------------------------------------------------
# LOAD PRETRAINED ENCODERS
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("  Loading Pretrained Encoders")
print("="*60)

# Temp VAE to load CNN weights
cnn_vae = VAE(vocab_size=VOCAB_SIZE, embedding_dim=128, hidden_dim=256, latent_dim=LATENT_DIM, max_len=MAX_LEN).to(DEVICE)
cnn_ckpt_path = os.path.join(CKPT_DIR, "vae_cnn.pt")
cnn_vae.load_state_dict(torch.load(cnn_ckpt_path, map_location=DEVICE))
pretrained_cnn_encoder = cnn_vae.encoder

# Temp VAE to load GNN weights
gnn_vae = VAE(vocab_size=VOCAB_SIZE, embedding_dim=128, hidden_dim=256, latent_dim=LATENT_DIM, max_len=MAX_LEN).to(DEVICE)
gnn_vae.encoder = GNNEncoder(node_dim=NODE_DIM, edge_dim=EDGE_DIM, hidden_dim=128, latent_dim=LATENT_DIM).to(DEVICE)
gnn_ckpt_path = os.path.join(CKPT_DIR, "vae_gnn.pt")
gnn_vae.load_state_dict(torch.load(gnn_ckpt_path, map_location=DEVICE))
pretrained_gnn_encoder = gnn_vae.encoder

print("  Pretrained encoders successfully loaded.")

# ---------------------------------------------------------------------------
# MODEL Setup
# ---------------------------------------------------------------------------
vae_hybrid = VAE(
    vocab_size    = VOCAB_SIZE,
    embedding_dim = 128,
    hidden_dim    = 256,
    latent_dim    = LATENT_DIM,
    max_len       = MAX_LEN,
).to(DEVICE)

# Assign HybridEncoder
vae_hybrid.encoder = HybridEncoder(
    cnn_encoder = pretrained_cnn_encoder,
    gnn_encoder = pretrained_gnn_encoder,
    latent_dim  = LATENT_DIM
).to(DEVICE)

print(f"  Hybrid VAE Total Params : {sum(p.numel() for p in vae_hybrid.parameters()):,}")

# Optimizer only updates fusion layers and decoder since encoder weights are frozen
trainable_params = filter(lambda p: p.requires_grad, vae_hybrid.parameters())
optimizer = optim.AdamW(trainable_params, lr=LR, weight_decay=WEIGHT_DECAY)

# ---------------------------------------------------------------------------
# TRAINING
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("  Training Hybrid VAE")
print("="*60)

train_vae(
    model            = vae_hybrid,
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
    save_path        = os.path.join(CKPT_DIR, "vae_hybrid.pt"),
)

torch.save(vae_hybrid.state_dict(), os.path.join(CKPT_DIR, "vae_hybrid.pt"))

print("\nHybrid VAE training complete")
