import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import json
import os
from tqdm import tqdm
from models.vae import VAE  # Ensure models/vae.py exists

# --- CONFIGURATION ---
# Set this to 'selfies' or 'smiles'
MODE = "selfies" 

# Hyperparameters
BATCH_SIZE = 64
LR = 1e-3
EPOCHS = 20  # 20 Epochs is the sweet spot for a 3-day MVP
LATENT_DIM = 128

# Device Handling (Autodetects your RTX 2050)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Paths
PROCESSED_DIR = "data/processed"
CHECKPOINT_DIR = "checkpoints"
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def loss_function(recon_x, x, mu, logvar, kld_weight):
    """
    Custom VAE Loss: Reconstruction + KL Divergence
    """
    # Reshape for CrossEntropy: [Batch * Seq, Vocab]
    B, L, V = recon_x.shape
    recon_x_flat = recon_x.view(B * L, V)
    x_flat = x.view(B * L)
    
    # 1. Reconstruction Loss (did we predict the right token?)
    # ignore_index=0 skips padding (so we don't learn to predict empty space)
    BCE = F.cross_entropy(recon_x_flat, x_flat, reduction='sum', ignore_index=0)
    
    # 2. KL Divergence (Regularization)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + (kld_weight * KLD), BCE, KLD

def train():
    print(f"🚀 Launching Training on: {DEVICE}")
    if DEVICE.type == 'cuda':
        print(f"   GPU: {torch.cuda.get_device_name(0)}")
    
    # 1. Load the Vocabulary & Data
    # (Created by preprocess.py from your train_molecules.csv)
    vocab_path = f"{PROCESSED_DIR}/{MODE}_vocab.json"
    data_path = f"{PROCESSED_DIR}/train_{MODE}.pt"
    
    if not os.path.exists(data_path):
        print(f"❌ Error: Could not find {data_path}. Did you run preprocess.py?")
        return

    # Load Vocab to get size
    with open(vocab_path, 'r') as f:
        vocab = json.load(f)
    # +3 handles <pad>, <sos>, <eos> which we added in preprocessing
    vocab_size = len(vocab) + 3 
    
    print(f"   Dataset: {MODE.upper()}")
    print(f"   Vocab Size: {vocab_size}")

    # Load Tensors
    data = torch.load(data_path)
    # Limit check: If you ran preprocess with MVP_LIMIT, this is already small.
    dataset = TensorDataset(data)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
    
    # 2. Initialize the Model
    model = VAE(vocab_size=vocab_size, latent_dim=LATENT_DIM).to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)
    
    # 3. Training Loop
    model.train()
    print("   Starting Training Loop...")
    
    for epoch in range(1, EPOCHS + 1):
        total_loss = 0
        total_kld = 0
        
        # Annealing: Slowly introduce KL term to prevent "Posterior Collapse"
        kld_weight = min(1.0, epoch / 5)
        
        progress = tqdm(loader, desc=f"Epoch {epoch}/{EPOCHS}", leave=False)
        
        for batch in progress:
            x = batch[0].to(DEVICE) # Shape: [Batch, Seq_Len]
            
            optimizer.zero_grad()
            
            # Forward Pass
            recon_x, mu, logvar = model(x)
            
            # Loss
            loss, bce, kld = loss_function(recon_x, x, mu, logvar, kld_weight)
            
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            total_kld += kld.item()
            
            # Update progress bar
            progress.set_postfix({"loss": loss.item() / BATCH_SIZE})

        # Epoch Summary
        avg_loss = total_loss / len(loader.dataset)
        avg_kld = total_kld / len(loader.dataset)
        print(f"Epoch {epoch}: Avg Loss = {avg_loss:.4f} (KLD: {avg_kld:.4f})")
        
        # Save Checkpoint (Overwrites previous to save space for MVP)
        torch.save(model.state_dict(), f"{CHECKPOINT_DIR}/vae_{MODE}_latest.pth")

    print(f"✅ Training Complete! Model saved to {CHECKPOINT_DIR}/vae_{MODE}_latest.pth")

if __name__ == "__main__":
    train()