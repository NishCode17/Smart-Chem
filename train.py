import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import json
import os
from tqdm import tqdm
from models.vae import VAE

# --- CONFIGURATION (The Control Panel) ---
# CHANGE THIS to "selfies" or "smiles" depending on which model you want to train
MODE = "selfies" 

BATCH_SIZE = 64
LR = 1e-3
EPOCHS = 10  # For MVP, 10-20 is enough. For final result, do 50.
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
processed_dir = "data/processed"
checkpoint_dir = "checkpoints"

os.makedirs(checkpoint_dir, exist_ok=True)

def loss_function(recon_x, x, mu, logvar, kld_weight):
    """
    The Custom ML Loss:
    1. Reconstruction (CrossEntropy): Did we spell the molecule right?
    2. KL Divergence: Is the latent space smooth?
    """
    # Flatten inputs for CrossEntropy
    # x: [batch, seq_len] -> [batch * seq_len]
    # recon_x: [batch, seq_len, vocab_size] -> [batch * seq_len, vocab_size]
    
    # We must reshape to (N, C) for CrossEntropy
    B, L, V = recon_x.shape
    recon_x_flat = recon_x.view(B * L, V)
    x_flat = x.view(B * L)
    
    # 1. Reconstruction Loss
    BCE = F.cross_entropy(recon_x_flat, x_flat, reduction='sum', ignore_index=0) # 0 is padding
    
    # 2. KL Divergence
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + (kld_weight * KLD), BCE, KLD

def train():
    print(f"🚀 Starting Training for: {MODE.upper()}")
    print(f"Device: {DEVICE}")

    # 1. Load Data & Vocab
    vocab_path = f"{processed_dir}/{MODE}_vocab.json"
    data_path = f"{processed_dir}/train_{MODE}.pt"
    
    with open(vocab_path, 'r') as f:
        vocab = json.load(f)
    vocab_size = len(vocab) + 3 # +3 for <pad>, <sos>, <eos> logic in preprocess
    
    print(f"Vocab Size: {vocab_size}")
    
    data = torch.load(data_path)
    # Create DataLoader
    dataset = TensorDataset(data)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
    
    # 2. Initialize Model
    model = VAE(vocab_size=vocab_size, latent_dim=128).to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)
    
    # 3. Training Loop
    model.train()
    for epoch in range(1, EPOCHS + 1):
        total_loss = 0
        total_bce = 0
        total_kld = 0
        
        # KL Annealing Strategy: Slowly increase KLD weight from 0.0 to 1.0
        # This prevents "Posterior Collapse" (The ML defense you will use)
        kld_weight = min(1.0, epoch / 5) 
        
        progress = tqdm(loader, desc=f"Epoch {epoch}/{EPOCHS}")
        
        for batch in progress:
            x = batch[0].to(DEVICE) # Input [batch, seq_len]
            
            optimizer.zero_grad()
            
            # Forward Pass
            recon_x, mu, logvar = model(x)
            
            # Loss Calculation
            loss, bce, kld = loss_function(recon_x, x, mu, logvar, kld_weight)
            
            # Backward Pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            total_bce += bce.item()
            total_kld += kld.item()
            
            progress.set_postfix({"Loss": loss.item() / BATCH_SIZE})

        avg_loss = total_loss / len(loader.dataset)
        print(f"====> Epoch: {epoch} Average loss: {avg_loss:.4f} (BCE: {total_bce/len(loader.dataset):.2f}, KLD: {total_kld/len(loader.dataset):.2f})")
        
        # Save Checkpoint
        if epoch % 5 == 0 or epoch == EPOCHS:
            save_path = f"{checkpoint_dir}/{MODE}_model_epoch_{epoch}.pth"
            torch.save(model.state_dict(), save_path)
            print(f"✅ Saved checkpoint: {save_path}")

if __name__ == "__main__":
    train()