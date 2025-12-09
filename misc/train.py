import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import json
import os
from tqdm import tqdm
from models.vae import VAE

# --- PRO CONFIGURATION ---
MODE = "selfies" 
BATCH_SIZE = 64  # Keep 64 for RTX 2050 (Safe for 4GB VRAM)
LR = 5e-4        # Lower learning rate for stable, long-term learning
EPOCHS = 100     # <--- THE BIG CHANGE
LATENT_DIM = 128
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

PROCESSED_DIR = "data/processed"
CHECKPOINT_DIR = "checkpoints"
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def loss_function(recon_x, x, mu, logvar, kld_weight):
    B, L, V = recon_x.shape
    recon_x_flat = recon_x.view(B * L, V)
    x_flat = x.view(B * L)
    
    BCE = F.cross_entropy(recon_x_flat, x_flat, reduction='sum', ignore_index=0)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + (kld_weight * KLD), BCE, KLD

def train():
    print(f"🚀 Launching PRO Training Sequence on {DEVICE}")
    
    # Load Data
    vocab_path = f"{PROCESSED_DIR}/{MODE}_vocab.json"
    data_path = f"{PROCESSED_DIR}/train_{MODE}.pt"
    
    with open(vocab_path, 'r') as f: vocab = json.load(f)
    vocab_size = len(vocab) + 3
    
    data = torch.load(data_path)
    # Shuffle is crucial for long training
    dataset = TensorDataset(data)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
    
    print(f"   Dataset Size: {len(data)} molecules")
    print(f"   Estimated Batches per Epoch: {len(loader)}")

    model = VAE(vocab_size=vocab_size, latent_dim=LATENT_DIM).to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)
    
    # Learning Rate Scheduler (Slows down as we get smarter)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

    best_loss = float('inf')
    model.train()
    
    for epoch in range(1, EPOCHS + 1):
        total_loss = 0
        # Slower annealing for 100 epochs (full KL weight at epoch 20)
        kld_weight = min(1.0, epoch / 20)
        
        progress = tqdm(loader, desc=f"Epoch {epoch}/{EPOCHS}", leave=False)
        
        for batch in progress:
            x = batch[0].to(DEVICE)
            optimizer.zero_grad()
            
            recon_x, mu, logvar = model(x)
            loss, bce, kld = loss_function(recon_x, x, mu, logvar, kld_weight)
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0) # Prevent crashes
            optimizer.step()
            
            total_loss += loss.item()
            progress.set_postfix({"loss": loss.item() / BATCH_SIZE})

        avg_loss = total_loss / len(loader.dataset)
        
        # Update LR based on performance
        current_lr = optimizer.param_groups[0]['lr']
        scheduler.step(avg_loss)
        
        print(f"Epoch {epoch}: Loss={avg_loss:.4f} | LR={current_lr:.6f}")
        
        # Save Best Model Logic
        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save(model.state_dict(), f"{CHECKPOINT_DIR}/vae_{MODE}_best.pth")
            print(f"   ⭐ New Best Model Saved! (Loss: {best_loss:.4f})")
            
        # Regular Backup
        if epoch % 10 == 0:
             torch.save(model.state_dict(), f"{CHECKPOINT_DIR}/vae_{MODE}_epoch_{epoch}.pth")

    print(f"✅ Training Complete. Best Loss: {best_loss:.4f}")

if __name__ == "__main__":
    train()