import torch
import numpy as np
import json
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.vae import VAE
from backend.chem_utils import get_mol_from_sequence
from rdkit import Chem

# Config
MODE = "selfies" 
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
CHECKPOINT_DIR = "checkpoints"
MODEL_PATH = f"{CHECKPOINT_DIR}/vae_{MODE}_best.pth"
VOCAB_PATH = f"data/processed/{MODE}_vocab.json"
DATA_PATH = f"data/processed/train_{MODE}.pt"

def slerp(val, low, high):
    """
    Spherical Linear Interpolation between two vectors.
    """
    low_norm = low / torch.norm(low, dim=1, keepdim=True)
    high_norm = high / torch.norm(high, dim=1, keepdim=True)
    
    # Calculate cosine of angle between vectors
    omega = torch.acos(torch.clamp(torch.matmul(low_norm, high_norm.t()), -1, 1))
    so = torch.sin(omega)
    
    # Handle case where vectors are too close (sin(omega) ~ 0)
    if so == 0:
        return (1.0 - val) * low + val * high
    
    return (torch.sin((1.0 - val) * omega) / so) * low + (torch.sin(val * omega) / so) * high

def interpolate():
    print(f"Latent Space Interpolation (SLERP) Analysis")
    print(f"   Device: {DEVICE}")

    # 1. Load resources
    with open(VOCAB_PATH, 'r') as f: vocab = json.load(f)
    # Invert vocabulary
    idx_to_token = {v: k for k, v in vocab.items()}
    # Handle special tokens
    idx_to_token[0] = ""; idx_to_token[1] = ""; idx_to_token[2] = ""
    
    vocab_size = len(vocab) + 3
    
    # 2. Load model
    model = VAE(vocab_size=vocab_size, latent_dim=128).to(DEVICE)
    if not os.path.exists(MODEL_PATH):
        print(f"❌ Model not found at {MODEL_PATH}")
        return
        
    model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    model.eval()
    print("VAE Model Loaded")

    # 3. Load random samples
    data = torch.load(DATA_PATH)
    idx1, idx2 = np.random.randint(0, len(data), 2)
    
    batch = torch.stack([data[idx1], data[idx2]]).to(DEVICE)
    
    # 4. Get latent vectors
    with torch.no_grad():
        # Encode input
        embedded = model.embedding(batch).permute(0, 2, 1)
        c1 = torch.nn.functional.relu(model.conv1(embedded))
        c2 = torch.nn.functional.relu(model.conv2(c1))
        c3 = torch.nn.functional.relu(model.conv3(c2))                 
        pooled = model.adaptive_pool(c3).squeeze(2)  
        mu = model.fc_mu(pooled)
        # Deterministic representation using mean
    
    z_start = mu[0].unsqueeze(0)
    z_end = mu[1].unsqueeze(0)

    print(f"\nInterpolating between:")
    print(f"   Start Mol Index: {idx1}")
    print(f"   End Mol Index:   {idx2}")

    # 5. Perform interpolation
    steps = 10
    print(f"\n--- Interpolation Path ({steps} steps) ---")
    
    for i in range(steps + 1):
        val = i / steps
        z_interp = slerp(val, z_start, z_end)
        
        # Decode
        with torch.no_grad():
            decoded_indices = model.decode(z_interp, DEVICE, temperature=0.1)
        
        # Convert to string
        seq = ""
        for idx in decoded_indices[0]:
            if idx.item() == 2: break # EOS
            if idx.item() > 2: seq += idx_to_token.get(idx.item(), "")
            
        # Check validity
        mol = get_mol_from_sequence(seq, mode=MODE)
        valid = "Valid" if mol else "Invalid"
        smiles = Chem.MolToSmiles(mol) if mol else "N/A"
        
        if i == 0: tag = "(START)"
        elif i == steps: tag = "(END)"
        else: tag = f"({int(val*100)}%)"

        print(f"{tag:<10} | {valid} | {smiles[:60]}...")

    print("\nAnalysis Complete. If intermediate points are mostly valid, the space is smooth.")

if __name__ == "__main__":
    interpolate()
