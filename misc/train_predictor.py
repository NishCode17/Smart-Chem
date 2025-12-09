import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from models.vae import VAE
from models.predictor import PropertyPredictor
from backend.chem_utils import get_mol_from_sequence
from rdkit.Chem import Descriptors, QED, rdMolDescriptors
from tqdm import tqdm
import json
import os

# CONFIG
MODE = "selfies"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
VAE_CHECKPOINT = f"checkpoints/vae_{MODE}_best.pth"
VOCAB_PATH = f"data/processed/{MODE}_vocab.json"
DATA_PATH = f"data/processed/train_{MODE}.pt"

def train_predictor():
    print(f"🚀 Training Precision Predictor on {DEVICE} (FULL DATASET)...")

    # 1. Load Resources
    with open(VOCAB_PATH, 'r') as f: vocab = json.load(f)
    vocab_size = len(vocab) + 3
    
    vae = VAE(vocab_size=vocab_size, latent_dim=128).to(DEVICE)
    if not os.path.exists(VAE_CHECKPOINT):
        print("❌ VAE Checkpoint not found.")
        return
    vae.load_state_dict(torch.load(VAE_CHECKPOINT, map_location=DEVICE))
    vae.eval()

    # 2. Generate Training Data
    print("   Generating Property Labels from FULL dataset...")
    raw_data = torch.load(DATA_PATH)
    # NO LIMIT - Use all data for maximum accuracy
    # raw_data = raw_data[:15000] <--- REMOVED
    
    # Batch size 64 for speed
    temp_loader = DataLoader(TensorDataset(raw_data), batch_size=64)
    
    X_train = []
    y_train = [] 

    idx_to_token = {v: k for k, v in vocab.items()}
    idx_to_token[0] = ""; idx_to_token[1] = ""; idx_to_token[2] = ""

    def decode(tensor):
        toks = []
        for idx in tensor:
            if idx.item() == 2: break
            if idx.item() > 2: toks.append(idx_to_token.get(idx.item(), ""))
        return "".join(toks)

    with torch.no_grad():
        for batch in tqdm(temp_loader, desc="Labeling Data"):
            inputs = batch[0].to(DEVICE)
            _, mu, _ = vae(inputs) 
            
            for i in range(len(inputs)):
                seq = decode(inputs[i])
                mol = get_mol_from_sequence(seq, mode=MODE)
                
                if mol:
                    q = QED.qed(mol)
                    l = Descriptors.MolLogP(mol)
                    # Use a cap for SAS/Complexity so it doesn't explode gradients
                    s = rdMolDescriptors.CalcNumRings(mol) + rdMolDescriptors.CalcNumRotatableBonds(mol)
                    
                    X_train.append(mu[i].cpu())
                    y_train.append([q, l, s])

    X_train = torch.stack(X_train).to(DEVICE)
    y_train = torch.tensor(y_train).float().to(DEVICE)
    
    print(f"   Dataset Ready: {len(X_train)} samples")

    # 3. Train with Weighted Loss
    predictor = PropertyPredictor().to(DEVICE)
    optimizer = optim.Adam(predictor.parameters(), lr=1e-3)

    # Weights: Focus heavily on QED (Index 0)
    loss_weights = torch.tensor([10.0, 1.0, 1.0]).to(DEVICE) 

    predictor.train()
    print("   Starting Training Loop...")
    
    # 100 Epochs for deep learning
    for epoch in range(101): 
        optimizer.zero_grad()
        preds = predictor(X_train)
        
        squared_diff = (preds - y_train) ** 2
        weighted_loss = (squared_diff * loss_weights).mean()
        
        weighted_loss.backward()
        optimizer.step()
        
        if epoch % 10 == 0:
            with torch.no_grad():
                mae = torch.abs(preds - y_train).mean(dim=0)
                print(f"   Epoch {epoch}: Loss={weighted_loss.item():.4f} | Avg Errors -> QED:{mae[0]:.3f}, LogP:{mae[1]:.3f}")

    torch.save(predictor.state_dict(), "checkpoints/multi_predictor.pth")
    print("✅ Precision Predictor Saved!")

if __name__ == "__main__":
    train_predictor()