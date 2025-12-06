from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware  # <--- NEW IMPORT
from pydantic import BaseModel
import torch
import json
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.vae import VAE
from backend.chem_utils import get_mol_from_sequence, calculate_properties

app = FastAPI(title="SmartChem API", version="1.0")

# --- NEW: ALLOW REACT TO CONNECT ---
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins (simplest for local dev)
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ... (Keep your existing Config/Model Loading code mostly the same) ...
MODE = "selfies"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
CHECKPOINT_PATH = "checkpoints/selfies_model_epoch_10.pth" # Check your file name!
VOCAB_PATH = f"data/processed/{MODE}_vocab.json"

model = None
idx_to_token = {}

# ... (Keep load_model and decode_tensor functions exactly as they were) ...
# I will skip pasting them to save space, assume they are unchanged. 
# Just make sure they are defined above.

# Redefine the class just in case
class GenerateRequest(BaseModel):
    num_molecules: int = 5

# Copy this Start-up logic if you deleted it, otherwise keep yours
@app.on_event("startup")
def load_model():
    global model, idx_to_token
    # ... (Your existing loading logic) ...
    # PASTE YOUR PREVIOUS LOAD_MODEL CODE HERE
    # Ensure you load vocab and model exactly as before.
    print("⚡ Loading SmartChem Model...")
    if not os.path.exists(VOCAB_PATH): return
    with open(VOCAB_PATH, 'r') as f: vocab = json.load(f)
    idx_to_token = {v: k for k, v in vocab.items()}
    idx_to_token[0] = ""; idx_to_token[1] = ""; idx_to_token[2] = ""
    vocab_size = len(vocab) + 3
    model = VAE(vocab_size=vocab_size, latent_dim=128).to(DEVICE)
    if os.path.exists(CHECKPOINT_PATH):
        model.load_state_dict(torch.load(CHECKPOINT_PATH, map_location=DEVICE))
        model.eval()
        print("✅ Model loaded!")

# ... (Keep decode_tensor) ...
def decode_tensor(tensor_seq):
    tokens = []
    for idx in tensor_seq:
        idx = idx.item()
        if idx == 2: break
        if idx > 2: tokens.append(idx_to_token.get(idx, ""))
    return "".join(tokens)

@app.post("/generate")
def generate(req: GenerateRequest):
    if model is None: raise HTTPException(503, "Model loading...")
    results = []
    with torch.no_grad():
        sampled_indices = model.sample(req.num_molecules, DEVICE)
        for i in range(req.num_molecules):
            seq = decode_tensor(sampled_indices[i])
            mol = get_mol_from_sequence(seq, mode=MODE)
            
            # Helper to generate image base64 directly in backend
            # This makes React logic much simpler
            props = calculate_properties(mol)
            
            # You can keep image generation in frontend (using RDKit JS) 
            # OR send base64 from backend. 
            # For 3 days, sending props is enough, we will render structure in React or just show text.
            # actually, let's update calculate_properties in chem_utils to be safe.
            
            results.append({"sequence": seq, "properties": props})
    return {"data": results}