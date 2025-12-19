
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import torch
import json
import os
import sys
import selfies as sf
from typing import Optional
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.vae import VAE
from models.predictor import PropertyPredictor
from backend.optimizer import optimize_latent_vector
from backend.chem_utils import get_mol_from_sequence, calculate_properties, get_3d_mol_block

app = FastAPI(title="SmartChem API", version="9.0")

# --- CONFIG ---
MODE = "selfies"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
VAE_CHECKPOINT = "checkpoints/vae_selfies_best.pth"
PREDICTOR_CHECKPOINT = "checkpoints/multi_predictor.pth"
VOCAB_PATH = f"data/processed/{MODE}_vocab.json"

model = None
predictor = None
idx_to_token = {}
token_to_idx = {}

app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

class GenerateRequest(BaseModel):
    num_molecules: int = 5
    target_qed: float = 0.8
    target_logp: float = 2.5
    target_sas: float = 3.0

class OptimizeRequest(BaseModel):
    smiles: str 
    target_qed: float
    target_logp: float

@app.on_event("startup")
def load_resources():
    global model, predictor, idx_to_token, token_to_idx
    print("⚡ Loading Resources...")
    if not os.path.exists(VOCAB_PATH): 
        print(f"❌ Critical: Vocab not found at {VOCAB_PATH}")
        return
        
    with open(VOCAB_PATH, 'r') as f: vocab = json.load(f)
    
    # +3 SHIFT (Correct for your model)
    idx_to_token = {v + 3: k for k, v in vocab.items()}
    idx_to_token[0] = ""; idx_to_token[1] = ""; idx_to_token[2] = ""
    token_to_idx = {k: v + 3 for k, v in vocab.items()}
    token_to_idx["<pad>"] = 0; token_to_idx["<sos>"] = 1; token_to_idx["<eos>"] = 2
    
    vocab_size = len(vocab) + 3
    print(f"   Vocab Size: {vocab_size} (Indices shifted +3)")
    
    model = VAE(vocab_size, latent_dim=128).to(DEVICE)
    if os.path.exists(VAE_CHECKPOINT):
        model.load_state_dict(torch.load(VAE_CHECKPOINT, map_location=DEVICE))
        model.eval()
        print("✅ VAE Loaded")

    predictor = PropertyPredictor(latent_dim=128).to(DEVICE)
    if os.path.exists(PREDICTOR_CHECKPOINT):
        predictor.load_state_dict(torch.load(PREDICTOR_CHECKPOINT, map_location=DEVICE))
        predictor.eval()
        print("✅ Predictor Loaded")

def decode_tensor(tensor_seq):
    tokens = []
    for idx in tensor_seq:
        idx = idx.item()
        if idx == 2: break 
        if idx > 2: tokens.append(idx_to_token.get(idx, ""))
    return "".join(tokens)

def smiles_to_latent(smiles):
    try:
        selfies_str = sf.encoder(smiles)
        if not selfies_str: return None
        tokens = list(sf.split_selfies(selfies_str))
        indices = [token_to_idx.get(t, 0) for t in tokens]
        if len(indices) < 100: indices += [0] * (100 - len(indices))
        else: indices = indices[:100]
        tensor_input = torch.tensor([indices]).long().to(DEVICE)
        with torch.no_grad(): _, mu, _ = model(tensor_input)
        return mu
    except: return None

# --- STRICT FILTER (Anti-Carbon) ---
def is_valid_candidate(mol):
    if mol is None: return False
    
    # 1. Size Rule (STRICT)
    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms < 5 or num_atoms > 50: return False
    
    # 2. Composition Rule (NO PURE CARBON)
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    unique_atoms = set(atoms)
    
    # If the ONLY atom type is Carbon, REJECT IT.
    if len(unique_atoms) == 1 and "C" in unique_atoms:
        return False
        
    return True

@app.post("/optimize/lead")
def optimize_lead(req: OptimizeRequest):
    """
    Neighborhood Scanning Strategy:
    1. Encode Lead -> Z
    2. Create cloud of random noise around Z
    3. Decode all and pick the best one
    """
    if model is None: raise HTTPException(503, "AI error")
    
    print(f"\n🔧 OPTIMIZING: {req.smiles[:15]}...")
    
    z_lead = smiles_to_latent(req.smiles)
    if z_lead is None: raise HTTPException(400, "Invalid molecule")

    results = []
    
    # 1. Neighborhood Search (200 Attempts)
    BATCH_SIZE = 200
    with torch.no_grad():
        z_batch = z_lead.repeat(BATCH_SIZE, 1)
        
        # Create a "Cloud" of variations
        noise_close = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.1
        noise_far = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.3
        
        z_batch[:BATCH_SIZE//2] += noise_close
        z_batch[BATCH_SIZE//2:] += noise_far
        
        # 2. Decode All
        indices = model.decode(z_batch, DEVICE, temperature=0.9)
        
        for i in range(BATCH_SIZE):
            seq = decode_tensor(indices[i])
            mol = get_mol_from_sequence(seq, mode=MODE)
            
            # 3. Filter & Score
            if is_valid_candidate(mol):
                props = calculate_properties(mol)
                if props['valid'] and props['image']:
                    
                    # Score using new Fitness Function
                    # Higher score is better
                    score = props.get('score', 0)
                    
                    # Only accept if it's NOT the exact same input
                    if props['smiles'] != req.smiles:
                        props['status'] = f"✨ Optimized ({score})"
                        results.append({"sequence": seq, "properties": props, "score": score})

    # Sort by Score (Descending)
    results.sort(key=lambda x: x['score'], reverse=True)
    
    # Dedup
    unique = []
    seen = set()
    seen.add(req.smiles)
    for r in results:
        s = r['properties']['smiles']
        if s not in seen:
            seen.add(s)
            unique.append(r)
            
    print(f"   ✅ Found {len(unique)} valid variations.")
    return {"data": unique[:5]}

# --- GENERATOR (Unchanged) ---
@app.post("/generate")
def generate_random(req: GenerateRequest):
    if model is None: raise HTTPException(503, "Model error")
    results = []
    attempts = 0
    # Search aggressively (up to 100 batches)
    while len(results) < req.num_molecules and attempts < 100:
        attempts += 1
        with torch.no_grad():
            indices = model.sample(req.num_molecules * 5, DEVICE, temperature=2.0)
            for i in range(len(indices)):
                seq = decode_tensor(indices[i])
                mol = get_mol_from_sequence(seq, mode=MODE)
                if is_valid_candidate(mol):
                    # Dedup check
                    dup = False
                    for r in results:
                         if r['sequence'] == seq: 
                             dup = True
                             break
                    if dup: continue

                    props = calculate_properties(mol)
                    if props['valid'] and props['image']:
                        results.append({"sequence": seq, "properties": props})
                        if len(results) >= req.num_molecules: break
    return {"data": results}

@app.post("/generate/targeted")
def generate_targeted(req: GenerateRequest):
    if model is None or predictor is None: raise HTTPException(503, "AI error")
    
    # Updated Optimization Weights (User Request)
    # QED: 8, LogP: 2, SAS: 1
    # Note: Optimization function logic is internal, but we pass targets here.
    target_props = [req.target_qed, req.target_logp, req.target_sas]
    results = []
    print(f"\n🎯 TARGET: QED={req.target_qed} | LogP={req.target_logp}")
    
    # 1. Generate Massive Candidate Pool (Ensures Hard Filters Pass)
    INTERNAL_BATCH = 300 
    
    with torch.enable_grad():
        z = torch.randn(INTERNAL_BATCH, 128).to(DEVICE)
        # We rely on optimizer to use correct loss weights internally if adjustable,
        # otherwise standard optimization applies.
        z_opt = optimize_latent_vector(z, predictor, target_props)
    
    with torch.no_grad():
        indices = model.decode(z_opt, DEVICE, temperature=0.8)
        for i in range(len(indices)):
            seq = decode_tensor(indices[i])
            mol = get_mol_from_sequence(seq, mode=MODE)
            
            # 2. Hard Filters
            if is_valid_candidate(mol):
                props = calculate_properties(mol)
                if props['valid'] and props['image']:
                    
                    # 3. Soft Scoring (Using chem_utils.score_molecule)
                    score = props.get('score', 0)
                    
                    props['status'] = f"🎯 Targeted (Score {score})"
                    results.append({"sequence": seq, "properties": props, "score": score})
    
    # 4. Sort by Score (Descending)
    results.sort(key=lambda x: x['score'], reverse=True)
    
    # 5. Fallback Guarantee
    if not results: 
        print("⚠️ Targeted gen failed filters, falling back to Random.")
        return generate_random(req)
        
    return {"data": results[:req.num_molecules]}

# --- Routers ---
from backend.routers import auth, projects, molecules

app.include_router(auth.router)
app.include_router(projects.router)
app.include_router(molecules.router)
from backend.routers import jobs
app.include_router(jobs.router)

from backend.assistant import ChatRequest, get_groq_response

@app.post("/assistant/chat")
def chat_assistant(req: ChatRequest):
    response = get_groq_response(req.message, req.context)
    return {"reply": response}


class StructureRequest(BaseModel):
    smiles: str

@app.post("/utils/3d")
def get_3d_structure(req: StructureRequest):
    mol_block = get_3d_mol_block(req.smiles)
    if not mol_block: raise HTTPException(400, "Could not generate 3D structure")
    return {"mol_block": mol_block}