import torch
import json
import os
import sys
import selfies as sf
from fastapi import HTTPException 

# Internal Project Imports
from models.vae import VAE
from models.predictor import PropertyPredictor
from backend.optimizer import optimize_latent_vector
from backend.chem_utils import get_mol_from_sequence, calculate_properties

# --- CONFIG ---
MODE = "selfies"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
VAE_CHECKPOINT = "checkpoints/vae_selfies_best.pth"
PREDICTOR_CHECKPOINT = "checkpoints/multi_predictor.pth"
VOCAB_PATH = f"data/processed/{MODE}_vocab.json"

# --- GLOBAL STATE ---
_model = None
_predictor = None
_idx_to_token = {}
_token_to_idx = {}

def load_resources():
    """
    Loads VAE, Predictor, and Vocab into global state.
    Idempotent: checks if already loaded.
    """
    global _model, _predictor, _idx_to_token, _token_to_idx
    
    if _model is not None:
        print("⚡ Resources already loaded.")
        return

    print("⚡ [Executor] Loading Resources...")
    if not os.path.exists(VOCAB_PATH): 
        print(f"❌ Critical: Vocab not found at {VOCAB_PATH}")
        return
        
    with open(VOCAB_PATH, 'r') as f: vocab = json.load(f)
    
    # +3 SHIFT
    _idx_to_token = {v + 3: k for k, v in vocab.items()}
    _idx_to_token[0] = ""; _idx_to_token[1] = ""; _idx_to_token[2] = ""
    _token_to_idx = {k: v + 3 for k, v in vocab.items()}
    _token_to_idx["<pad>"] = 0; _token_to_idx["<sos>"] = 1; _token_to_idx["<eos>"] = 2
    
    vocab_size = len(vocab) + 3
    
    _model = VAE(vocab_size, latent_dim=128).to(DEVICE)
    if os.path.exists(VAE_CHECKPOINT):
        _model.load_state_dict(torch.load(VAE_CHECKPOINT, map_location=DEVICE))
        _model.eval()
        print("✅ VAE Loaded")

    _predictor = PropertyPredictor(latent_dim=128).to(DEVICE)
    if os.path.exists(PREDICTOR_CHECKPOINT):
        _predictor.load_state_dict(torch.load(PREDICTOR_CHECKPOINT, map_location=DEVICE))
        _predictor.eval()
        print("✅ Predictor Loaded")

# --- INTERNAL HELPERS ---
def _decode_tensor(tensor_seq):
    tokens = []
    for idx in tensor_seq:
        idx = idx.item()
        if idx == 2: break 
        if idx > 2: tokens.append(_idx_to_token.get(idx, ""))
    return "".join(tokens)

def _smiles_to_latent(smiles):
    try:
        selfies_str = sf.encoder(smiles)
        if not selfies_str: return None
        tokens = list(sf.split_selfies(selfies_str))
        indices = [_token_to_idx.get(t, 0) for t in tokens]
        if len(indices) < 100: indices += [0] * (100 - len(indices))
        else: indices = indices[:100]
        tensor_input = torch.tensor([indices]).long().to(DEVICE)
        with torch.no_grad(): _, mu, _ = _model(tensor_input)
        return mu
    except: return None

def _is_valid_candidate(mol):
    if mol is None: return False
    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms < 5 or num_atoms > 50: return False
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    unique_atoms = set(atoms)
    if len(unique_atoms) == 1 and "C" in unique_atoms:
        return False
    return True

# --- EXPORTED LOGIC ---

def run_lead_optimization(smiles: str):
    """
    Core Logic for OPTIMIZE_LEAD task.
    Returns list of dicts: {sequence, properties, score}
    """
    if _model is None: load_resources()
    
    print(f"\n🔧 [Executor] OPTIMIZING: {smiles[:15]}...")
    
    z_lead = _smiles_to_latent(smiles)
    if z_lead is None:
        raise ValueError("Invalid molecule (encoding failed)")
    
    BATCH_SIZE = 200
    results = []
    
    with torch.no_grad():
        z_batch = z_lead.repeat(BATCH_SIZE, 1)
        
        # Cloud Generation
        noise_close = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.1
        noise_far = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.3
        z_batch[:BATCH_SIZE//2] += noise_close
        z_batch[BATCH_SIZE//2:] += noise_far
        
        # Decode
        indices = _model.decode(z_batch, DEVICE, temperature=0.9)
        
        unique_smiles = {smiles}
        
        for i in range(BATCH_SIZE):
            seq = _decode_tensor(indices[i])
            mol = get_mol_from_sequence(seq, mode=MODE)
            
            if _is_valid_candidate(mol):
                props = calculate_properties(mol)
                if props['valid'] and props['image']:
                    score = props.get('score', 0)
                    s_out = props['smiles']
                    
                    if s_out not in unique_smiles:
                        unique_smiles.add(s_out)
                        props['status'] = f"✨ Optimized ({score})"
                        results.append({
                            "sequence": seq, 
                            "properties": props, 
                            "score": score
                        })
    
    results.sort(key=lambda x: x['score'], reverse=True)
    return results[:5]

def run_random_generation(num_molecules: int):
    """
    Core Logic for GENERATE_RANDOM task.
    """
    if _model is None: load_resources()
    
    results = []
    attempts = 0
    while len(results) < num_molecules and attempts < 100:
        attempts += 1
        with torch.no_grad():
            indices = _model.sample(num_molecules * 5, DEVICE, temperature=2.0)
            for i in range(len(indices)):
                seq = _decode_tensor(indices[i])
                mol = get_mol_from_sequence(seq, mode=MODE)
                if _is_valid_candidate(mol):
                    # Dedup check
                    if any(r['sequence'] == seq for r in results): continue 

                    props = calculate_properties(mol)
                    if props['valid'] and props['image']:
                        results.append({"sequence": seq, "properties": props})
                        if len(results) >= num_molecules: break
    return results

def run_targeted_generation(num_molecules: int, target_qed: float, target_logp: float, target_sas: float):
    """
    Core Logic for GENERATE_TARGETED task.
    """
    if _model is None: load_resources()
    if _predictor is None: raise ValueError("Predictor not loaded")

    target_props = [target_qed, target_logp, target_sas]
    results = []
    
    # 1. Generate Candidate Pool
    INTERNAL_BATCH = 300 
    
    with torch.enable_grad():
        z = torch.randn(INTERNAL_BATCH, 128).to(DEVICE)
        z_opt = optimize_latent_vector(z, _predictor, target_props)
    
    with torch.no_grad():
        indices = _model.decode(z_opt, DEVICE, temperature=0.8)
        for i in range(len(indices)):
            seq = _decode_tensor(indices[i])
            mol = get_mol_from_sequence(seq, mode=MODE)
            
            if _is_valid_candidate(mol):
                props = calculate_properties(mol)
                if props['valid'] and props['image']:
                    score = props.get('score', 0)
                    props['status'] = f"🎯 Targeted (Score {score})"
                    results.append({"sequence": seq, "properties": props, "score": score})
    
    # Sort
    results.sort(key=lambda x: x['score'], reverse=True)
    
    # Fallback
    if not results: 
        return run_random_generation(num_molecules)
        
    return results[:num_molecules]
