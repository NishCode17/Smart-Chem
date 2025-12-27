
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
    # Use Shared Executor to load resources (idempotent)
    import backend.ml_executor as ml_exec
    ml_exec.load_resources()

# --- STRICT FILTER IS NOW IN ml_executor (Imported internally) ---

@app.post("/optimize/lead")
def optimize_lead(req: OptimizeRequest):
    """
    Synchronous Endpoint (Legacy)
    Delegates to shared executor.
    """
    try:
        import backend.ml_executor as ml_exec
        results = ml_exec.run_lead_optimization(req.smiles)
        print(f"   ✅ Found {len(results)} valid variations.")
        return {"data": results}
    except Exception as e:
        raise HTTPException(500, str(e))

# --- GENERATOR ---
@app.post("/generate")
def generate_random(req: GenerateRequest):
    try:
        import backend.ml_executor as ml_exec
        results = ml_exec.run_random_generation(req.num_molecules)
        return {"data": results}
    except Exception as e:
        raise HTTPException(500, str(e))

@app.post("/generate/targeted")
def generate_targeted(req: GenerateRequest):
    try:
        import backend.ml_executor as ml_exec
        results = ml_exec.run_targeted_generation(
            req.num_molecules, 
            req.target_qed, 
            req.target_logp, 
            req.target_sas
        )
        return {"data": results}
    except Exception as e:
        raise HTTPException(500, str(e))


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