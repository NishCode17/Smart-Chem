"""
run_targeted_eval.py
====================
Experimental evaluation of the targeted molecule generation module for SmartChem.

Evaluates three encoder types (CNN, GNN, Hybrid) across 9 target property
settings using gradient-based latent space optimisation toward QED / LogP / SAS.

Outputs
-------
  evaluation/logs/targeted_generation_results.csv   — per-molecule rows
  evaluation/logs/targeted_generation_summary.csv   — aggregated (encoder × setting)
  evaluation/logs/encoder_comparison_stats.csv      — Mann-Whitney U test results

Usage
-----
  python run_targeted_eval.py
  python run_targeted_eval.py --n_vectors 100 --n_runs 3 --setting T1
"""

import os
import sys
import csv
import json
import math
import time
import uuid
import argparse
import warnings
import itertools
from copy import deepcopy
from dataclasses import dataclass, field, asdict
from typing import Optional

import torch
import torch.optim as optim
import numpy as np

warnings.filterwarnings("ignore")

# ── project root on path ──────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# ── rdkit / selfies ───────────────────────────────────────────────────────────
from rdkit import Chem
from rdkit.Chem import QED as rdQED, Descriptors, DataStructs, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import selfies as sf

# ── project imports ───────────────────────────────────────────────────────────
from models.vae            import VAE
from models.gnn_encoder    import GNNEncoder
from models.hybrid_encoder import HybridEncoder
from models.predictor      import PropertyPredictor
from data.graph_dataset    import smiles_to_graph

try:
    from sascorer import calculateScore as _sas_fn
except ImportError:
    import importlib.util as _ilu
    _s = _ilu.spec_from_file_location("sascorer", os.path.join(ROOT, "sascorer.py"))
    _m = _ilu.module_from_spec(_s)
    _s.loader.exec_module(_m)
    _sas_fn = _m.calculateScore

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
DEVICE        = torch.device("cuda" if torch.cuda.is_available() else "cpu")
CKPT_DIR      = os.path.join(ROOT, "checkpoints")
PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
OUT_DIR       = os.path.join(ROOT, "evaluation", "logs")

LATENT_DIM = 128
NODE_DIM   = 17
EDGE_DIM   = 6

os.makedirs(OUT_DIR, exist_ok=True)

# ── Target property settings ──────────────────────────────────────────────────
TARGET_SETTINGS = {
    "T1": {"desc": "High drug-likeness",        "qed": 0.85, "logp": 2.5, "sas": 2.5, "w_logp": 0.5, "w_sas": 0.5},
    "T2": {"desc": "Moderate drug-likeness",    "qed": 0.65, "logp": 3.0, "sas": 3.5, "w_logp": 0.5, "w_sas": 0.5},
    "T3": {"desc": "Low drug-likeness",         "qed": 0.35, "logp": 4.5, "sas": 5.0, "w_logp": 0.5, "w_sas": 0.5},
    "T4": {"desc": "High QED, high LogP",       "qed": 0.80, "logp": 5.0, "sas": 3.0, "w_logp": 0.5, "w_sas": 0.5},
    "T5": {"desc": "High QED, low SAS",         "qed": 0.80, "logp": 2.5, "sas": 2.0, "w_logp": 0.5, "w_sas": 0.5},
    "T6": {"desc": "Low LogP (hydrophilic)",    "qed": 0.60, "logp": 0.5, "sas": 3.5, "w_logp": 1.0, "w_sas": 0.5},
    "T7": {"desc": "High LogP (lipophilic)",    "qed": 0.50, "logp": 5.5, "sas": 4.0, "w_logp": 1.0, "w_sas": 0.5},
    "T8": {"desc": "Easy synthesis",            "qed": 0.70, "logp": 3.0, "sas": 2.0, "w_logp": 0.5, "w_sas": 1.0},
    "T9": {"desc": "Challenging synthesis",     "qed": 0.60, "logp": 3.5, "sas": 6.5, "w_logp": 0.5, "w_sas": 1.0},
}

TOLERANCES = {"qed": 0.05, "logp": 0.3, "sas": 0.3}

# ── FilterCatalog setup ───────────────────────────────────────────────────────
def _build_filter_catalogs():
    cats = {}
    for name, pname in [
        ("pains",  FilterCatalogParams.FilterCatalogs.PAINS),
        ("brenk",  FilterCatalogParams.FilterCatalogs.BRENK),
        ("nih",    FilterCatalogParams.FilterCatalogs.NIH),
    ]:
        params = FilterCatalogParams()
        params.AddCatalog(pname)
        cats[name] = FilterCatalog(params)
    return cats

FILTER_CATS = _build_filter_catalogs()


# ─────────────────────────────────────────────────────────────────────────────
# CHEMISTRY UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def decode_tokens_to_smiles(token_ids: list, itos: dict) -> Optional[str]:
    """Convert integer token list → SMILES via SELFIES."""
    try:
        tokens = [itos[t] for t in token_ids
                  if t in itos and itos[t] not in ("<pad>", "<sos>", "<eos>")]
        selfies_str = "".join(tokens)
        smiles = sf.decoder(selfies_str)
        return smiles if smiles else None
    except Exception:
        return None


def canonical(smiles: Optional[str]) -> Optional[str]:
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def rdkit_properties(smiles: str):
    """Returns (qed, logp, sas) or (None, None, None) on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    try:
        qed  = rdQED.qed(mol)
        logp = Descriptors.MolLogP(mol)
        sas  = float(_sas_fn(mol))
        return qed, logp, sas
    except Exception:
        return None, None, None


def lipinski_check(smiles: str):
    """Returns (compliant: bool, violations: int) or (None, None)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mw  = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    viols = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return viols == 0, viols


def run_filters(smiles: str):
    """Returns dict of filter results for a valid molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {k: None for k in ["pains_count","pains_flag","brenk_count","brenk_flag","nih_count","nih_flag","filter_pass","lipinski_compliant","lipinski_violations"]}
    results = {}
    total_alerts = 0
    for name in ["pains", "brenk", "nih"]:
        entries = FILTER_CATS[name].GetMatches(mol)
        cnt = len(entries)
        results[f"{name}_count"] = cnt
        results[f"{name}_flag"]  = cnt > 0
        total_alerts += cnt
    lip_ok, lip_viol = lipinski_check(smiles)
    results["lipinski_compliant"]  = lip_ok
    results["lipinski_violations"] = lip_viol
    results["filter_pass"] = (lip_ok is True) and (total_alerts == 0)
    return results


def _infer_vocab_size(ckpt_path: str) -> int:
    """Read vocab size from checkpoint embedding weight shape."""
    with open(ckpt_path, "rb") as fh:
        sd = torch.load(fh, map_location="cpu")
    # CNN/Hybrid: encoder.embedding.weight; GNN: embedding.weight
    for key in ["encoder.embedding.weight", "embedding.weight"]:
        if key in sd:
            return sd[key].shape[0]
    raise RuntimeError("Cannot infer vocab_size from: " + ckpt_path)


def _load_state(model, path):
    with open(path, "rb") as fh:
        model.load_state_dict(torch.load(fh, map_location=DEVICE))


def morgan_fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def tanimoto(fp1, fp2) -> float:
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def tanimoto_max_sim(fp, training_fps) -> float:
    if not training_fps or fp is None:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(fp, training_fps)
    return float(max(sims))


def pairwise_diversity(fps) -> float:
    if len(fps) < 2:
        return 0.0
    n = len(fps)
    total = 0.0
    count = 0
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        total += sum(sims)
        count += i
    avg_sim = total / count if count > 0 else 0.0
    return 1.0 - avg_sim


# -----------------------------------------------------------------------------
# MODEL LOADING
# -----------------------------------------------------------------------------

def load_vocab():
    path = os.path.join(PROCESSED_DIR, "selfies_vocab.json")
    with open(path) as f:
        stoi = json.load(f)
    itos = {v: k for k, v in stoi.items()}
    return stoi, itos


def load_cnn_vae(_unused_vocab_size: int = 0) -> VAE:
    ckpt = os.path.join(CKPT_DIR, "vae_cnn.pt")
    vocab_size = _infer_vocab_size(ckpt)
    model = VAE(vocab_size=vocab_size, embedding_dim=128,
                hidden_dim=256, latent_dim=LATENT_DIM, max_len=100).to(DEVICE)
    _load_state(model, ckpt)
    model.eval()
    print("  [CNN VAE]    loaded: vae_cnn.pt  (vocab=" + str(vocab_size) + ")")
    return model, vocab_size


def _build_legacy_gnn_encoder():
    """Build the 2-layer GINEConv encoder matching the vae_gnn.pt checkpoint."""
    from torch_geometric.nn import GINEConv, global_mean_pool, global_max_pool
    import torch.nn as nn
    class _LegacyGNNEncoder(nn.Module):
        def __init__(self, node_dim=17, edge_dim=6, hid=128, lat=128):
            super().__init__()
            self.node_proj = nn.Linear(node_dim, hid)
            mlp1 = nn.Sequential(nn.Linear(hid, hid), nn.ReLU(), nn.Linear(hid, hid))
            self.conv1 = GINEConv(mlp1, edge_dim=edge_dim)
            self.bn1 = nn.BatchNorm1d(hid)
            mlp2 = nn.Sequential(nn.Linear(hid, hid), nn.ReLU(), nn.Linear(hid, hid))
            self.conv2 = GINEConv(mlp2, edge_dim=edge_dim)
            self.bn2 = nn.BatchNorm1d(hid)
            pool_dim = hid * 2
            self.mu_net     = nn.Sequential(nn.Linear(pool_dim, 128), nn.ReLU(), nn.Linear(128, lat))
            self.logvar_net = nn.Sequential(nn.Linear(pool_dim, 128), nn.ReLU(), nn.Linear(128, lat))
        def forward(self, x, edge_index, edge_attr, batch):
            import torch.nn.functional as F
            h = self.node_proj(x)
            h_in = h; h = self.conv1(h, edge_index, edge_attr); h = self.bn1(h); h = F.relu(h) + h_in
            h_in = h; h = self.conv2(h, edge_index, edge_attr); h = self.bn2(h); h = F.relu(h) + h_in
            h_g = torch.cat([global_mean_pool(h, batch), global_max_pool(h, batch)], dim=1)
            return self.mu_net(h_g), self.logvar_net(h_g).clamp(-6.0, 2.0)
    return _LegacyGNNEncoder(node_dim=NODE_DIM, edge_dim=EDGE_DIM, hid=128, lat=LATENT_DIM).to(DEVICE)


def load_gnn_vae(_unused_vocab_size: int = 0) -> VAE:
    ckpt = os.path.join(CKPT_DIR, "vae_gnn.pt")
    vocab_size = _infer_vocab_size(ckpt)
    model = VAE(vocab_size=vocab_size, embedding_dim=128,
                hidden_dim=256, latent_dim=LATENT_DIM, max_len=100).to(DEVICE)
    model.encoder = _build_legacy_gnn_encoder()
    _load_state(model, ckpt)
    model.eval()
    print("  [GNN VAE]    loaded: vae_gnn.pt  (vocab=" + str(vocab_size) + ", legacy 2-conv)")
    return model, vocab_size


def load_hybrid_vae(_unused_vocab_size: int = 0) -> VAE:
    ckpt = os.path.join(CKPT_DIR, "vae_hybrid.pt")
    vocab_size = _infer_vocab_size(ckpt)
    model = VAE(vocab_size=vocab_size, embedding_dim=128,
                hidden_dim=256, latent_dim=LATENT_DIM, max_len=100).to(DEVICE)
    _gnn = _build_legacy_gnn_encoder()
    model.encoder = HybridEncoder(cnn_encoder=model.encoder,
                                  gnn_encoder=_gnn,
                                  latent_dim=LATENT_DIM).to(DEVICE)
    _load_state(model, ckpt)
    model.eval()
    print("  [Hybrid VAE] loaded: vae_hybrid.pt  (vocab=" + str(vocab_size) + ")")
    return model, vocab_size


def load_predictor(ckpt_name="multi_predictor.pth") -> PropertyPredictor:
    pred = PropertyPredictor(latent_dim=LATENT_DIM).to(DEVICE)
    # Try multi_predictor first (3-output), then fallback options
    for fname in [ckpt_name, "predictor.pt"]:
        path = os.path.join(CKPT_DIR, fname)
        if os.path.exists(path):
            try:
                with open(path, "rb") as fh:
                    sd = torch.load(fh, map_location=DEVICE)
                pred.load_state_dict(sd)
                try:
                    print("  [Predictor]  loaded: " + fname)
                except Exception:
                    pass
                pred.eval()
                return pred
            except Exception as e:
                try:
                    print("  [Predictor]  " + fname + " failed: " + str(e)[:80])
                except Exception:
                    pass
    try:
        print("  [Predictor]  WARNING - no valid checkpoint; using random weights")
    except Exception:
        pass
    pred.eval()
    return pred


def load_training_smiles() -> list:
    path = os.path.join(PROCESSED_DIR, "train_smiles.txt")
    if not os.path.exists(path):
        print("  [Training SMILES] Not found — novelty check disabled")
        return []
    with open(path) as f:
        smiles = [l.strip() for l in f if l.strip()]
    print(f"  [Training SMILES] {len(smiles):,} molecules loaded")
    return smiles


# ─────────────────────────────────────────────────────────────────────────────
# GRADIENT-BASED LATENT OPTIMISATION
# ─────────────────────────────────────────────────────────────────────────────

def optimise_latent(z0: torch.Tensor,
                    predictor: PropertyPredictor,
                    target: dict,
                    max_steps: int = 200,
                    lr: float = 0.05) -> tuple:
    """
    Gradient-based latent vector optimisation toward target properties.

    Returns
    -------
    z_opt      : optimised latent vector (1, latent_dim)
    init_pred  : (qed, logp, sas) predicted at z0
    final_pred : (qed, logp, sas) predicted at z_opt
    n_steps    : steps executed
    final_loss : loss at z_opt
    early_stop : whether early stopping fired
    """
    z = z0.clone().detach().requires_grad_(True)
    optimizer = optim.Adam([z], lr=lr, betas=(0.9, 0.999))

    w_qed  = 1.0
    w_logp = target.get("w_logp", 0.5)
    w_sas  = target.get("w_sas",  0.5)

    t_qed  = target["qed"]
    t_logp = target["logp"]
    t_sas  = target["sas"]

    with torch.no_grad():
        p0 = predictor(z0).squeeze()
        init_pred = (p0[0].item(), p0[1].item(), p0[2].item())

    early_stop = False
    n_steps    = 0
    final_loss = 0.0

    for step in range(max_steps):
        optimizer.zero_grad()
        preds = predictor(z).squeeze()       # (3,)
        loss = (w_qed  * (preds[0] - t_qed)  ** 2 +
                w_logp * (preds[1] - t_logp) ** 2 +
                w_sas  * (preds[2] - t_sas)  ** 2)
        loss.backward()
        optimizer.step()
        n_steps    = step + 1
        final_loss = loss.item()

        # Early stop if within tolerance on all three properties
        with torch.no_grad():
            p = predictor(z).squeeze()
            if (abs(p[0].item() - t_qed)  < TOLERANCES["qed"]  and
                abs(p[1].item() - t_logp) < TOLERANCES["logp"] and
                abs(p[2].item() - t_sas)  < TOLERANCES["sas"]):
                early_stop = True
                break

    with torch.no_grad():
        pf = predictor(z).squeeze()
        final_pred = (pf[0].item(), pf[1].item(), pf[2].item())

    return z.detach(), init_pred, final_pred, n_steps, final_loss, early_stop


# ─────────────────────────────────────────────────────────────────────────────
# DECODE LATENT → SMILES
# ─────────────────────────────────────────────────────────────────────────────

def decode_z(model: VAE, z: torch.Tensor, itos: dict,
             temperature: float = 1.0) -> Optional[str]:
    """Decode a single latent vector to canonical SMILES."""
    model.eval()
    with torch.no_grad():
        z_batch = z.unsqueeze(0) if z.dim() == 1 else z
        tokens  = model._decode_loop(z_batch, DEVICE, temperature)
        token_list = tokens[0].cpu().tolist()
    raw = decode_tokens_to_smiles(token_list, itos)
    return canonical(raw)


# ─────────────────────────────────────────────────────────────────────────────
# PER-MOLECULE ROW BUILDER
# ─────────────────────────────────────────────────────────────────────────────

ROW_FIELDS = [
    "experiment_id", "encoder_type", "target_setting_id", "target_qed",
    "target_logp", "target_sas", "latent_vector_id", "optimization_run_id",
    "random_seed", "initial_pred_qed", "initial_pred_logp", "initial_pred_sas",
    "final_pred_qed", "final_pred_logp", "final_pred_sas",
    "delta_qed", "delta_logp_improvement", "delta_sas_improvement",
    "optimization_steps", "final_opt_loss", "early_stopped",
    "decoded_smiles", "canonical_smiles", "is_valid", "is_unique", "is_novel",
    "rdkit_qed", "rdkit_logp", "rdkit_sas",
    "qed_predictor_error", "logp_predictor_error", "sas_predictor_error",
    "property_hit", "tanimoto_max_sim",
    "pains_count", "pains_flag", "brenk_count", "brenk_flag",
    "nih_count", "nih_flag", "filter_pass",
    "lipinski_compliant", "lipinski_violations",
    "decoding_mode", "heavy_atom_count", "molecular_weight",
]


def _signed_improvement(final_val, init_val, target_val) -> float:
    """Reduction in absolute distance to target = improvement."""
    return abs(init_val - target_val) - abs(final_val - target_val)


def build_row(experiment_id, encoder_type, setting_id, setting,
              vec_id, run_id, seed,
              init_pred, final_pred, n_steps, final_loss, early_stop,
              smiles, can_smi, is_valid, is_unique, is_novel,
              rdk_qed, rdk_logp, rdk_sas, tanimoto,
              filter_results, decoding_mode) -> dict:
    row = {
        "experiment_id":       experiment_id,
        "encoder_type":        encoder_type,
        "target_setting_id":   setting_id,
        "target_qed":          setting["qed"],
        "target_logp":         setting["logp"],
        "target_sas":          setting["sas"],
        "latent_vector_id":    vec_id,
        "optimization_run_id": run_id,
        "random_seed":         seed,
        "initial_pred_qed":    round(init_pred[0], 5),
        "initial_pred_logp":   round(init_pred[1], 5),
        "initial_pred_sas":    round(init_pred[2], 5),
        "final_pred_qed":      round(final_pred[0], 5),
        "final_pred_logp":     round(final_pred[1], 5),
        "final_pred_sas":      round(final_pred[2], 5),
        "delta_qed":                  round(final_pred[0] - init_pred[0], 5),
        "delta_logp_improvement":     round(_signed_improvement(final_pred[1], init_pred[1], setting["logp"]), 5),
        "delta_sas_improvement":      round(_signed_improvement(final_pred[2], init_pred[2], setting["sas"]),  5),
        "optimization_steps":  n_steps,
        "final_opt_loss":      round(final_loss, 6),
        "early_stopped":       early_stop,
        "decoded_smiles":      smiles or "",
        "canonical_smiles":    can_smi or "",
        "is_valid":            is_valid,
        "is_unique":           is_unique,
        "is_novel":            is_novel,
        "rdkit_qed":           round(rdk_qed,  5) if rdk_qed  is not None else "",
        "rdkit_logp":          round(rdk_logp, 5) if rdk_logp is not None else "",
        "rdkit_sas":           round(rdk_sas,  5) if rdk_sas  is not None else "",
        "decoding_mode":       decoding_mode,
    }

    # Predictor calibration errors
    if is_valid and rdk_qed is not None:
        row["qed_predictor_error"]  = round(abs(final_pred[0] - rdk_qed),  5)
        row["logp_predictor_error"] = round(abs(final_pred[1] - rdk_logp), 5)
        row["sas_predictor_error"]  = round(abs(final_pred[2] - rdk_sas),  5)
        hit = (abs(rdk_qed  - setting["qed"])  < TOLERANCES["qed"]  and
               abs(rdk_logp - setting["logp"]) < TOLERANCES["logp"] and
               abs(rdk_sas  - setting["sas"])  < TOLERANCES["sas"])
        row["property_hit"] = hit
    else:
        row["qed_predictor_error"]  = ""
        row["logp_predictor_error"] = ""
        row["sas_predictor_error"]  = ""
        row["property_hit"]         = ""

    row["tanimoto_max_sim"] = round(tanimoto, 5) if tanimoto is not None else ""

    # Filter results
    for k in ["pains_count","pains_flag","brenk_count","brenk_flag",
              "nih_count","nih_flag","filter_pass","lipinski_compliant","lipinski_violations"]:
        row[k] = filter_results.get(k, "")

    # Molecular descriptors
    if is_valid and can_smi:
        mol = Chem.MolFromSmiles(can_smi)
        if mol:
            row["heavy_atom_count"] = mol.GetNumHeavyAtoms()
            row["molecular_weight"] = round(Descriptors.MolWt(mol), 3)
        else:
            row["heavy_atom_count"] = ""
            row["molecular_weight"] = ""
    else:
        row["heavy_atom_count"] = ""
        row["molecular_weight"] = ""

    return row


# ─────────────────────────────────────────────────────────────────────────────
# PER (ENCODER × SETTING) EVALUATION LOOP
# ─────────────────────────────────────────────────────────────────────────────

def run_encoder_setting(
    encoder_type:   str,
    model:          VAE,
    predictor:      PropertyPredictor,
    z_bank:         torch.Tensor,    # (n_vectors, latent_dim)  — shared across encoders
    setting_id:     str,
    setting:        dict,
    n_runs:         int,
    itos:           dict,
    train_smiles_set: set,
    train_fps:      list,
    experiment_id:  str,
    csv_writer,
    decoding_mode:  str = "greedy",
) -> list:
    """
    Run all (vector, run) combinations for one encoder × setting.
    Returns list of row dicts (for summary computation).
    """
    temperature = 1.0 if decoding_mode == "greedy" else 0.8
    n_vectors   = z_bank.shape[0]
    all_rows    = []
    seen_smiles = set()

    t0 = time.time()
    total = n_vectors * n_runs
    done  = 0
    # Suppress RDKit MorganGenerator deprecation warnings in output
    import sys as _sys

    for vec_id in range(n_vectors):
        z0 = z_bank[vec_id].unsqueeze(0).to(DEVICE)   # (1, latent_dim)

        for run_id in range(1, n_runs + 1):
            seed = vec_id * 1000 + run_id
            torch.manual_seed(seed)

            # ── Optimise ────────────────────────────────────────────────────
            z_opt, init_pred, final_pred, n_steps, final_loss, early_stop = \
                optimise_latent(z0, predictor, setting, max_steps=200, lr=0.05)

            # ── Decode ──────────────────────────────────────────────────────
            smiles  = decode_z(model, z_opt.squeeze(0), itos, temperature)
            can_smi = canonical(smiles)
            is_valid = can_smi is not None

            # ── Uniqueness / Novelty ─────────────────────────────────────────
            is_unique = is_novel = False
            if is_valid:
                is_unique = can_smi not in seen_smiles
                if is_unique:
                    seen_smiles.add(can_smi)
                is_novel = can_smi not in train_smiles_set

            # ── RDKit properties ─────────────────────────────────────────────
            rdk_qed = rdk_logp = rdk_sas = None
            if is_valid:
                rdk_qed, rdk_logp, rdk_sas = rdkit_properties(can_smi)

            # ── Tanimoto to training set ─────────────────────────────────────
            tanimoto_sim = None
            if is_valid and is_novel and train_fps:
                fp = morgan_fp(can_smi)
                if fp is not None:
                    tanimoto_sim = tanimoto_max_sim(fp, train_fps)

            # ── Toxicity filters ─────────────────────────────────────────────
            filter_res = {}
            if is_valid:
                filter_res = run_filters(can_smi)

            # ── Build row ────────────────────────────────────────────────────
            row = build_row(
                experiment_id, encoder_type, setting_id, setting,
                vec_id + 1, run_id, seed,
                init_pred, final_pred, n_steps, final_loss, early_stop,
                smiles, can_smi, is_valid, is_unique, is_novel,
                rdk_qed, rdk_logp, rdk_sas, tanimoto_sim,
                filter_res, decoding_mode,
            )
            all_rows.append(row)
            csv_writer.writerow(row)

            done += 1
            if done % 50 == 0 or done == total:
                elapsed = time.time() - t0
                rate    = done / elapsed if elapsed > 0 else 0.001
                remain  = (total - done) / rate if rate > 0 else 0
                n_valid_so_far = sum(1 for r in all_rows if r['is_valid'] is True)
                print("    [{} | {}] {}/{}  valid={}  ETA {:.0f}s".format(
                    encoder_type, setting_id, done, total, n_valid_so_far, remain), end="\r")

    print()
    return all_rows


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY AGGREGATION
# ─────────────────────────────────────────────────────────────────────────────

def pearson_r(xs, ys) -> float:
    xs, ys = np.array(xs), np.array(ys)
    if len(xs) < 2:
        return float("nan")
    return float(np.corrcoef(xs, ys)[0, 1])


def safe_mean(vals) -> float:
    vals = [v for v in vals if v not in (None, "") and not (isinstance(v, float) and math.isnan(v))]
    return float(np.mean(vals)) if vals else float("nan")


def aggregate_rows(rows: list, encoder_type: str, setting_id: str,
                   setting: dict, experiment_id: str) -> dict:
    valid   = [r for r in rows if r["is_valid"] is True]
    unique  = [r for r in valid if r["is_unique"] is True]
    novel   = [r for r in unique if r["is_novel"] is True]
    passed  = [r for r in valid if r["filter_pass"] is True]

    n_total   = len(rows)
    n_valid   = len(valid)
    n_unique  = len(unique)
    n_novel   = len(novel)
    n_pass    = len(passed)

    # Pairwise diversity
    fps = []
    for r in unique:
        fp = morgan_fp(r["canonical_smiles"]) if r["canonical_smiles"] else None
        if fp is not None:
            fps.append(fp)
    diversity = pairwise_diversity(fps) if len(fps) >= 2 else float("nan")

    # Stats helpers
    def pct(sub): return len(sub) / n_total if n_total > 0 else 0.0
    def uq_pct(sub): return len(sub) / n_unique if n_unique > 0 else 0.0
    def vld_pct(sub): return len(sub) / n_valid if n_valid > 0 else 0.0

    # Predictor calibration
    pred_qed  = [r["rdkit_qed"]  for r in valid if r["rdkit_qed"]  != ""]
    pred_logp = [r["rdkit_logp"] for r in valid if r["rdkit_logp"] != ""]
    pred_sas  = [r["rdkit_sas"]  for r in valid if r["rdkit_sas"]  != ""]
    fpq = [r["final_pred_qed"]  for r in valid if r["rdkit_qed"]  != ""]
    fpl = [r["final_pred_logp"] for r in valid if r["rdkit_logp"] != ""]
    fps2= [r["final_pred_sas"]  for r in valid if r["rdkit_sas"]  != ""]

    return {
        "experiment_id":        experiment_id,
        "encoder_type":         encoder_type,
        "target_setting_id":    setting_id,
        "target_qed":           setting["qed"],
        "target_logp":          setting["logp"],
        "target_sas":           setting["sas"],
        "n_total_attempts":     n_total,
        "n_valid":              n_valid,
        "n_unique_valid":       n_unique,
        "n_novel_valid":        n_novel,
        "n_filter_pass":        n_pass,
        "validity_rate":        round(pct(valid),  4),
        "uniqueness_rate":      round(n_unique / n_valid if n_valid > 0 else 0.0, 4),
        "novelty_rate":         round(n_novel / n_unique if n_unique > 0 else 0.0, 4),
        "mean_pairwise_diversity": round(diversity, 4),
        "mean_tanimoto_max_sim": round(safe_mean([r["tanimoto_max_sim"] for r in novel]), 4),
        "mean_delta_qed":        round(safe_mean([r["delta_qed"] for r in rows]), 5),
        "mean_delta_logp":       round(safe_mean([r["delta_logp_improvement"] for r in rows]), 5),
        "mean_delta_sas":        round(safe_mean([r["delta_sas_improvement"] for r in rows]), 5),
        "mean_final_pred_qed":   round(safe_mean([r["final_pred_qed"] for r in rows]), 5),
        "mean_final_pred_logp":  round(safe_mean([r["final_pred_logp"] for r in rows]), 5),
        "mean_final_pred_sas":   round(safe_mean([r["final_pred_sas"] for r in rows]), 5),
        "mean_rdkit_qed":        round(safe_mean(pred_qed),  4),
        "mean_rdkit_logp":       round(safe_mean(pred_logp), 4),
        "mean_rdkit_sas":        round(safe_mean(pred_sas),  4),
        "property_hit_rate":     round(safe_mean([r["property_hit"] for r in valid if r["property_hit"] != ""]), 4),
        "mean_qed_pred_error":   round(safe_mean([r["qed_predictor_error"]  for r in valid if r["qed_predictor_error"]  != ""]), 5),
        "mean_logp_pred_error":  round(safe_mean([r["logp_predictor_error"] for r in valid if r["logp_predictor_error"] != ""]), 5),
        "mean_sas_pred_error":   round(safe_mean([r["sas_predictor_error"]  for r in valid if r["sas_predictor_error"]  != ""]), 5),
        "mean_opt_steps":        round(safe_mean([r["optimization_steps"] for r in rows]), 2),
        "mean_final_opt_loss":   round(safe_mean([r["final_opt_loss"] for r in rows]), 6),
        "early_stop_rate":       round(safe_mean([r["early_stopped"] for r in rows]), 4),
        "lipinski_compliance_rate": round(vld_pct([r for r in valid if r["lipinski_compliant"] is True]), 4),
        "pains_clean_rate":      round(vld_pct([r for r in valid if r["pains_flag"] is False]), 4),
        "brenk_clean_rate":      round(vld_pct([r for r in valid if r["brenk_flag"] is False]), 4),
        "nih_clean_rate":        round(vld_pct([r for r in valid if r["nih_flag"]   is False]), 4),
        "pearson_r_qed":         round(pearson_r(fpq,  pred_qed),  4),
        "pearson_r_logp":        round(pearson_r(fpl,  pred_logp), 4),
        "pearson_r_sas":         round(pearson_r(fps2, pred_sas),  4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MANN-WHITNEY COMPARISON
# ─────────────────────────────────────────────────────────────────────────────

def mann_whitney_comparison(all_rows_by_enc: dict, setting_ids: list,
                            experiment_id: str) -> list:
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        print("  [scipy not found] Skipping Mann-Whitney tests")
        return []

    COMPARE_METRICS = [
        ("validity_rate",       lambda r: float(r["is_valid"]) if isinstance(r["is_valid"], bool) else float("nan")),
        ("delta_qed",           lambda r: float(r["delta_qed"]) if r["delta_qed"] != "" else float("nan")),
        ("delta_logp_improvement", lambda r: float(r["delta_logp_improvement"]) if r["delta_logp_improvement"] != "" else float("nan")),
        ("property_hit",        lambda r: float(r["property_hit"]) if r["property_hit"] not in ("", None) else float("nan")),
        ("tanimoto_max_sim",    lambda r: float(r["tanimoto_max_sim"]) if r["tanimoto_max_sim"] not in ("", None) else float("nan")),
    ]

    encoders = list(all_rows_by_enc.keys())
    pairs = list(itertools.combinations(encoders, 2))
    n_tests = len(pairs) * len(COMPARE_METRICS) * len(setting_ids)

    stat_rows = []
    for s_id in setting_ids:
        enc_data = {}
        for enc in encoders:
            enc_data[enc] = [r for r in all_rows_by_enc[enc] if r["target_setting_id"] == s_id]

        for metric_name, extractor in COMPARE_METRICS:
            for enc_a, enc_b in pairs:
                a_vals = [extractor(r) for r in enc_data[enc_a]]
                b_vals = [extractor(r) for r in enc_data[enc_b]]
                a_clean = [v for v in a_vals if not math.isnan(v)]
                b_clean = [v for v in b_vals if not math.isnan(v)]
                if len(a_clean) < 2 or len(b_clean) < 2:
                    continue
                try:
                    stat, p = mannwhitneyu(a_clean, b_clean, alternative="two-sided")
                    p_bonf = min(p * n_tests, 1.0)
                    stat_rows.append({
                        "experiment_id":    experiment_id,
                        "target_setting_id": s_id,
                        "metric_name":      metric_name,
                        "encoder_a":        enc_a,
                        "encoder_b":        enc_b,
                        "mean_encoder_a":   round(np.nanmean(a_vals), 5),
                        "mean_encoder_b":   round(np.nanmean(b_vals), 5),
                        "mw_u_statistic":   round(stat, 3),
                        "p_value_raw":      round(p, 6),
                        "p_value_bonferroni": round(p_bonf, 6),
                        "significant":      p_bonf < 0.05,
                    })
                except Exception:
                    pass
    return stat_rows


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_vectors",  type=int, default=50,
                        help="Number of initial latent vectors (default 50, paper uses 500)")
    parser.add_argument("--n_runs",     type=int, default=3,
                        help="Optimisation runs per vector (default 3, paper uses 10)")
    parser.add_argument("--setting",    type=str, default=None,
                        help="Single target setting to run (e.g. T1). Default: all")
    parser.add_argument("--seed",       type=int, default=42)
    parser.add_argument("--encoders",   type=str, default="CNN,GNN,Hybrid",
                        help="Comma-separated list of encoders to evaluate")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    experiment_id = str(uuid.uuid4())[:8]
    encoder_list  = [e.strip() for e in args.encoders.split(",")]

    settings_to_run = (
        {args.setting: TARGET_SETTINGS[args.setting]} if args.setting
        else TARGET_SETTINGS
    )

    print("\n" + "="*65)
    print("  SmartChem - Targeted Generation Evaluation")
    print("  Experiment ID : " + experiment_id)
    print("  Device        : " + str(DEVICE))
    print("  Latent vectors: " + str(args.n_vectors) + "  Runs/vector: " + str(args.n_runs))
    print("  Encoders      : " + str(encoder_list))
    print("  Settings      : " + str(list(settings_to_run.keys())))
    print("="*65 + "\n")

    # ── Load support data ──────────────────────────────────────────────────
    stoi, itos = load_vocab()
    # Note: actual vocab_size is inferred per-checkpoint in load_*_vae
    # stoi/itos used only for decoding — use the JSON vocab

    train_smiles_raw = load_training_smiles()
    train_smiles_set = set(train_smiles_raw)

    print("  Computing training fingerprints for Tanimoto (may take ~60s)…")
    train_fps = []
    if train_smiles_raw:
        for smi in train_smiles_raw[:20_000]:   # cap at 20k for speed
            fp = morgan_fp(smi)
            if fp is not None:
                train_fps.append(fp)
        print(f"  Training fingerprints: {len(train_fps):,}")

    predictor = load_predictor()

    # ── Fixed z_bank (shared across encoders per setting for fair comparison) ─
    torch.manual_seed(args.seed)
    z_bank = torch.randn(args.n_vectors, LATENT_DIM)   # CPU — moved to GPU per vector

    # ── CSV writers ───────────────────────────────────────────────────────────
    results_path  = os.path.join(OUT_DIR, "targeted_generation_results.csv")
    summary_path  = os.path.join(OUT_DIR, "targeted_generation_summary.csv")
    stats_path    = os.path.join(OUT_DIR, "encoder_comparison_stats.csv")

    results_need_header = not os.path.exists(results_path) or os.path.getsize(results_path) == 0
    summary_need_header = not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0

    results_file  = open(results_path,  "a", newline="", encoding="utf-8")
    summary_file  = open(summary_path,  "a", newline="", encoding="utf-8")

    results_writer = csv.DictWriter(results_file, fieldnames=ROW_FIELDS, extrasaction="ignore")
    if results_need_header:
        results_writer.writeheader()

    SUMMARY_FIELDS = [
        "experiment_id","encoder_type","target_setting_id","target_qed","target_logp","target_sas",
        "n_total_attempts","n_valid","n_unique_valid","n_novel_valid","n_filter_pass",
        "validity_rate","uniqueness_rate","novelty_rate","mean_pairwise_diversity",
        "mean_tanimoto_max_sim","mean_delta_qed","mean_delta_logp","mean_delta_sas",
        "mean_final_pred_qed","mean_final_pred_logp","mean_final_pred_sas",
        "mean_rdkit_qed","mean_rdkit_logp","mean_rdkit_sas","property_hit_rate",
        "mean_qed_pred_error","mean_logp_pred_error","mean_sas_pred_error",
        "mean_opt_steps","mean_final_opt_loss","early_stop_rate",
        "lipinski_compliance_rate","pains_clean_rate","brenk_clean_rate","nih_clean_rate",
        "pearson_r_qed","pearson_r_logp","pearson_r_sas",
    ]
    summary_writer = csv.DictWriter(summary_file, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
    if summary_need_header:
        summary_writer.writeheader()

    # ── Main evaluation loop ──────────────────────────────────────────────────
    all_rows_by_enc = {e: [] for e in encoder_list}

    for encoder_type in encoder_list:
        print("\n" + "-"*65)
        print("  Loading " + encoder_type + " VAE")
        print("-"*65)

        if encoder_type == "CNN":
            model, _ = load_cnn_vae()
        elif encoder_type == "GNN":
            model, _ = load_gnn_vae()
        elif encoder_type == "Hybrid":
            model, _ = load_hybrid_vae()
        else:
            print("  Unknown encoder '" + encoder_type + "' - skipping")
            continue

        for setting_id, setting in settings_to_run.items():
            print("\n  >> Encoder=" + encoder_type + "  Setting=" + setting_id + " [" + setting['desc'] + "]")
            print("     Target QED=" + str(setting['qed']) + "  LogP=" + str(setting['logp']) + "  SAS=" + str(setting['sas']))

            rows = run_encoder_setting(
                encoder_type=encoder_type,
                model=model,
                predictor=predictor,
                z_bank=z_bank,
                setting_id=setting_id,
                setting=setting,
                n_runs=args.n_runs,
                itos=itos,
                train_smiles_set=train_smiles_set,
                train_fps=train_fps,
                experiment_id=experiment_id,
                csv_writer=results_writer,
                decoding_mode="sample",    # temperature=0.8 — larger, more diverse molecules
            )
            results_file.flush()

            all_rows_by_enc[encoder_type].extend(rows)

            # Aggregate summary
            summary = aggregate_rows(rows, encoder_type, setting_id,
                                     setting, experiment_id)
            summary_writer.writerow(summary)
            summary_file.flush()

            n_v = summary["n_valid"]
            n_t = summary["n_total_attempts"]
            print("    Done | Valid={}/{} ({:.1f}%)  Novel={:.1f}%  Diversity={:.3f}  Accepted={}".format(
                n_v, n_t, summary['validity_rate']*100,
                summary['novelty_rate']*100,
                summary['mean_pairwise_diversity'],
                summary['n_filter_pass']
            ))

        del model
        torch.cuda.empty_cache() if DEVICE.type == "cuda" else None

    results_file.close()
    summary_file.close()

    # ── Statistical comparison ────────────────────────────────────────────────
    print("\n  Computing encoder comparison statistics...")
    stat_rows = mann_whitney_comparison(
        all_rows_by_enc, list(settings_to_run.keys()), experiment_id)
    if stat_rows:
        STAT_FIELDS = ["experiment_id","target_setting_id","metric_name",
                       "encoder_a","encoder_b","mean_encoder_a","mean_encoder_b",
                       "mw_u_statistic","p_value_raw","p_value_bonferroni","significant"]
        with open(stats_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=STAT_FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(stat_rows)
        print(f"  Stats written → {stats_path}")

    print("\n" + "="*65)
    print("  EVALUATION COMPLETE")
    print("  Results  -> " + results_path)
    print("  Summary  -> " + summary_path)
    print("  Stats    -> " + stats_path)
    print("="*65 + "\n")


if __name__ == "__main__":
    main()
