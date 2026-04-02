"""
evaluation/eval_logger.py

Evaluation logger for SmartChem.
"""

import os
import csv
import json
import threading
from datetime import datetime

# Directory layout
_ROOT = os.path.join(os.path.dirname(__file__))  # evaluation/
LOGS_DIR = os.path.join(_ROOT, "logs")
PLOTS_DIR = os.path.join(_ROOT, "plots")

# Log file paths
VAE_LOG        = os.path.join(LOGS_DIR, "vae_training_log.csv")
PREDICTOR_LOG  = os.path.join(LOGS_DIR, "predictor_log.csv")
OPTIM_LOG      = os.path.join(LOGS_DIR, "optimization_log.csv")
VALIDITY_JSON  = os.path.join(LOGS_DIR, "validity_stats.json")
PREDICTOR_TRUE_PRED_LOG = os.path.join(LOGS_DIR, "predictor_true_pred_log.csv")

# Internal helpers
_lock = threading.Lock()


def _ensure_dirs():
    """Create log/plot directories if they do not exist yet."""
    os.makedirs(LOGS_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)


def _append_csv(filepath: str, row: dict, fieldnames: list):
    """
    Thread-safe append of one row to a CSV file.
    Writes the header automatically on the first call.
    """
    _ensure_dirs()
    with _lock:
        file_exists = os.path.isfile(filepath) and os.path.getsize(filepath) > 0
        with open(filepath, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row)


# Public API

def log_vae_epoch(epoch: int, bce_loss: float, kl_loss: float, total_loss: float):
    """
    Append one row to evaluation/logs/vae_training_log.csv.

    Call once per training epoch inside the VAE training loop.

    Parameters
    ----------
    epoch      : 1-based epoch number
    bce_loss   : Binary-Cross-Entropy component of VAE loss
    kl_loss    : KL-Divergence component
    total_loss : bce_loss + kl_loss (or weighted variant)
    """
    fields = ["timestamp", "epoch", "bce_loss", "kl_loss", "total_loss"]
    row = {
        "timestamp": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "epoch": epoch,
        "bce_loss": round(float(bce_loss), 6),
        "kl_loss": round(float(kl_loss), 6),
        "total_loss": round(float(total_loss), 6),
    }
    _append_csv(VAE_LOG, row, fields)


def log_predictor_epoch(epoch: int, mse_qed: float, mse_logp: float, mse_sas: float):
    """
    Append one row to evaluation/logs/predictor_log.csv.

    Call once per training epoch inside the PropertyPredictor training loop.

    Parameters
    ----------
    epoch    : 1-based epoch number
    mse_qed  : MSE loss for QED head
    mse_logp : MSE loss for LogP head
    mse_sas  : MSE loss for SAS head
    """
    fields = ["timestamp", "epoch", "mse_qed", "mse_logp", "mse_sas"]
    row = {
        "timestamp": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "epoch": epoch,
        "mse_qed": round(float(mse_qed), 6),
        "mse_logp": round(float(mse_logp), 6),
        "mse_sas": round(float(mse_sas), 6),
    }
    _append_csv(PREDICTOR_LOG, row, fields)


def log_predictor_sample(true_qed: float, pred_qed: float,
                          true_logp: float, pred_logp: float,
                          true_sas: float, pred_sas: float):
    """
    Append one per-molecule true/predicted row to predictor_true_pred_log.csv.
    Used by plot_metrics.py to generate scatter plots.
    """
    fields = ["timestamp", "true_qed", "pred_qed",
              "true_logp", "pred_logp", "true_sas", "pred_sas"]
    row = {
        "timestamp": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "true_qed":  round(float(true_qed),  4),
        "pred_qed":  round(float(pred_qed),  4),
        "true_logp": round(float(true_logp), 4),
        "pred_logp": round(float(pred_logp), 4),
        "true_sas":  round(float(true_sas),  4),
        "pred_sas":  round(float(pred_sas),  4),
    }
    _append_csv(PREDICTOR_TRUE_PRED_LOG, row, fields)


def log_optimization_step(step: int, qed: float, logp: float, l2_distance: float):
    """
    Append one row to evaluation/logs/optimization_log.csv.

    Call once per gradient-descent step inside optimize_latent_vector().

    Parameters
    ----------
    step        : 0-based step index
    qed         : Mean predicted QED across all latent vectors in this step
    logp        : Mean predicted LogP
    l2_distance : Mean L2 distance from the original (seed) latent vector
    """
    fields = ["timestamp", "step", "qed", "logp", "l2_distance"]
    row = {
        "timestamp":   datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "step":        step,
        "qed":         round(float(qed),         4),
        "logp":        round(float(logp),        4),
        "l2_distance": round(float(l2_distance), 4),
    }
    _append_csv(OPTIM_LOG, row, fields)


def log_validity_stats(total_generated: int,
                       valid_selfies: int,
                       passed_rdkit: int,
                       passed_lipinski: int):
    """
    Write / overwrite evaluation/logs/validity_stats.json with a summary dict.

    Parameters
    ----------
    total_generated  : raw number of SELFIES sequences decoded
    valid_selfies    : sequences that parsed to a non-None RDKit mol
    passed_rdkit     : mols that additionally passed _is_valid_candidate()
    passed_lipinski  : mols that additionally satisfied Lipinski Rule-of-5
    """
    _ensure_dirs()
    stats = {
        "timestamp":        datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "total_generated":  total_generated,
        "valid_selfies":    valid_selfies,
        "passed_rdkit":     passed_rdkit,
        "passed_lipinski":  passed_lipinski,
        "validity_rate":    round(valid_selfies  / max(total_generated, 1), 4),
        "rdkit_pass_rate":  round(passed_rdkit   / max(total_generated, 1), 4),
        "lipinski_rate":    round(passed_lipinski / max(total_generated, 1), 4),
    }
    with _lock:
        with open(VALIDITY_JSON, "w", encoding="utf-8") as fh:
            json.dump(stats, fh, indent=2)
    print(f"[EvalLogger] Validity stats written → {VALIDITY_JSON}")
    return stats
