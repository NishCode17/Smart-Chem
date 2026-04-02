"""
scripts/run_training.py

Evaluation training pipeline entry point over preprocessed dataset.
Outputs are saved to `evaluation/logs` and `evaluation/plots`.
"""

import os
import sys
import json
import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm

# Add project root to path
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from models.vae import VAE
from models.predictor import PropertyPredictor
from evaluation.eval_logger import (
    log_vae_epoch,
    log_predictor_epoch,
    log_predictor_sample,
)

# Configuration
MODE          = "selfies"
BATCH_SIZE    = 64
VAE_EPOCHS    = 8
PRED_EPOCHS   = 8
LR_VAE        = 5e-4
LR_PRED       = 1e-3
LATENT_DIM    = 128
MAX_LEN       = 100

DEVICE        = torch.device("cuda" if torch.cuda.is_available() else "cpu")
PROCESSED_DIR = os.path.join(ROOT, "data", "processed")
CKPT_DIR      = os.path.join(ROOT, "checkpoints", "eval_run")
PROD_CKPT_DIR = os.path.join(ROOT, "checkpoints")

os.makedirs(CKPT_DIR, exist_ok=True)


def vae_loss(recon_x, x, mu, logvar, kld_weight):
    B, L, V = recon_x.shape
    recon_flat = recon_x.view(B * L, V)
    x_flat     = x.view(B * L)
    bce = F.cross_entropy(recon_flat, x_flat, reduction="sum", ignore_index=0)
    kld = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return bce + kld_weight * kld, bce, kld


def train_vae(vocab_size, loader):
    print(f"\n{'='*60}")
    print(f"  STAGE 1 — VAE Training  ({VAE_EPOCHS} epochs, {DEVICE})")
    print(f"{'='*60}")

    model = VAE(vocab_size=vocab_size, latent_dim=LATENT_DIM).to(DEVICE)

    # Initialize checkpoint if available
    prod_ckpt = os.path.join(PROD_CKPT_DIR, f"vae_{MODE}_best.pth")
    if os.path.exists(prod_ckpt):
        model.load_state_dict(torch.load(prod_ckpt, map_location=DEVICE))
        print(f"  ✅ Warm-started from {prod_ckpt}")
    else:
        print("  ℹ️  No existing checkpoint — training from scratch")

    optimizer = optim.Adam(model.parameters(), lr=LR_VAE)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", patience=3, factor=0.5
    )

    model.train()
    best_loss = float("inf")

    for epoch in range(1, VAE_EPOCHS + 1):
        # KL annealing
        kld_weight = min(1.0, epoch / 4.0)

        running_bce   = 0.0
        running_kld   = 0.0
        running_total = 0.0
        n_samples     = 0

        pbar = tqdm(loader, desc=f"  VAE Epoch {epoch}/{VAE_EPOCHS}", leave=False,
                    ncols=90)
        for (x,) in pbar:
            x = x.to(DEVICE)
            optimizer.zero_grad()

            recon_x, mu, logvar = model(x)
            loss, bce, kld = vae_loss(recon_x, x, mu, logvar, kld_weight)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            bs = x.size(0)
            running_bce   += bce.item()
            running_kld   += kld.item()
            running_total += loss.item()
            n_samples     += bs
            pbar.set_postfix({"loss": f"{loss.item()/bs:.4f}"})

        # Calculate epoch averages
        avg_bce   = running_bce   / n_samples
        avg_kld   = running_kld   / n_samples
        avg_total = running_total / n_samples

        scheduler.step(avg_total)
        lr_now = optimizer.param_groups[0]["lr"]

        # Validation check
        assert not (torch.isnan(torch.tensor(avg_total)) or avg_total > 1e8), \
            f"Loss exploded at epoch {epoch}: {avg_total}"

        print(f"  VAE  Epoch {epoch:02d}/{VAE_EPOCHS} | "
              f"BCE={avg_bce:.4f}  KL={avg_kld:.4f}  "
              f"Total={avg_total:.4f}  LR={lr_now:.2e}")

        # Log evaluation metrics
        log_vae_epoch(
            epoch=epoch,
            bce_loss=avg_bce,
            kl_loss=avg_kld,
            total_loss=avg_total,
        )

        # Save checkpoint
        if avg_total < best_loss:
            best_loss = avg_total
            torch.save(
                model.state_dict(),
                os.path.join(CKPT_DIR, f"vae_{MODE}_eval_best.pth"),
            )

    torch.save(
        model.state_dict(),
        os.path.join(CKPT_DIR, f"vae_{MODE}_eval_final.pth"),
    )
    print(f"\n  ✅ VAE training complete.  Best loss: {best_loss:.4f}")
    print(f"     Weights saved → {CKPT_DIR}/")
    return model


def build_predictor_dataset(vae, loader):
    """
    Run the VAE encoder over all training molecules to get μ vectors,
    then compute true QED / LogP / SAS labels with RDKit.
    Returns (X_tensor, y_tensor) on CPU.
    """
    from rdkit.Chem import QED as rdQED, Descriptors, rdMolDescriptors
    from backend.chem_utils import get_mol_from_sequence
    import json

    with open(os.path.join(PROCESSED_DIR, f"{MODE}_vocab.json")) as fh:
        vocab = json.load(fh)

    # Build token mapping
    idx_to_token = {v + 3: k for k, v in vocab.items()}
    idx_to_token[0] = ""
    idx_to_token[1] = ""
    idx_to_token[2] = ""

    def decode(tensor_row):
        tokens = []
        for idx in tensor_row:
            idx = idx.item()
            if idx == 2:
                break
            if idx > 2:
                tokens.append(idx_to_token.get(idx, ""))
        return "".join(tokens)

    X_list, y_list = [], []
    vae.eval()

    with torch.no_grad():
        for (x,) in tqdm(loader, desc="  Labelling molecules", leave=False, ncols=90):
            x = x.to(DEVICE)
            _, mu, _ = vae(x)

            for i in range(x.size(0)):
                seq = decode(x[i])
                mol = get_mol_from_sequence(seq, mode=MODE)
                if mol is None:
                    continue

                try:
                    q = rdQED.qed(mol)
                    l = Descriptors.MolLogP(mol)
                    # Calculate metrics
                    s = (rdMolDescriptors.CalcNumRings(mol) +
                         rdMolDescriptors.CalcNumRotatableBonds(mol))
                    X_list.append(mu[i].cpu())
                    y_list.append([q, l, float(s)])
                except Exception:
                    continue

    if not X_list:
        raise RuntimeError("No valid molecules found — check your dataset / vocab.")

    X = torch.stack(X_list)
    y = torch.tensor(y_list, dtype=torch.float32)
    print(f"  ✅ Predictor dataset: {len(X)} samples  "
          f"| QED μ={y[:,0].mean():.3f}  LogP μ={y[:,1].mean():.3f}")
    return X, y


def train_predictor(vae, loader):
    print(f"\n{'='*60}")
    print(f"  STAGE 2 — Predictor Training  ({PRED_EPOCHS} epochs)")
    print(f"{'='*60}")

    X, y = build_predictor_dataset(vae, loader)
    X, y = X.to(DEVICE), y.to(DEVICE)

    predictor  = PropertyPredictor(latent_dim=LATENT_DIM).to(DEVICE)
    optimizer  = optim.Adam(predictor.parameters(), lr=LR_PRED)
    weights    = torch.tensor([10.0, 1.0, 1.0], device=DEVICE)

    pred_loader = DataLoader(
        TensorDataset(X, y), batch_size=256, shuffle=True
    )

    predictor.train()

    for epoch in range(1, PRED_EPOCHS + 1):
        epoch_mse_qed  = 0.0
        epoch_mse_logp = 0.0
        epoch_mse_sas  = 0.0
        n_batches = 0

        for x_batch, y_batch in pred_loader:
            optimizer.zero_grad()
            preds = predictor(x_batch)

            # Compute loss
            mse_qed  = F.mse_loss(preds[:, 0], y_batch[:, 0])
            mse_logp = F.mse_loss(preds[:, 1], y_batch[:, 1])
            mse_sas  = F.mse_loss(preds[:, 2], y_batch[:, 2])

            # Compute weighted loss
            sq_diff      = (preds - y_batch) ** 2
            weighted_loss = (sq_diff * weights).mean()

            weighted_loss.backward()
            optimizer.step()

            epoch_mse_qed  += mse_qed.item()
            epoch_mse_logp += mse_logp.item()
            epoch_mse_sas  += mse_sas.item()
            n_batches      += 1

        avg_mse_qed  = epoch_mse_qed  / n_batches
        avg_mse_logp = epoch_mse_logp / n_batches
        avg_mse_sas  = epoch_mse_sas  / n_batches

        print(f"  Pred Epoch {epoch:02d}/{PRED_EPOCHS} | "
              f"MSE_QED={avg_mse_qed:.5f}  "
              f"MSE_LogP={avg_mse_logp:.5f}  "
              f"MSE_SAS={avg_mse_sas:.5f}")

        # Log evaluation metrics
        log_predictor_epoch(
            epoch=epoch,
            mse_qed=avg_mse_qed,
            mse_logp=avg_mse_logp,
            mse_sas=avg_mse_sas,
        )

    # Log samples
    print("  Logging per-sample true/pred values for scatter plot...")
    predictor.eval()
    with torch.no_grad():
        # Use a single pass over the full dataset
        all_preds = predictor(X)
        n = min(len(X), 500)   # log up to 500 samples (keeps CSV manageable)
        for i in range(n):
            log_predictor_sample(
                true_qed=y[i, 0].item(),  pred_qed=all_preds[i, 0].item(),
                true_logp=y[i, 1].item(), pred_logp=all_preds[i, 1].item(),
                true_sas=y[i, 2].item(),  pred_sas=all_preds[i, 2].item(),
            )

    torch.save(
        predictor.state_dict(),
        os.path.join(CKPT_DIR, "predictor_eval.pth"),
    )
    print(f"\n  ✅ Predictor training complete.  Weights → {CKPT_DIR}/")
    return predictor


def run_inference_validation():
    """
    Boot ml_executor (which loads the PRODUCTION checkpoints) and trigger one
    random generation + one targeted optimisation.  This populates:
      • evaluation/logs/optimization_log.csv
      • evaluation/logs/validity_stats.json
    """
    print(f"\n{'='*60}")
    print("  STAGE 3 — Inference Validation")
    print(f"{'='*60}")

    try:
        from backend.ml_executor import run_random_generation, run_targeted_generation
    except Exception as e:
        print(f"  ⚠️  Could not import ml_executor: {e}")
        print("       Skipping inference validation.")
        return

    # Random generation
    print("  Generating 5 random molecules...")
    try:
        results = run_random_generation(5)
        print(f"  ✅ Random generation: {len(results)} valid molecules returned")
    except Exception as e:
        print(f"  ⚠️  Random generation failed: {e}")

    # Targeted optimization
    print("  Running targeted optimisation (QED=0.8, LogP=2.0, SAS=2.0)...")
    try:
        results = run_targeted_generation(
            num_molecules=3,
            target_qed=0.8,
            target_logp=2.0,
            target_sas=2.0,
        )
        print(f"  ✅ Targeted optimisation: {len(results)} molecules returned")
    except Exception as e:
        print(f"  ⚠️  Targeted optimisation failed: {e}")


def main():
    print(f"\n{'#'*60}")
    print("  SmartChem — Evaluation Training Run")
    print(f"  Device : {DEVICE}")
    print(f"  VAE epochs : {VAE_EPOCHS}   Predictor epochs : {PRED_EPOCHS}")
    print(f"{'#'*60}")

    # Load data
    vocab_path = os.path.join(PROCESSED_DIR, f"{MODE}_vocab.json")
    data_path  = os.path.join(PROCESSED_DIR, f"train_{MODE}.pt")

    with open(vocab_path) as fh:
        vocab = json.load(fh)
    vocab_size = len(vocab) + 3       # +3 for <pad>, <sos>, <eos>

    data    = torch.load(data_path, weights_only=True)
    dataset = TensorDataset(data)
    loader  = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True,
                         drop_last=True, num_workers=0)

    print(f"\n  Dataset  : {len(data)} molecules × {data.shape[1]} tokens")
    print(f"  Vocab    : {vocab_size} tokens")
    print(f"  Batches  : {len(loader)} per epoch  (batch_size={BATCH_SIZE})")

    # Stage 1: VAE
    vae = train_vae(vocab_size, loader)

    # Stage 2: Predictor
    full_loader = DataLoader(
        TensorDataset(data), batch_size=BATCH_SIZE, shuffle=False, num_workers=0
    )
    train_predictor(vae, full_loader)

    # Stage 3: Inference
    run_inference_validation()

    # Stage 4: Plots
    print(f"\n{'='*60}")
    print("  STAGE 4 — Regenerating Evaluation Plots")
    print(f"{'='*60}")
    import subprocess
    plot_script = os.path.join(ROOT, "evaluation", "analysis", "plot_metrics.py")
    subprocess.run([sys.executable, plot_script], check=False)

    print(f"\n{'#'*60}")
    print("  ✅ All done!  Check evaluation/logs/ and evaluation/plots/")
    print(f"{'#'*60}\n")


if __name__ == "__main__":
    main()
