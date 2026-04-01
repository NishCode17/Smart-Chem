"""
evaluation/training_hooks.py
=============================
Copy-paste snippets showing exactly WHICH lines to add to your VAE and
Predictor training scripts to activate evaluation logging.

This file is documentation-as-code — it is NOT imported by any other module.
Run it on its own to verify the logger API is importable:

    python evaluation/training_hooks.py
"""

# ── Verify imports ────────────────────────────────────────────────────────────
from evaluation.eval_logger import (
    log_vae_epoch,
    log_predictor_epoch,
    log_predictor_sample,
    log_optimization_step,
    log_validity_stats,
)

print("✅ eval_logger API imported successfully.\n")

# ─────────────────────────────────────────────────────────────────────────────
# HOW TO INSTRUMENT  vae_trainer.py
# ─────────────────────────────────────────────────────────────────────────────
VAE_TRAINING_SNIPPET = """
# --- ADD AT TOP of vae_trainer.py ---
from evaluation.eval_logger import log_vae_epoch

# --- INSIDE epoch loop ---
for epoch in range(1, NUM_EPOCHS + 1):
    for batch in dataloader:
        optimizer.zero_grad()
        recon, mu, logvar = model(batch)

        bce  = F.binary_cross_entropy_with_logits(recon, target, reduction='sum')
        kl   = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        loss = bce + beta * kl
        loss.backward()
        optimizer.step()

    # ← ADD THIS LINE (one call per epoch)
    log_vae_epoch(epoch=epoch,
                  bce_loss=bce.item() / len(dataloader.dataset),
                  kl_loss=kl.item()  / len(dataloader.dataset),
                  total_loss=loss.item() / len(dataloader.dataset))
"""

# ─────────────────────────────────────────────────────────────────────────────
# HOW TO INSTRUMENT  predictor_trainer.py
# ─────────────────────────────────────────────────────────────────────────────
PREDICTOR_TRAINING_SNIPPET = """
# --- ADD AT TOP of predictor_trainer.py ---
from evaluation.eval_logger import log_predictor_epoch, log_predictor_sample

# --- INSIDE epoch loop ---
for epoch in range(1, NUM_EPOCHS + 1):
    # Training pass
    for z_batch, y_batch in dataloader:
        preds = predictor(z_batch)
        loss_qed  = F.mse_loss(preds[:, 0], y_batch[:, 0])
        loss_logp = F.mse_loss(preds[:, 1], y_batch[:, 1])
        loss_sas  = F.mse_loss(preds[:, 2], y_batch[:, 2])
        loss      = loss_qed + loss_logp + loss_sas
        loss.backward()
        optimizer.step()

        # ← Log per-sample true/pred for scatter plot
        for i in range(len(z_batch)):
            log_predictor_sample(
                true_qed=y_batch[i,0].item(), pred_qed=preds[i,0].item(),
                true_logp=y_batch[i,1].item(), pred_logp=preds[i,1].item(),
                true_sas=y_batch[i,2].item(),  pred_sas=preds[i,2].item(),
            )

    # ← ADD THIS (one call per epoch)
    log_predictor_epoch(epoch=epoch,
                        mse_qed=loss_qed.item(),
                        mse_logp=loss_logp.item(),
                        mse_sas=loss_sas.item())
"""

print("VAE training snippet:")
print(VAE_TRAINING_SNIPPET)
print("\nPredictor training snippet:")
print(PREDICTOR_TRAINING_SNIPPET)
print("\nOptimisation and validity logging are already wired into:")
print("  • backend/optimizer.py      (eval_log=True flag)")
print("  • backend/ml_executor.py    (auto-logged on every inference call)")
