# SmartChem — Evaluation Pipeline

This directory contains a **fully self-contained** evaluation, logging, and
visualisation pipeline for the SmartChem ML system.  
It is intentionally kept **decoupled** from core model code so it can be run,
modified, or removed without disrupting training or inference.

---

## Directory Layout

```
evaluation/
├── eval_logger.py          ← Central logging API (import this in training scripts)
├── logs/
│   ├── vae_training_log.csv        ← Per-epoch VAE loss (auto-written during training)
│   ├── predictor_log.csv           ← Per-epoch predictor MSE (auto-written during training)
│   ├── predictor_true_pred_log.csv ← Per-sample true/pred values (auto-written)
│   ├── optimization_log.csv        ← Per-step latent optimisation metrics (auto-written during inference)
│   └── validity_stats.json         ← Molecule generation quality summary (auto-written during inference)
├── plots/
│   ├── vae_loss.png                ← VAE BCE / KL / Total loss vs epoch
│   ├── predictor_scatter.png       ← True vs Predicted QED scatter plot
│   ├── predictor_training_loss.png ← MSE per property head vs epoch
│   ├── optimization_curve.png      ← QED & L2 trajectory over optimisation steps
│   └── validity_stats.png          ← Bar chart of generation / validity statistics
└── analysis/
    └── plot_metrics.py             ← Standalone script that generates all plots
```

---

## Log File Reference

| File | Written by | Key columns |
|---|---|---|
| `vae_training_log.csv` | VAE training loop | `epoch`, `bce_loss`, `kl_loss`, `total_loss` |
| `predictor_log.csv` | Predictor training loop | `epoch`, `mse_qed`, `mse_logp`, `mse_sas` |
| `predictor_true_pred_log.csv` | Predictor training loop | `true_qed/logp/sas`, `pred_qed/logp/sas` |
| `optimization_log.csv` | `backend/optimizer.py` (via `eval_log=True`) | `step`, `qed`, `logp`, `l2_distance` |
| `validity_stats.json` | `backend/ml_executor.py` | `total_generated`, `valid_selfies`, `passed_rdkit`, `passed_lipinski` |

---

## Running the Analysis Script

> **Requirements:** `pandas`, `matplotlib` (already in `requirements.txt`)

From the **project root**:

```bash
python evaluation/analysis/plot_metrics.py
```

The script:
1. Reads each log file (falls back to **simulated data** gracefully if a log is absent).
2. Generates all five plots with a dark-mode, publication-ready aesthetic.
3. Saves them to `evaluation/plots/`.

No Model weights, RDKit, or PyTorch are needed to run this script.

---

## Plot Descriptions

### `vae_loss.png`
Three loss curves (BCE, KL divergence, Total) plotted against training epoch.  
Shows the balance between reconstruction fidelity and latent-space regularisation.  
**Healthy training**: BCE decreases monotonically; KL rises then stabilises.

### `predictor_scatter.png`
Scatter plot of true vs predicted QED values.  
The diagonal line (`y = x`) represents perfect prediction.  
R² is annotated on the chart.  
**Good predictor**: points cluster tightly around the diagonal.

### `predictor_training_loss.png`
MSE loss curves for each of the three property heads (QED, LogP, SAS) across epochs.  
Useful for diagnosing which property the predictor finds hardest to learn.

### `optimization_curve.png`
Two-panel chart:
- **Top**: Mean predicted QED over gradient-descent steps — shows property improvement.
- **Bottom**: Mean L2 distance of the optimised latent vector from the seed — quantifies exploration radius.

### `validity_stats.png`
Bar chart with four categories:
1. **Generated** — total decoded SELFIES sequences
2. **Valid (SELFIES)** — sequences that parsed to a valid RDKit molecule
3. **Passed RDKit** — mols that additionally passed structural filters (`_is_valid_candidate`)
4. **Passed Lipinski** — mols passing Rule-of-5 (≤ 1 violation)

Percentage rates are annotated inside each bar.

---

## Integration with Training Scripts

Call the logger functions at appropriate points in your training loops:

```python
from evaluation.eval_logger import (
    log_vae_epoch,
    log_predictor_epoch,
    log_predictor_sample,
)

# Inside VAE training loop:
log_vae_epoch(epoch=epoch_num, bce_loss=bce.item(),
              kl_loss=kl.item(), total_loss=loss.item())

# Inside predictor training loop:
log_predictor_epoch(epoch=epoch_num,
                    mse_qed=mse_q.item(),
                    mse_logp=mse_l.item(),
                    mse_sas=mse_s.item())

# Per sample during predictor validation:
log_predictor_sample(true_qed=y[i,0].item(), pred_qed=yhat[i,0].item(),
                     true_logp=y[i,1].item(), pred_logp=yhat[i,1].item(),
                     true_sas=y[i,2].item(),  pred_sas=yhat[i,2].item())
```

Optimisation logs and validity stats are written **automatically** by `backend/optimizer.py`
and `backend/ml_executor.py` whenever the inference pipeline runs.

---

## Academic Justification

| Metric | Why it matters |
|---|---|
| BCE + KL loss | Quantitative evidence that the VAE learns a meaningful latent space |
| Predictor R² | Validates that property prediction is reliable enough to guide optimisation |
| Optimisation trajectory | Demonstrates that gradient ascent in latent space meaningfully improves QED |
| Validity rates | Reports generation quality — critical for drug-discovery benchmarks |

---

## Reproducibility

- All CSV logs are append-only and timestamped so runs accumulate cleanly.
- `validity_stats.json` is **overwritten** on each inference run (latest snapshot).
- The analysis script uses a fixed random seed (`numpy.random.seed(42 / 7)`) for any
  simulated data, ensuring deterministic placeholders.
