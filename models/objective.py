import torch

def objective(z, predictor):
    preds = predictor(z)  # (B, 3)
    qed = preds[:, 0]
    sas = preds[:, 2]
    return qed - 0.1 * sas