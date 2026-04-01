import torch

def optimize_latent_vector(z, predictor, target_props, steps=75, lr=0.02,
                           eval_log=False):
    """
    Lead Optimization Mode.
    Wider Clamp [-3.0, 3.0] allows existing molecules to survive optimization
    without being crushed, while still preventing explosion.

    Parameters
    ----------
    z            : Initial latent tensor (batch_size × latent_dim)
    predictor    : PropertyPredictor model
    target_props : [target_qed, target_logp, target_sas]
    steps        : Number of gradient-descent steps (default 75)
    lr           : Adam learning rate (default 0.02)
    eval_log     : If True, emit per-step rows to evaluation/logs/optimization_log.csv
    """
    z_opt = z.clone().detach().requires_grad_(True)
    optimizer = torch.optim.Adam([z_opt], lr=lr)

    target_tensor = torch.tensor(target_props).float().to(z.device)
    # Weights: [QED=10, LogP=1, SAS=1]
    weights = torch.tensor([10.0, 1.0, 1.0]).to(z.device)

    # Lazy import so that optimizer.py stays independent of eval infrastructure
    if eval_log:
        try:
            from evaluation.eval_logger import log_optimization_step
        except ImportError:
            eval_log = False

    for i in range(steps):
        optimizer.zero_grad()
        preds = predictor(z_opt)

        # Target Loss
        diff_sq = (preds - target_tensor.repeat(z.shape[0], 1)) ** 2
        weighted_dist = (diff_sq * weights).sum(dim=1).mean()

        # Anchor Loss (Stay close to the lead molecule!)
        dist_from_seed = torch.norm(z_opt - z, p=2)
        anchor_penalty = 0.5 * dist_from_seed

        loss = weighted_dist + anchor_penalty
        loss.backward()
        optimizer.step()

        # WIDENED CLAMP (The Fix)
        # [-3, 3] covers 99.7% of the latent space.
        # 1.5 was too tight and destroying valid leads.
        with torch.no_grad():
            z_opt.data.clamp_(-3.0, 3.0)

        # ── Evaluation logging (non-blocking, best-effort) ──────────────────
        if eval_log:
            with torch.no_grad():
                mean_preds = predictor(z_opt).mean(dim=0)   # shape [3]
                mean_qed   = mean_preds[0].item()
                mean_logp  = mean_preds[1].item()
                l2_dist    = torch.norm(z_opt - z, p=2).item() / z.shape[0]
            try:
                log_optimization_step(
                    step=i,
                    qed=mean_qed,
                    logp=mean_logp,
                    l2_distance=l2_dist,
                )
            except Exception as _log_err:
                pass  # Never let logging crash the optimization

    return z_opt.detach()