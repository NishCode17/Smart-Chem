import torch

def optimize_latent_vector(z, predictor, target_props, steps=75, lr=0.02): 
    """
    Lead Optimization Mode.
    Wider Clamp [-3.0, 3.0] allows existing molecules to survive optimization
    without being crushed, while still preventing explosion.
    """
    z_opt = z.clone().detach().requires_grad_(True)
    optimizer = torch.optim.Adam([z_opt], lr=lr)
    
    target_tensor = torch.tensor(target_props).float().to(z.device)
    # Weights: [QED=10, LogP=1, SAS=1]
    weights = torch.tensor([10.0, 1.0, 1.0]).to(z.device)
    
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
        
    return z_opt.detach()