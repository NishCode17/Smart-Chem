import torch
import torch.nn as nn
import torch.nn.functional as F


def reconstruction_loss(logits, targets):
    """
    Cross-entropy reconstruction loss.

    Args:
        logits:  (B, L, vocab_size) — raw decoder output
        targets: (B, L)             — integer token indices

    Returns:
        Scalar mean loss over the batch.
    """
    B, L, V = logits.shape
    logits_flat  = logits.view(B * L, V)
    targets_flat = targets.view(B * L)
    return F.cross_entropy(logits_flat, targets_flat, ignore_index=0, reduction="mean")


def kl_divergence(mu, logvar):
    """
    KL divergence against N(0, I).

    L_KL = -0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)

    Returns:
        Scalar mean KL loss over the batch.
    """
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
    return kl.mean()


def vae_loss(logits, targets, mu, logvar, beta=0.1):
    """
    Combined VAE loss: reconstruction + beta * KL.

    Args:
        logits:  (B, L, vocab_size)
        targets: (B, L)
        mu:      (B, latent_dim)
        logvar:  (B, latent_dim)
        beta:    KL weighting factor (default 0.1)

    Returns:
        total_loss, recon_loss, kl_loss (all scalars)
    """
    recon = reconstruction_loss(logits, targets)
    kl    = kl_divergence(mu, logvar)
    total = recon + beta * kl
    return total, recon, kl
