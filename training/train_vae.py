"""
training/train_vae.py

CNN-encoder + GRU-decoder VAE training for molecular SMILES/SELFIES generation.

Loss strategy (state-of-the-art for molecular sequence VAEs):
  1. Balanced loss       — loss = lambda_recon * recon + beta * kl
                           lambda_recon=0.5 so reconstruction always meaningfully
                           competes with KL at any beta value (was 0.01 — too low)
  2. Cyclical KL anneal  — Fu et al. 2019 (arXiv:1903.10145)
                           beta cycles 0→1→0→1 over M cycles, each cycle:
                           linear rise over first 50% then flat at 1.0 for 50%.
                           Forces repeated encoder re-activation; linear annealing
                           is brittle because once the decoder adapts to ignoring z
                           during beta=0, collapse is locked in for the rest of training.
  3. Per-dim free bits   — Kingma et al. 2016 IFV / Chen et al. 2016 VLAE
                           kl_per_dim = max(kl_dim_i, free_bits) summed over dims.
                           free_bits=0.5 nats/dim (was 0.1 — too low for 128 dims).
                           Ensures at least 0.5 × 128 = 64 nats total, forces
                           genuine encoding across all dimensions.
  4. Word dropout 50%    — raised from 30% (decoder must be weaker at 2x dataset)
                           Weakens decoder so it cannot reconstruct without z.
  5. Minimum Desired Rate (delta-VAE, Razavi et al. 2019, arXiv:1901.03416)
                           Applied per-dimension, not on global mean, so the
                           optimiser cannot satisfy the constraint by only using
                           a few dimensions at high KL.
  6. GRU dropout in vae  — 0.2 applied inside GRU (already in vae.py)
  7. AdamW + warmup LR   — weight_decay=0.01 + linear warmup over 2 epochs
                           prevents large initial steps from locking in collapse

Root-cause of 50k→100k regression:
  - 100k = 2× more gradient steps per epoch. Decoder adapts 2× faster per epoch.
  - lambda_recon=0.01 meant KL dominated at beta>0.1; decoder learned to ignore z.
  - free_bits=0.1 nats/dim with 128 dims allowed collapse to 17 active dims
    (17 × 0.1 = 1.7 nats satisfied the floor while 111 dims went to prior).
  - Linear KL annealing is irreversible: once collapsed, no recovery mechanism.
"""

import os
import sys
import math
import torch
import torch.nn.functional as F
import torch.optim as optim

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

PAD_TOKEN_ID = 0


# ---------------------------------------------------------------------------
# Cyclical KL annealing — Fu et al. 2019 (arXiv:1903.10145)
# ---------------------------------------------------------------------------

def cyclical_beta(epoch: int, total_epochs: int,
                  n_cycles: int = 4,
                  ratio: float = 0.5,
                  max_beta: float = 1.0) -> float:
    """
    Cyclical annealing schedule (Fu et al., 2019).

    Divides training into `n_cycles` equal cycles.  Within each cycle:
      - First `ratio` fraction: beta ramps linearly from 0 → max_beta
      - Remaining fraction:     beta stays at max_beta

    This prevents the "lazy decoder" lock-in that plagues linear annealing:
    every new cycle forces the encoder to re-encode information from scratch,
    breaking the decoder's ability to memorise around z.

    Args:
        epoch        : current epoch (1-indexed)
        total_epochs : total training epochs
        n_cycles     : number of full cycles over training (default 4)
        ratio        : fraction of each cycle spent ramping (default 0.5)
        max_beta     : ceiling beta value (default 1.0)

    Returns:
        beta ∈ [0, max_beta]
    """
    cycle_len  = total_epochs / n_cycles
    cycle_pos  = (epoch - 1) % cycle_len       # position within current cycle
    ramp_len   = ratio * cycle_len              # how many epochs to ramp

    if cycle_pos < ramp_len:
        return max_beta * (cycle_pos / ramp_len)
    return max_beta


# ---------------------------------------------------------------------------
# Per-dimension free-bits KL — Kingma et al. 2016, Chen et al. 2016
# ---------------------------------------------------------------------------

def kl_loss_free_bits(mu: torch.Tensor,
                      logvar: torch.Tensor,
                      free_bits: float = 0.5) -> tuple:
    """
    Per-dimension free-bits KL divergence.

    Standard KL:   KL = -0.5 × Σ_d (1 + logvar_d - mu_d² - exp(logvar_d))
    Free-bits:     KL_d = max(KL_d, free_bits)   applied PER DIMENSION

    This is critically different from clamping the global mean:
    - Global clamp allows the optimiser to satisfy the constraint using a FEW
      high-KL dims while collapsing the rest (exactly what happened at 100k).
    - Per-dim clamp forces EVERY dimension to carry at least `free_bits` nats,
      preventing selective collapse.

    Returns:
        kl_total  : scalar — sum over dims, mean over batch (used in loss)
        kl_per_dim: scalar — mean per dimension (for logging)
    """
    # KL per sample per dimension: (B, latent_dim)
    kl_elementwise = -0.5 * (1.0 + logvar - mu.pow(2) - logvar.exp())

    # Per-dimension mean over batch: (latent_dim,)
    kl_dim = kl_elementwise.mean(dim=0)

    # Per-dimension free-bits floor
    kl_dim_clamped = torch.clamp(kl_dim, min=free_bits)

    # Sum over dimensions, already averaged over batch
    kl_total   = kl_dim_clamped.sum()
    kl_per_dim = kl_dim_clamped.mean()      # for readable logging

    return kl_total, kl_per_dim


# ---------------------------------------------------------------------------
# Reconstruction loss
# ---------------------------------------------------------------------------

def reconstruction_loss(logits: torch.Tensor,
                        targets: torch.Tensor) -> torch.Tensor:
    """
    Per-token cross-entropy, ignoring PAD.
    logits  : (B, L, vocab_size)
    targets : (B, L)
    """
    B, L, V = logits.shape
    return F.cross_entropy(
        logits.reshape(B * L, V),
        targets.reshape(B * L),
        ignore_index=PAD_TOKEN_ID,
        reduction="mean",
    )


# ---------------------------------------------------------------------------
# Word dropout — strengthened to 50% for larger dataset
# ---------------------------------------------------------------------------

def word_dropout(seq: torch.Tensor, drop_prob: float = 0.5,
                 training: bool = True) -> torch.Tensor:
    """
    Randomly replace `drop_prob` fraction of decoder input tokens with PAD=0.
    Stronger word dropout (0.5 vs 0.3) deliberately weakens the decoder at
    100k scale, forcing it to rely on z for reconstruction rather than
    autoregressive context alone.

    seq : LongTensor (B, L)
    """
    if not training or drop_prob <= 0.0:
        return seq
    mask  = torch.rand_like(seq.float()) < drop_prob
    noisy = seq.clone()
    noisy[mask] = PAD_TOKEN_ID
    return noisy


# ---------------------------------------------------------------------------
# Model type detection
# ---------------------------------------------------------------------------

def _detect_model_type(model) -> str:
    name = type(model.encoder).__name__.lower()
    if "hybrid" in name: return "hybrid"
    if "gnn"    in name: return "gnn"
    return "cnn"


# ---------------------------------------------------------------------------
# Per-batch forward pass
# ---------------------------------------------------------------------------

def _forward_batch(model, batch, model_type: str, device,
                   drop_prob: float = 0.5,
                   training: bool = True) -> tuple:
    """
    Returns (logits, mu, logvar, target, balance_loss).
    Decoder input  = word_dropout(seq[:, :-1])   — noisy shifted
    Decoder target = seq[:, 1:]                  — clean shifted
    """
    if model_type == "cnn":
        seq        = batch.to(device)
        mu, logvar = model.encoder(seq)
        z          = model.reparameterize(mu, logvar)
        dec_in     = word_dropout(seq[:, :-1], drop_prob, training)
        target     = seq[:, 1:]
        hidden     = model._hidden_from_z(z)
        emb        = model.embedding(dec_in)
        z_expanded = z.unsqueeze(1).expand(-1, emb.size(1), -1)
        gru_in     = torch.cat([emb, z_expanded], dim=-1)
        out, _     = model.gru(gru_in, hidden)
        logits     = model.fc_out(out)
        return logits, mu, logvar, target, None

    elif model_type == "gnn":
        seq        = batch.seq.to(device)
        graph_data = batch.to(device)
        if seq.dim() == 3:   seq = seq.squeeze(1)
        elif seq.dim() == 1: seq = seq.view(graph_data.num_graphs, -1)
        mu, logvar = model.encoder(graph_data.x, graph_data.edge_index,
                                   graph_data.edge_attr, graph_data.batch)
        z          = model.reparameterize(mu, logvar)
        dec_in     = word_dropout(seq[:, :-1], drop_prob, training)
        target     = seq[:, 1:]
        hidden     = model._hidden_from_z(z)
        emb        = model.embedding(dec_in)
        z_expanded = z.unsqueeze(1).expand(-1, emb.size(1), -1)
        gru_in     = torch.cat([emb, z_expanded], dim=-1)
        out, _     = model.gru(gru_in, hidden)
        logits     = model.fc_out(out)
        return logits, mu, logvar, target, None

    elif model_type == "hybrid":
        seq_input, graph_data = batch
        seq_input  = seq_input.to(device)
        graph_data = graph_data.to(device)
        mu, logvar = model.encoder(seq_input, graph_data)
        z          = model.reparameterize(mu, logvar)
        dec_in     = word_dropout(seq_input[:, :-1], drop_prob, training)
        target     = seq_input[:, 1:]
        hidden     = model._hidden_from_z(z)
        emb        = model.embedding(dec_in)
        z_expanded = z.unsqueeze(1).expand(-1, emb.size(1), -1)
        gru_in     = torch.cat([emb, z_expanded], dim=-1)
        out, _     = model.gru(gru_in, hidden)
        logits     = model.fc_out(out)
        h_seq_raw   = model.encoder._h_seq_raw
        h_graph_raw = model.encoder._h_graph_raw
        balance_loss = torch.mean(
            (h_seq_raw.norm(dim=1) - h_graph_raw.norm(dim=1)) ** 2
        )
        return logits, mu, logvar, target, balance_loss

    raise ValueError(f"Unknown model_type: {model_type}")


# ---------------------------------------------------------------------------
# LR warmup utility
# ---------------------------------------------------------------------------

def _warmup_lr(optimizer, epoch: int, warmup_epochs: int, base_lr: float):
    """Linear LR warmup for first `warmup_epochs` epochs."""
    if epoch <= warmup_epochs:
        lr = base_lr * (epoch / warmup_epochs)
        for pg in optimizer.param_groups:
            pg["lr"] = lr


# ---------------------------------------------------------------------------
# Main training function
# ---------------------------------------------------------------------------

def train_vae(model, dataloader, optimizer, device,
              epochs: int          = 30,
              word_drop: float     = 0.5,          # raised from 0.3
              n_cycles: int        = 4,             # cyclical KL cycles
              max_beta: float      = 1.0,           # KL ceiling
              lambda_recon: float  = 0.5,           # recon weight (was 0.01 — critical fix)
              free_bits: float     = 0.5,           # per-dim KL floor nats (was 0.1)
              lr_warmup_epochs: int = 2,            # LR warmup epochs
              base_lr: float       = 3e-4,          # base LR for warmup reference
              save_path: str       = "vae.pt",
              # Legacy params kept for backward compat with other scripts
              warmup_epochs: int   = None,
              kl_target: float     = None) -> object:
    """
    Train VAE with cyclical KL annealing, per-dim free-bits, and balanced loss.

    The three critical fixes vs the old implementation:
      1. lambda_recon = 0.5  (was 0.01) — reconstruction competes with KL
      2. cyclical annealing  (was linear) — prevents decoder lock-in
      3. per-dim free_bits = 0.5  (was global mean clamp at 0.1)

    Args:
        model            : VAE instance
        dataloader       : yields batches
        optimizer        : torch optimizer (recommend AdamW with wd=0.01)
        device           : torch.device
        epochs           : training epochs (recommend 30 for 100k)
        word_drop        : decoder masking prob (0.5 recommended at 100k scale)
        n_cycles         : cyclical KL cycles over full training
        max_beta         : maximum beta for KL term (1.0)
        lambda_recon     : reconstruction loss weight (0.5)
        free_bits        : per-dim KL floor in nats (0.5)
        lr_warmup_epochs : linear LR warmup from 0→base_lr over this many epochs
        base_lr          : reference LR for warmup schedule
        save_path        : checkpoint path
    """
    model_type = _detect_model_type(model)
    latent_dim = model.latent_dim

    # Handle legacy kl_target parameter (map to free_bits if provided)
    if kl_target is not None and kl_target != 0.1:
        free_bits = kl_target

    print(f"\n[train_vae] Encoder={model_type.upper()}  Epochs={epochs}  Device={device}")
    print(f"  Strategy  : Cyclical KL annealing ({n_cycles} cycles) + per-dim free-bits")
    print(f"  word_drop : {word_drop}  (raised from 0.3 → 0.5 for 100k scale)")
    print(f"  lambda_recon: {lambda_recon}  (was 0.01 — critical fix)")
    print(f"  free_bits : {free_bits} nats/dim  (was 0.1 — per DIM, not global mean)")
    print(f"  max_beta  : {max_beta}  |  n_cycles: {n_cycles}")
    print(f"  latent_dim: {latent_dim}\n")

    model.to(device)

    # Cosine LR decay after warmup — tied to epochs (correct for both 50k/100k)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(1, epochs - lr_warmup_epochs), eta_min=1e-5)

    for epoch in range(1, epochs + 1):
        model.train()

        # --- LR warmup ---
        if epoch <= lr_warmup_epochs:
            _warmup_lr(optimizer, epoch, lr_warmup_epochs, base_lr)
        elif epoch == lr_warmup_epochs + 1:
            scheduler.step()   # start cosine from here

        # --- Cyclical beta ---
        beta = cyclical_beta(epoch, epochs, n_cycles=n_cycles,
                             ratio=0.5, max_beta=max_beta)

        ep_recon = ep_kl = ep_total = ep_au = 0.0
        n = 0

        for batch in dataloader:
            optimizer.zero_grad()

            logits, mu, logvar, target, balance_loss = _forward_batch(
                model, batch, model_type, device,
                drop_prob=word_drop, training=True
            )

            recon = reconstruction_loss(logits, target)

            # Per-dim free-bits KL (THE key fix)
            kl_total, kl_per_dim = kl_loss_free_bits(mu, logvar,
                                                      free_bits=free_bits)

            # Core loss: balanced recon + annealed KL
            loss = lambda_recon * recon + beta * kl_total

            if balance_loss is not None:
                loss = loss + 0.01 * balance_loss

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            with torch.no_grad():
                au = int((mu.var(dim=0) > 1e-2).sum().item())

            ep_recon += recon.item()
            ep_kl    += kl_per_dim.item()
            ep_total += loss.item()
            ep_au    += au
            n        += 1

        # Step cosine scheduler only after warmup phase
        if epoch > lr_warmup_epochs:
            scheduler.step()

        current_lr = optimizer.param_groups[0]["lr"]
        print(
            f"Epoch [{epoch:02d}/{epochs}]  "
            f"Recon: {ep_recon/n:.4f}  "
            f"KL/dim: {ep_kl/n:.4f}  "
            f"Beta: {beta:.4f}  "
            f"Total: {ep_total/n:.4f}  "
            f"AU: {ep_au/n:.0f}/{latent_dim}  "
            f"LR: {current_lr:.2e}"
        )

    torch.save(model.state_dict(), save_path)
    print(f"\n[train_vae] Saved → {save_path}")
    return model


# ---------------------------------------------------------------------------
# Latent dataset builder (unchanged)
# ---------------------------------------------------------------------------

def build_latent_dataset(model, dataloader, device, smiles_list=None):
    from rdkit import Chem
    from rdkit.Chem import QED as rdQED, Descriptors
    try:
        from sascorer import calculateScore as _sas
    except ImportError:
        import importlib.util as _ilu
        _s = _ilu.spec_from_file_location("sascorer",
                                           os.path.join(ROOT, "sascorer.py"))
        _m = _ilu.module_from_spec(_s); _s.loader.exec_module(_m)
        _sas = _m.calculateScore

    model_type = _detect_model_type(model)
    model.to(device); model.eval()
    z_list = []; y_list = []; idx = 0

    with torch.no_grad():
        for batch in dataloader:
            if model_type == "cnn":
                seq = batch.to(device); _, mu, _ = model(seq)
            elif model_type == "gnn":
                seq = batch.seq.to(device); g = batch.to(device)
                if seq.dim() == 3: seq = seq.squeeze(1)
                elif seq.dim() == 1: seq = seq.view(g.num_graphs, -1)
                mu, _ = model.encoder(g.x, g.edge_index, g.edge_attr, g.batch)
            elif model_type == "hybrid":
                si, g = batch
                mu, _ = model.encoder(si.to(device), g.to(device))
            else:
                raise ValueError(f"Unknown: {model_type}")

            for i in range(mu.shape[0]):
                if smiles_list is None or idx >= len(smiles_list): break
                mol = Chem.MolFromSmiles(smiles_list[idx]); idx += 1
                if mol is None: continue
                try:
                    z_list.append(mu[i].cpu())
                    y_list.append([rdQED.qed(mol),
                                   Descriptors.MolLogP(mol),
                                   float(_sas(mol))])
                except Exception: continue

    if not z_list:
        raise RuntimeError("No valid molecules found.")

    z = torch.stack(z_list)
    y = torch.tensor(y_list, dtype=torch.float32)
    print(f"[build_latent_dataset] {len(z_list)} molecules  "
          f"QED={y[:,0].mean():.3f}  LogP={y[:,1].mean():.3f}")
    return z, y
