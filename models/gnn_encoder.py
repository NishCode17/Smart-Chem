"""
models/gnn_encoder.py

GNN encoder for molecular graph → latent space.

Architecture (3-layer GINEConv, matching original spec):
  - Initial node projection → hidden_dim
  - 3 × GINEConv blocks with BatchNorm + residuals
  - Global mean + max pooling (concatenated)
  - Two-layer MLP heads → mu, logvar
  - logvar clamped to [-6, 2]  (matches CNN encoder — prevents NaN / extreme variance)
  - logvar_net zero-initialised → starts near unit Gaussian prior

Fixes vs. original (large-dataset robustness):
  1. 3rd GINEConv layer added  — matches spec; deeper receptive field for 100k diversity
  2. logvar clamped to [-6, 2] — at 100k with more gradient steps, unclamped logvar
                                  can collapse or explode (exp(logvar) → 0 or inf)
  3. logvar_net last layer zero-init — mirrors CNN encoder; stable start at large scale
  4. Removed F.normalize before MLP — hard unit-sphere constraint kills gradient
                                       diversity when minibatch distribution is wide
                                       (100k); replaced with LayerNorm inside the MLP
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINEConv, global_mean_pool, global_max_pool


class GNNEncoder(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=128, latent_dim=128):
        super(GNNEncoder, self).__init__()

        # Initial node projection to match hidden_dim for residual connections
        self.node_proj = nn.Linear(node_dim, hidden_dim)

        # ── GINEConv layer 1 ──────────────────────────────────────────────
        mlp1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv1 = GINEConv(mlp1, edge_dim=edge_dim)
        self.bn1   = nn.BatchNorm1d(hidden_dim)

        # ── GINEConv layer 2 ──────────────────────────────────────────────
        mlp2 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv2 = GINEConv(mlp2, edge_dim=edge_dim)
        self.bn2   = nn.BatchNorm1d(hidden_dim)

        # ── GINEConv layer 3 (NEW — matches 3-layer spec) ─────────────────
        mlp3 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv3 = GINEConv(mlp3, edge_dim=edge_dim)
        self.bn3   = nn.BatchNorm1d(hidden_dim)

        # ── Latent heads ──────────────────────────────────────────────────
        # FIX: LayerNorm inside MLP replaces pre-MLP F.normalize.
        # F.normalize hard-constrains all inputs to the unit sphere; at large
        # dataset scale the wide minibatch distribution fights this constraint,
        # slowing convergence and hurting latent diversity.
        pool_dim = hidden_dim * 2   # mean + max pool concatenated

        self.mu_net = nn.Sequential(
            nn.Linear(pool_dim, 128),
            nn.LayerNorm(128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )
        self.logvar_net = nn.Sequential(
            nn.Linear(pool_dim, 128),
            nn.LayerNorm(128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

        self._init_weights()

    def _init_weights(self):
        """Xavier init for all linear layers; zero-init last logvar layer → unit Gaussian start."""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

        # FIX: zero-init logvar output layer — mirrors CNN encoder.
        # Ensures logvar ≈ 0 at training start (unit Gaussian prior),
        # preventing early KL spikes that can lock the model into collapse.
        last_logvar_layer = self.logvar_net[-1]   # nn.Linear(128, latent_dim)
        nn.init.zeros_(last_logvar_layer.weight)
        nn.init.zeros_(last_logvar_layer.bias)

    def forward(self, x, edge_index, edge_attr, batch):
        # ── Initial projection ────────────────────────────────────────────
        h = self.node_proj(x)

        # ── Conv block 1 with residual ────────────────────────────────────
        h_in = h
        h    = self.conv1(h, edge_index, edge_attr)
        h    = self.bn1(h)
        h    = F.relu(h)
        h    = h + h_in

        # ── Conv block 2 with residual ────────────────────────────────────
        h_in = h
        h    = self.conv2(h, edge_index, edge_attr)
        h    = self.bn2(h)
        h    = F.relu(h)
        h    = h + h_in

        # ── Conv block 3 with residual (NEW) ──────────────────────────────
        h_in = h
        h    = self.conv3(h, edge_index, edge_attr)
        h    = self.bn3(h)
        h    = F.relu(h)
        h    = h + h_in

        # ── Graph pooling (mean + max) ────────────────────────────────────
        h_mean  = global_mean_pool(h, batch)
        h_max   = global_max_pool(h, batch)
        h_graph = torch.cat([h_mean, h_max], dim=1)   # (B, hidden_dim * 2)

        # NOTE: F.normalize removed here — see module docstring for rationale.
        # LayerNorm is applied inside mu_net / logvar_net instead.

        # ── Latent projections ────────────────────────────────────────────
        mu     = self.mu_net(h_graph)
        logvar = self.logvar_net(h_graph).clamp(-6.0, 2.0)   # FIX: numerical guard

        return mu, logvar
