"""
GNN encoder for molecular graph → latent space.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINEConv, global_mean_pool, global_max_pool


class GNNEncoder(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=128, latent_dim=128):
        super(GNNEncoder, self).__init__()

        # Initial projection
        self.node_proj = nn.Linear(node_dim, hidden_dim)

        # Conv layer 1
        mlp1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv1 = GINEConv(mlp1, edge_dim=edge_dim)
        self.bn1   = nn.BatchNorm1d(hidden_dim)

        # Conv layer 2
        mlp2 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv2 = GINEConv(mlp2, edge_dim=edge_dim)
        self.bn2   = nn.BatchNorm1d(hidden_dim)

        # Conv layer 3
        mlp3 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv3 = GINEConv(mlp3, edge_dim=edge_dim)
        self.bn3   = nn.BatchNorm1d(hidden_dim)

        # Latent heads
        pool_dim = hidden_dim * 2

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
        # Init weights
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

        last_logvar_layer = self.logvar_net[-1]
        nn.init.zeros_(last_logvar_layer.weight)
        nn.init.zeros_(last_logvar_layer.bias)

    def forward(self, x, edge_index, edge_attr, batch):
        # Initial projection
        h = self.node_proj(x)

        # Conv block 1
        h_in = h
        h    = self.conv1(h, edge_index, edge_attr)
        h    = self.bn1(h)
        h    = F.relu(h)
        h    = h + h_in

        # Conv block 2
        h_in = h
        h    = self.conv2(h, edge_index, edge_attr)
        h    = self.bn2(h)
        h    = F.relu(h)
        h    = h + h_in

        # Conv block 3
        h_in = h
        h    = self.conv3(h, edge_index, edge_attr)
        h    = self.bn3(h)
        h    = F.relu(h)
        h    = h + h_in

        # Graph pooling
        h_mean  = global_mean_pool(h, batch)
        h_max   = global_max_pool(h, batch)
        h_graph = torch.cat([h_mean, h_max], dim=1)   # (B, hidden_dim * 2)

        # Latent projections
        mu     = self.mu_net(h_graph)
        logvar = self.logvar_net(h_graph).clamp(-6.0, 2.0)

        return mu, logvar
