import torch
import torch.nn as nn
import torch.nn.functional as F


class HybridEncoder(nn.Module):
    def __init__(self, cnn_encoder, gnn_encoder, latent_dim=128):
        super(HybridEncoder, self).__init__()

        self.cnn_encoder = cnn_encoder
        self.gnn_encoder = gnn_encoder

        # Freeze base encoders
        for p in self.cnn_encoder.parameters():
            p.requires_grad = False
        for p in self.gnn_encoder.parameters():
            p.requires_grad = False

        # Projection layers
        self.seq_proj   = nn.Linear(latent_dim, 64)
        self.graph_proj = nn.Linear(latent_dim, 64)

        # Fusion gate
        self.gate_layer = nn.Linear(128, 128)

        # Fusion MLP
        self.fusion = nn.Sequential(
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

        # Latent heads
        self.mu_layer     = nn.Linear(latent_dim, latent_dim)
        self.logvar_layer = nn.Linear(latent_dim, latent_dim)

        # Cached features
        self._h_seq_raw   = None
        self._h_graph_raw = None

    def forward(self, seq_input, graph_data):
        # Extract features
        with torch.no_grad():
            mu_seq, _   = self.cnn_encoder(seq_input)
            mu_graph, _ = self.gnn_encoder(
                graph_data.x,
                graph_data.edge_index,
                graph_data.edge_attr,
                graph_data.batch,
            )

        # Project
        h_seq   = F.relu(self.seq_proj(mu_seq))       # (B, 64)
        h_graph = F.relu(self.graph_proj(mu_graph))   # (B, 64)

        # Gated fusion
        gate           = torch.sigmoid(self.gate_layer(torch.cat([h_seq, h_graph], dim=1)))
        g_seq, g_graph = gate.chunk(2, dim=1)          # each (B, 64)

        h_seq_gated   = g_seq   * h_seq
        h_graph_gated = g_graph * h_graph

        # Store pre-norm
        self._h_seq_raw   = h_seq_gated
        self._h_graph_raw = h_graph_gated

        # Normalize
        h_seq_n   = F.normalize(h_seq_gated,   dim=-1)
        h_graph_n = F.normalize(h_graph_gated, dim=-1)

        # Concat + Fusion
        h_cat   = torch.cat([h_seq_n, h_graph_n], dim=1)
        h_fused = self.fusion(h_cat)

        # Latent heads
        mu     = self.mu_layer(h_fused)
        logvar = self.logvar_layer(h_fused)

        return mu, logvar
