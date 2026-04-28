"""
models/cnn_encoder.py

Upgraded CNN encoder for molecular sequence → latent space.

Architecture:
  - Embedding → 4 dilated Conv1d blocks with residual connections
  - Layer Norm + GELU (more stable gradients than ReLU)
  - AdaptiveMaxPool + AdaptiveAvgPool (richer pooling)
  - Two-layer MLP projection → mu, logvar
  - logvar clamped to [-6, 2] (prevents NaN / extreme variance)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class _ConvBlock(nn.Module):
    """Dilated Conv1d → LayerNorm → GELU with residual."""
    def __init__(self, in_ch, out_ch, kernel_size, dilation=1):
        super().__init__()
        pad = (kernel_size - 1) * dilation // 2
        self.conv = nn.Conv1d(in_ch, out_ch, kernel_size,
                              padding=pad, dilation=dilation)
        self.norm = nn.LayerNorm(out_ch)
        self.proj = nn.Conv1d(in_ch, out_ch, 1) if in_ch != out_ch else nn.Identity()

    def forward(self, x):
        # x : (B, C, L)
        out = self.conv(x)
        out = self.norm(out.transpose(1, 2)).transpose(1, 2)
        out = F.gelu(out)
        return out + self.proj(x)


class CNNEncoder(nn.Module):
    def __init__(self, vocab_size, embedding_dim=128, latent_dim=128,
                 hidden_dim=256):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embedding_dim,
                                      padding_idx=0)

        # 4 dilated conv blocks with increasing dilation for larger receptive field
        self.blocks = nn.Sequential(
            _ConvBlock(embedding_dim, hidden_dim, kernel_size=3, dilation=1),
            _ConvBlock(hidden_dim,    hidden_dim, kernel_size=3, dilation=2),
            _ConvBlock(hidden_dim,    hidden_dim, kernel_size=3, dilation=4),
            _ConvBlock(hidden_dim,    hidden_dim, kernel_size=3, dilation=8),
        )

        # Dual pooling — captures both peak and average signal
        self.max_pool = nn.AdaptiveMaxPool1d(1)
        self.avg_pool = nn.AdaptiveAvgPool1d(1)

        # MLP projection: hidden_dim*2 → latent_dim (concat of max+avg)
        pool_dim = hidden_dim * 2
        self.fc1    = nn.Linear(pool_dim, hidden_dim)
        self.norm1  = nn.LayerNorm(hidden_dim)
        self.fc_mu     = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

        self._init_weights()

    def _init_weights(self):
        """Xavier init for all linear/conv layers; zero-init logvar bias."""
        for m in self.modules():
            if isinstance(m, (nn.Linear, nn.Conv1d)):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
        # Start logvar near zero → unit Gaussian prior at training start
        nn.init.zeros_(self.fc_logvar.weight)
        nn.init.zeros_(self.fc_logvar.bias)

    def forward(self, x):
        """
        x : LongTensor (B, L)
        Returns mu, logvar — each (B, latent_dim)
        """
        emb = self.embedding(x).permute(0, 2, 1)   # (B, embed, L)
        h   = self.blocks(emb)                       # (B, hidden, L)

        # Dual pooling → concat → (B, hidden*2)
        h_max = self.max_pool(h).squeeze(2)
        h_avg = self.avg_pool(h).squeeze(2)
        h_cat = torch.cat([h_max, h_avg], dim=1)

        # MLP
        h1 = F.gelu(self.norm1(self.fc1(h_cat)))

        mu     = self.fc_mu(h1)
        logvar = self.fc_logvar(h1).clamp(-6.0, 2.0)  # numerical safety

        return mu, logvar
