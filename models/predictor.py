import torch
import torch.nn as nn

class PropertyPredictor(nn.Module):
    def __init__(self, latent_dim=128):
        super(PropertyPredictor, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 3)
        )

    def forward(self, z):
        return self.net(z)