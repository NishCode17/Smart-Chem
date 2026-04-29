import torch
import torch.nn as nn
import torch.nn.functional as F
from models.cnn_encoder import CNNEncoder


class VAE(nn.Module):
    def __init__(self, vocab_size, embedding_dim=128, hidden_dim=256,
                 latent_dim=128, max_len=100):
        super(VAE, self).__init__()
        self.max_len    = max_len
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim
        self.num_layers = 1

        # Encoder
        self.encoder = CNNEncoder(vocab_size, embedding_dim, latent_dim,
                                  hidden_dim=hidden_dim)

        # Decoder
        self.embedding_dim = embedding_dim
        self.embedding     = nn.Embedding(vocab_size, embedding_dim, padding_idx=0)
        self.decoder_input = nn.Linear(latent_dim, hidden_dim * self.num_layers)
        self.gru    = nn.GRU(embedding_dim + latent_dim, hidden_dim,
                             num_layers=self.num_layers,
                             batch_first=True, dropout=0.0)
        self.fc_out = nn.Linear(hidden_dim, vocab_size)

        self._init_decoder()

    def _init_decoder(self):
        for name, p in self.gru.named_parameters():
            if "weight" in name:
                nn.init.orthogonal_(p)
            elif "bias" in name:
                nn.init.zeros_(p)
        nn.init.xavier_uniform_(self.fc_out.weight)
        nn.init.zeros_(self.fc_out.bias)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std) * 1.1
        return mu + eps * std

    def _hidden_from_z(self, z):
        """Project z → (num_layers, B, hidden_dim) GRU initial state."""
        h = self.decoder_input(z)                              # (B, hidden * num_layers)
        h = h.view(z.size(0), self.num_layers, self.hidden_dim)  # (B, L, hidden)
        return h.permute(1, 0, 2).contiguous()                 # (L, B, hidden)

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z          = self.reparameterize(mu, logvar)
        hidden     = self._hidden_from_z(z)
        emb        = self.embedding(x)                              # (B, L, emb_dim)
        z_expanded = z.unsqueeze(1).repeat(1, emb.size(1), 1)      # (B, L, latent_dim)
        dec_in     = torch.cat([emb, z_expanded], dim=-1)          # (B, L, emb+latent)
        out, _     = self.gru(dec_in, hidden)
        return self.fc_out(out), mu, logvar

    def sample(self, num_samples, device, temperature=1.0):
        z = torch.randn(num_samples, self.latent_dim).to(device)
        return self._decode_loop(z, device, temperature)

    def decode(self, z, device, temperature=0.8):
        return self._decode_loop(z, device, temperature)

    def _decode_loop(self, z, device, temperature):
        hidden    = self._hidden_from_z(z)
        input_seq = torch.ones(z.shape[0], 1, dtype=torch.long, device=device)
        generated = []
        for _ in range(self.max_len):
            emb        = self.embedding(input_seq)             # (B, 1, emb_dim)
            z_step     = z.unsqueeze(1)                        # (B, 1, latent_dim)
            dec_in     = torch.cat([emb, z_step], dim=-1)     # (B, 1, emb+latent)
            out, hidden = self.gru(dec_in, hidden)
            logits      = self.fc_out(out)
            probs       = F.softmax(logits.squeeze(1) / temperature, dim=-1)
            next_token  = torch.multinomial(probs, 1)
            generated.append(next_token)
            input_seq   = next_token
        return torch.cat(generated, dim=1)