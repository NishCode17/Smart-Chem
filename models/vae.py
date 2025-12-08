import torch
import torch.nn as nn
import torch.nn.functional as F

class VAE(nn.Module):
    def __init__(self, vocab_size, embedding_dim=128, hidden_dim=256, latent_dim=128, max_len=100):
        super(VAE, self).__init__()
        self.max_len = max_len
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim

        # Encoder
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.conv1 = nn.Conv1d(embedding_dim, 9, kernel_size=9)
        self.conv2 = nn.Conv1d(9, 9, kernel_size=9)
        self.conv3 = nn.Conv1d(9, 10, kernel_size=11)
        self.flatten = nn.Flatten()
        self.adaptive_pool = nn.AdaptiveMaxPool1d(1) 
        self.fc_mu = nn.Linear(10, latent_dim)
        self.fc_logvar = nn.Linear(10, latent_dim)

        # Decoder
        self.decoder_input = nn.Linear(latent_dim, hidden_dim)
        self.gru = nn.GRU(embedding_dim, hidden_dim, num_layers=3, batch_first=True)
        self.fc_out = nn.Linear(hidden_dim, vocab_size)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        embedded = self.embedding(x).permute(0, 2, 1)
        c1 = F.relu(self.conv1(embedded))
        c2 = F.relu(self.conv2(c1))
        c3 = F.relu(self.conv3(c2))                 
        pooled = self.adaptive_pool(c3).squeeze(2)  
        mu = self.fc_mu(pooled)
        logvar = self.fc_logvar(pooled)
        z = self.reparameterize(mu, logvar)
        hidden = self.decoder_input(z).unsqueeze(0).repeat(3, 1, 1)
        output, _ = self.gru(embedded.permute(0, 2, 1), hidden)
        return self.fc_out(output), mu, logvar

    def sample(self, num_samples, device, temperature=1.5):
        """Random Generation with Diversity"""
        z = torch.randn(num_samples, self.latent_dim).to(device)
        return self._decode_loop(z, device, temperature)

    def decode(self, z, device, temperature=0.8):
        """Targeted Generation with Stability"""
        return self._decode_loop(z, device, temperature)

    def _decode_loop(self, z, device, temperature):
        # Shared logic to prevent code duplication
        hidden = self.decoder_input(z).unsqueeze(0).repeat(3, 1, 1)
        input_seq = torch.ones(z.shape[0], 1).long().to(device)
        generated_indices = []
        
        for _ in range(self.max_len):
            embedded = self.embedding(input_seq)
            output, hidden = self.gru(embedded, hidden)
            logits = self.fc_out(output)
            
            # FORCE SAMPLING (Fixes Carbon Chains)
            probs = F.softmax(logits.squeeze(1) / temperature, dim=1)
            next_token = torch.multinomial(probs, 1)
            
            generated_indices.append(next_token)
            input_seq = next_token
            
        return torch.cat(generated_indices, dim=1)