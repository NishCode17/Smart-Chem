import torch
import torch.nn as nn
import torch.nn.functional as F

class VAE(nn.Module):
    def __init__(self, vocab_size, embedding_dim=128, hidden_dim=256, latent_dim=128, max_len=100):
        super(VAE, self).__init__()
        self.max_len = max_len
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim

        # --- ENCODER (CNN) ---
        # "We used 3 convolutional layers..." (Gomez-Bombarelli et al.)
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        
        self.conv1 = nn.Conv1d(embedding_dim, 9, kernel_size=9)
        self.conv2 = nn.Conv1d(9, 9, kernel_size=9)
        self.conv3 = nn.Conv1d(9, 10, kernel_size=11)
        
        self.flatten = nn.Flatten()
        
        # Calculate flat size after convs (depends on max_len)
        # For max_len=100 and kernel sizes above, we estimate linear input:
        # This is a simplification; in production we calculate exact shapes.
        # Let's use an adaptive pool to force a fixed size for simplicity in MVP.
        self.adaptive_pool = nn.AdaptiveMaxPool1d(1) 
        
        # Latent Space Projections
        self.fc_mu = nn.Linear(10, latent_dim)      # Mean
        self.fc_logvar = nn.Linear(10, latent_dim)  # Log Variance

        # --- DECODER (RNN/GRU) ---
        self.decoder_input = nn.Linear(latent_dim, hidden_dim)
        self.gru = nn.GRU(embedding_dim, hidden_dim, num_layers=3, batch_first=True)
        self.fc_out = nn.Linear(hidden_dim, vocab_size)

    def reparameterize(self, mu, logvar):
        """
        The 'Reparameterization Trick':
        z = mu + sigma * epsilon
        Allows backprop through random sampling.
        """
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        # x shape: [batch, seq_len]
        
        # 1. Encode
        embedded = self.embedding(x)                # [batch, seq, embed]
        embedded = embedded.permute(0, 2, 1)        # [batch, embed, seq] (Conv1d expects channels first)
        
        c1 = F.relu(self.conv1(embedded))
        c2 = F.relu(self.conv2(c1))
        c3 = F.relu(self.conv3(c2))                 # [batch, 10, reduced_seq]
        
        pooled = self.adaptive_pool(c3).squeeze(2)  # [batch, 10]
        
        mu = self.fc_mu(pooled)
        logvar = self.fc_logvar(pooled)
        
        # 2. Latent Sample
        z = self.reparameterize(mu, logvar)
        
        # 3. Decode
        # We need to expand z to be the initial hidden state or input
        # Here we treat z as the initial hidden state for the GRU
        hidden = self.decoder_input(z).unsqueeze(0).repeat(3, 1, 1) # [3 layers, batch, hidden]
        
        # For training, we use Teacher Forcing (feed actual input)
        # Note: In strict VAE, we usually decode from z. 
        # For this MVP 3-day build, we feed embedded input to GRU and condition hidden on z.
        output, _ = self.gru(embedded.permute(0, 2, 1), hidden)
        
        prediction = self.fc_out(output) # [batch, seq, vocab]
        
        return prediction, mu, logvar

    def sample(self, num_samples, device):
        """Generates new molecules using Temperature Sampling (Creativity)"""
        z = torch.randn(num_samples, self.latent_dim).to(device)
        
        # Initialize hidden state
        hidden = self.decoder_input(z).unsqueeze(0).repeat(3, 1, 1)
        
        # Start token (<sos>)
        input_seq = torch.ones(num_samples, 1).long().to(device) 
        
        generated_indices = []
        
        for _ in range(self.max_len):
            embedded = self.embedding(input_seq)
            output, hidden = self.gru(embedded, hidden)
            logits = self.fc_out(output) # [batch, 1, vocab]
            
            # --- THE FIX: SAMPLING INSTEAD OF ARGMAX ---
            # Apply softmax to turn logits into probabilities
            probs = F.softmax(logits.squeeze(1), dim=1)
            
            # Sample from the distribution (Roll the dice!)
            next_token = torch.multinomial(probs, 1)
            
            generated_indices.append(next_token)
            
            input_seq = next_token # Feed as input for next step
            
        return torch.cat(generated_indices, dim=1)