import pandas as pd
import selfies as sf
import torch
import json
import os
from tqdm import tqdm

# --- Configuration ---
# Update these paths to match where you put your downloaded files
TRAIN_PATH = "data/raw/train_molecules.csv" 
OUTPUT_DIR = "data/processed"

# Exact column names from your CSV
COL_SMILES = "smiles"   
COL_SELFIES = "SELFIES" 

MAX_LEN = 100           
MVP_LIMIT = 150000       # <--- THIS IS THE LINE YOU ASKED FOR

os.makedirs(OUTPUT_DIR, exist_ok=True)

def build_vocab(sequences, is_selfies):
    """
    Scans the entire list of molecules to create a dictionary of unique tokens.
    """
    vocab_set = set()
    for seq in tqdm(sequences, desc=f"Building {'SELFIES' if is_selfies else 'SMILES'} Vocab"):
        if is_selfies:
            try:
                # SELFIES splits by bracket tokens: [C][Branch1]...
                tokens = list(sf.split_selfies(seq))
                vocab_set.update(tokens)
            except:
                continue
        else:
            # SMILES splits by character
            vocab_set.update(list(seq))
            
    vocab = sorted(list(vocab_set))
    
    # Define special tokens
    # <pad>: 0 (Padding)
    # <sos>: 1 (Start of Sequence)
    # <eos>: 2 (End of Sequence)
    token_to_idx = { "<pad>": 0, "<sos>": 1, "<eos>": 2 }
    
    # Map real tokens starting from index 3
    for i, token in enumerate(vocab):
        token_to_idx[token] = i + 3
        
    return token_to_idx

def text_to_tensor(sequence, vocab, max_len, is_selfies):
    """
    Converts a single molecule string into a list of integers (Tensor).
    """
    try:
        if is_selfies:
            tokens = list(sf.split_selfies(sequence))
        else:
            tokens = list(sequence)
            
        # Add Start (<sos>) and End (<eos>) tokens
        indices = [vocab["<sos>"]] + [vocab[t] for t in tokens if t in vocab] + [vocab["<eos>"]]
        
        # Pad or Truncate to fixed length
        if len(indices) < max_len:
            indices += [vocab["<pad>"]] * (max_len - len(indices))
        else:
            indices = indices[:max_len]
            
        return indices
    except:
        return None

def process():
    print(f"Loading Training Data from {TRAIN_PATH}...")
    
    # --- THIS IS HOW YOU TAKE ONLY 20k-30k ENTRIES ---
    # We use 'nrows' to read only the top N rows from the CSV.
    # This keeps everything fast for your 3-day build.
    df = pd.read_csv(TRAIN_PATH, nrows=MVP_LIMIT)
    
    print(f"✅ Successfully loaded top {len(df)} rows.")

    # 1. Filter by length (Speed optimization)
    # We discard huge molecules that would slow down training
    mask = (df[COL_SMILES].str.len() < MAX_LEN) & (df[COL_SELFIES].str.len() < MAX_LEN)
    df = df[mask]
    print(f"Dataset size after length filtering: {len(df)}")

    # 2. Build Vocabularies (The Dictionaries)
    # We learn the "language" of chemistry from the training data
    smiles_vocab = build_vocab(df[COL_SMILES].tolist(), is_selfies=False)
    selfies_vocab = build_vocab(df[COL_SELFIES].tolist(), is_selfies=True)
    
    print(f"SMILES Vocab Size: {len(smiles_vocab)}")
    print(f"SELFIES Vocab Size: {len(selfies_vocab)}")

    # 3. Convert Text to Tensors (Numbers)
    print("Converting SMILES to Tensors...")
    smiles_tensors = [text_to_tensor(s, smiles_vocab, MAX_LEN, False) for s in df[COL_SMILES]]
    # Remove failures
    smiles_tensors = [t for t in smiles_tensors if t is not None]
    
    print("Converting SELFIES to Tensors...")
    selfies_tensors = [text_to_tensor(s, selfies_vocab, MAX_LEN, True) for s in df[COL_SELFIES]]
    selfies_tensors = [t for t in selfies_tensors if t is not None]

    # 4. Save Artifacts to Disk
    print("Saving processed files to 'data/processed/'...")
    
    # Save Dictionaries (JSON) - Your Backend/API needs these later!
    with open(f"{OUTPUT_DIR}/smiles_vocab.json", "w") as f: 
        json.dump(smiles_vocab, f)
    with open(f"{OUTPUT_DIR}/selfies_vocab.json", "w") as f: 
        json.dump(selfies_vocab, f)

    # Save Data (PyTorch Files) - Your Training Script needs these!
    torch.save(torch.tensor(smiles_tensors, dtype=torch.long), f"{OUTPUT_DIR}/train_smiles.pt")
    torch.save(torch.tensor(selfies_tensors, dtype=torch.long), f"{OUTPUT_DIR}/train_selfies.pt")
    
    print("✅ Module 1 Complete. Data is ready for the Model.")

if __name__ == "__main__":
    process()