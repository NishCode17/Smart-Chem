"""
data/preprocess_zinc.py

Builds aligned, split-ready datasets from raw SMILES for:
  - Sequence   (tokenised SELFIES tensor)
  - Graph      (PyG Data with .seq attached)
  - Hybrid     (seq_tensor, graph_data) pairs

Reuses the existing tokeniser from misc/preprocess.py and the existing
graph builder from data/graph_dataset.py.  No model code is touched.

Usage (standalone):
    python data/preprocess_zinc.py \\
        --input  data/raw/train_molecules.csv \\
        --output data/processed \\
        --target 100000

Usage (as library):
    from data.preprocess_zinc import build_full_dataset
    splits = build_full_dataset("data/raw/train_molecules.csv", target_size=100000)
    # splits["train_seq"], splits["train_graph"], ...
"""

import os
import sys
import json
import argparse
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data, Batch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from rdkit import Chem
import selfies as sf

from data.zinc_loader import load_smiles
from data.graph_dataset import smiles_to_graph

# Reuse tokeniser helpers from the original preprocessor
sys.path.insert(0, os.path.join(ROOT, "misc"))
from preprocess import build_vocab, text_to_tensor


# ---------------------------------------------------------------------------
# SECTION 4 — CLEANING & FILTERING
# ---------------------------------------------------------------------------

def clean_smiles(smiles_list: list, target_size: int = 100_000,
                 max_atoms: int = 60) -> list:
    """
    Deduplicate → validate with RDKit → filter by atom count → stop at target.

    Args:
        smiles_list : raw SMILES (may contain duplicates / invalid entries)
        target_size : stop collecting after this many valid molecules
        max_atoms   : discard molecules with more atoms than this

    Returns:
        Cleaned list of canonical SMILES (length ≤ target_size)
    """
    seen   = set()
    clean  = []

    print(f"[clean_smiles] Scanning {len(smiles_list):,} raw entries "
          f"(target={target_size:,}, max_atoms={max_atoms}) …")

    for smi in smiles_list:
        if len(clean) >= target_size:
            break

        smi = smi.strip()
        if not smi or smi in seen:
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        if mol.GetNumAtoms() > max_atoms:
            continue

        # Canonicalise so duplicates written differently are caught
        canon = Chem.MolToSmiles(mol, canonical=True)
        if canon in seen:
            continue

        seen.add(smi)
        seen.add(canon)
        clean.append(canon)

    print(f"[clean_smiles] Valid molecules collected : {len(clean):,}")
    return clean


# ---------------------------------------------------------------------------
# SECTION 5 — CONVERT TO SELFIES
# ---------------------------------------------------------------------------

def smiles_to_selfies_list(smiles_list: list) -> tuple:
    """
    Convert SMILES → SELFIES, keeping track of which indices succeeded.

    Returns:
        (selfies_list, valid_indices)
        selfies_list  : list of SELFIES strings (same length as valid_indices)
        valid_indices : indices into smiles_list for which conversion succeeded
    """
    selfies_out   = []
    valid_indices = []

    for i, smi in enumerate(smiles_list):
        try:
            sel = sf.encoder(smi)
            if sel:
                selfies_out.append(sel)
                valid_indices.append(i)
        except Exception:
            continue

    print(f"[smiles_to_selfies] {len(smiles_list):,} in  →  "
          f"{len(selfies_out):,} out  "
          f"({len(smiles_list) - len(selfies_out):,} failures skipped)")
    return selfies_out, valid_indices


# ---------------------------------------------------------------------------
# SECTION 6 — TOKENISATION
# ---------------------------------------------------------------------------

def tokenise(selfies_list: list, smiles_list: list,
             max_len: int = 100) -> tuple:
    """
    Build vocab and convert SELFIES to integer tensors.

    Returns:
        seq_tensor : LongTensor (N, max_len)
        vocab      : token → index dict
    """
    print("[tokenise] Building SELFIES vocabulary …")
    vocab = build_vocab(selfies_list, is_selfies=True)

    print(f"[tokenise] Vocab size = {len(vocab)}  |  max_len = {max_len}")
    print("[tokenise] Encoding sequences …")

    tensors = []
    kept    = []
    for i, sel in enumerate(selfies_list):
        row = text_to_tensor(sel, vocab, max_len, is_selfies=True)
        if row is not None:
            tensors.append(row)
            kept.append(i)

    seq_tensor = torch.tensor(tensors, dtype=torch.long)
    print(f"[tokenise] Tensor shape : {tuple(seq_tensor.shape)}")
    return seq_tensor, vocab, kept


# ---------------------------------------------------------------------------
# SECTION 7 — GRAPH DATASET
# ---------------------------------------------------------------------------

class AlignedGraphDataset:
    """
    Holds a list of PyG Data objects with .seq already attached.
    Indexable — supports slicing for train/test split.
    """
    def __init__(self, data_list: list):
        self.data_list = data_list

    def __len__(self):
        return len(self.data_list)

    def __getitem__(self, idx):
        return self.data_list[idx]


def build_graph_dataset(smiles_list: list,
                        seq_tensor: torch.Tensor) -> AlignedGraphDataset:
    """
    Build graph objects and attach the corresponding sequence tensor row.

    Args:
        smiles_list : canonical SMILES (length N)
        seq_tensor  : LongTensor (N, max_len) — aligned with smiles_list

    Returns:
        AlignedGraphDataset of length N  (invalid graphs are DROP NOT skipped
        — caller must pass already-validated SMILES)
    """
    assert len(smiles_list) == len(seq_tensor), (
        f"SMILES ({len(smiles_list)}) and tensor ({len(seq_tensor)}) "
        f"must have equal length."
    )

    data_list = []
    skipped   = 0

    for i, smi in enumerate(smiles_list):
        graph = smiles_to_graph(smi)
        if graph is None:
            skipped += 1
            continue
        # Attach sequence so GNN DataLoader can give it to the decoder.
        # Store as (1, max_len) so PyG batches B graphs → (B, max_len).
        # 1D tensors (max_len,) would be concatenated to (B*max_len,) — wrong.
        graph.seq = seq_tensor[i].unsqueeze(0)   # shape (1, max_len)
        graph.smiles_idx = i                     # traceability
        data_list.append(graph)

    if skipped:
        print(f"[build_graph_dataset] Dropped {skipped} molecules "
              f"(smiles_to_graph returned None)")

    print(f"[build_graph_dataset] Graph dataset size : {len(data_list):,}")
    return AlignedGraphDataset(data_list)


# ---------------------------------------------------------------------------
# SECTION 8 — HYBRID DATASET
# ---------------------------------------------------------------------------

class HybridDataset(Dataset):
    """
    Returns (seq_tensor_row, graph_data) pairs.
    Both items share the same molecule index — strict alignment guaranteed.
    """
    def __init__(self, seq_tensor: torch.Tensor,
                 graph_dataset: AlignedGraphDataset):
        assert len(seq_tensor) == len(graph_dataset), (
            f"Seq ({len(seq_tensor)}) and graph ({len(graph_dataset)}) "
            f"must have equal length."
        )
        self.seq   = seq_tensor
        self.graph = graph_dataset

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, idx):
        return self.seq[idx], self.graph[idx]


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def build_full_dataset(input_path:   str,
                       output_dir:   str  = None,
                       target_size:  int  = 100_000,
                       max_atoms:    int  = 60,
                       max_len:      int  = 100,
                       train_ratio:  float = 0.8) -> dict:
    """
    End-to-end pipeline: raw file → aligned, split datasets.

    On the FIRST run: reads the CSV, cleans, converts to SELFIES, tokenises,
    builds graphs, and saves artefacts to output_dir.

    On SUBSEQUENT RUNS: detects the saved artefacts in output_dir and loads
    them directly — skipping the entire CSV → clean → tokenise pipeline.
    Only the in-memory graph/hybrid datasets are rebuilt from the cached SMILES.

    Args:
        input_path  : path to CSV or TXT file with SMILES
        output_dir  : directory for cached artefacts (train_selfies.pt, etc.)
        target_size : target number of clean molecules (only used on first run)
        max_atoms   : heavy-atom filter limit (only used on first run)
        max_len     : tokenised sequence length (must match VAE max_len)
        train_ratio : train/test split fraction

    Returns:
        dict with keys:
            train_seq, test_seq          : LongTensor subsets
            train_graph, test_graph      : AlignedGraphDataset subsets
            train_hybrid, test_hybrid    : HybridDataset subsets
            train_smiles, test_smiles    : list[str]
            vocab                        : token → index
    """

    # ------------------------------------------------------------------
    # CACHE CHECK — load artefacts if they already exist
    # ------------------------------------------------------------------
    if output_dir:
        pt_path     = os.path.join(output_dir, "train_selfies.pt")
        vocab_path  = os.path.join(output_dir, "selfies_vocab.json")
        smiles_path = os.path.join(output_dir, "train_smiles.txt")

        if all(os.path.isfile(p) for p in [pt_path, vocab_path, smiles_path]):
            print(f"\n{'='*60}")
            print("  ✅ Cached artefacts found — skipping preprocessing")
            print(f"     {pt_path}")
            print(f"     {vocab_path}")
            print(f"     {smiles_path}")
            print(f"{'='*60}")

            seq_tensor  = torch.load(pt_path, weights_only=True)
            with open(vocab_path) as fh:
                vocab = json.load(fh)
            with open(smiles_path) as fh:
                final_smiles = [l.strip() for l in fh if l.strip()]

            N = len(seq_tensor)
            assert len(final_smiles) == N, (
                f"Cached SMILES ({len(final_smiles)}) and tensor ({N}) mismatch. "
                f"Delete data/processed/ artefacts and re-run to rebuild."
            )
            print(f"  Loaded  : {N:,} molecules  |  vocab {len(vocab)} tokens")

            # Skip straight to graph building (always in-memory)
            graph_dataset = build_graph_dataset(final_smiles, seq_tensor)

            # Re-align if any SMILES failed smiles_to_graph
            valid_pos = set(g.smiles_idx for g in graph_dataset.data_list)
            if len(valid_pos) < N:
                final_smiles = [final_smiles[i] for i in range(N) if i in valid_pos]
                seq_tensor   = seq_tensor[[i for i in range(N) if i in valid_pos]]
                N = len(final_smiles)
                graph_dataset = build_graph_dataset(final_smiles, seq_tensor)

            n_train = int(N * train_ratio)
            n_test  = N - n_train

            train_seq    = seq_tensor[:n_train]
            test_seq     = seq_tensor[n_train:]
            train_smiles = final_smiles[:n_train]
            test_smiles  = final_smiles[n_train:]
            train_graph  = AlignedGraphDataset(graph_dataset.data_list[:n_train])
            test_graph   = AlignedGraphDataset(graph_dataset.data_list[n_train:])
            train_hybrid = HybridDataset(train_seq, train_graph)
            test_hybrid  = HybridDataset(test_seq,  test_graph)

            print(f"  Train   : {n_train:,}  |  Test (held-out) : {n_test:,}\n")
            return {
                "train_seq": train_seq, "test_seq": test_seq,
                "train_graph": train_graph, "test_graph": test_graph,
                "train_hybrid": train_hybrid, "test_hybrid": test_hybrid,
                "train_smiles": train_smiles, "test_smiles": test_smiles,
                "vocab": vocab,
            }

    # ------------------------------------------------------------------
    # 1. Load raw SMILES (only reached on first run)
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Loading raw SMILES from  {input_path}")
    print(f"{'='*60}")
    raw = load_smiles(input_path)
    print(f"  Raw entries loaded       : {len(raw):,}")


    # ------------------------------------------------------------------
    # 2. Clean & filter
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  Cleaning & filtering")
    print(f"{'='*60}")
    clean = clean_smiles(raw, target_size=target_size, max_atoms=max_atoms)

    if len(clean) < target_size:
        print(f"  ⚠️  Only {len(clean):,} valid molecules found "
              f"(target was {target_size:,}). Proceeding with available data.")
    else:
        print(f"  ✅ Target met: {len(clean):,} molecules")

    # ------------------------------------------------------------------
    # 3. SMILES → SELFIES
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  Converting to SELFIES")
    print(f"{'='*60}")
    selfies_list, valid_idx = smiles_to_selfies_list(clean)

    # Keep only SMILES that survived SELFIES conversion (strict alignment)
    aligned_smiles = [clean[i] for i in valid_idx]

    # ------------------------------------------------------------------
    # 4. Tokenise
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  Tokenising")
    print(f"{'='*60}")
    seq_tensor, vocab, kept_idx = tokenise(selfies_list, aligned_smiles,
                                           max_len=max_len)

    # Final aligned SMILES (after tokeniser dropped any remaining failures)
    final_smiles = [aligned_smiles[i] for i in kept_idx]
    N = len(final_smiles)

    print(f"\n  Final aligned dataset size : {N:,}")
    assert N == len(seq_tensor), "Alignment broken between tensor and SMILES"

    # ------------------------------------------------------------------
    # 5. Graph dataset (uses same positions as seq_tensor)
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  Building graph dataset")
    print(f"{'='*60}")
    graph_dataset = build_graph_dataset(final_smiles, seq_tensor)

    # Some SMILES may still fail smiles_to_graph → re-align everything
    # by filtering through the indices that produced valid graphs
    valid_graph_positions = [
        g.smiles_idx for g in graph_dataset.data_list
    ]

    if len(valid_graph_positions) < N:
        pos_set = set(valid_graph_positions)
        final_smiles = [final_smiles[i] for i in range(N) if i in pos_set]
        seq_tensor   = seq_tensor[[i for i in range(N) if i in pos_set]]
        N = len(final_smiles)
        print(f"  After graph filter : {N:,} molecules remain")

        # Rebuild graphs with fresh .smiles_idx after re-indexing
        graph_dataset = build_graph_dataset(final_smiles, seq_tensor)

    # ------------------------------------------------------------------
    # 6. Train / test split  (deterministic — NO shuffle)
    # ------------------------------------------------------------------
    n_train = int(N * train_ratio)
    n_test  = N - n_train

    print(f"\n{'='*60}")
    print(f"  Train/Test split  ({int(train_ratio*100)}/{int((1-train_ratio)*100)})")
    print(f"{'='*60}")
    print(f"  Train : {n_train:,}")
    print(f"  Test  : {n_test:,}")

    train_seq    = seq_tensor[:n_train]
    test_seq     = seq_tensor[n_train:]

    train_smiles = final_smiles[:n_train]
    test_smiles  = final_smiles[n_train:]

    train_graph  = AlignedGraphDataset(graph_dataset.data_list[:n_train])
    test_graph   = AlignedGraphDataset(graph_dataset.data_list[n_train:])

    train_hybrid = HybridDataset(train_seq, train_graph)
    test_hybrid  = HybridDataset(test_seq,  test_graph)

    # ------------------------------------------------------------------
    # 7. Optional: save artefacts
    # ------------------------------------------------------------------
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        # Sequence tensor (full — same convention as original pipeline)
        torch.save(seq_tensor,
                   os.path.join(output_dir, "train_selfies.pt"))

        # Vocab
        with open(os.path.join(output_dir, "selfies_vocab.json"), "w") as fh:
            json.dump(vocab, fh)

        # SMILES text file (all, aligned with tensor)
        with open(os.path.join(output_dir, "train_smiles.txt"), "w") as fh:
            for s in final_smiles:
                fh.write(s + "\n")

        print(f"\n  Artefacts saved → {output_dir}/")
        print(f"    train_selfies.pt     ({N:,} × {max_len})")
        print(f"    selfies_vocab.json   ({len(vocab)} tokens)")
        print(f"    train_smiles.txt     ({N:,} SMILES)")

    # ------------------------------------------------------------------
    # 8. Sanity check
    # ------------------------------------------------------------------
    assert len(train_seq)    == len(train_smiles) == len(train_graph) == len(train_hybrid)
    assert len(test_seq)     == len(test_smiles)  == len(test_graph)  == len(test_hybrid)

    print(f"\n{'='*60}")
    print("  ✅ Dataset pipeline complete")
    print(f"  Total molecules   : {N:,}")
    print(f"  Train             : {n_train:,}")
    print(f"  Test (held-out)   : {n_test:,}")
    print(f"{'='*60}\n")

    return {
        "train_seq":    train_seq,
        "test_seq":     test_seq,
        "train_graph":  train_graph,
        "test_graph":   test_graph,
        "train_hybrid": train_hybrid,
        "test_hybrid":  test_hybrid,
        "train_smiles": train_smiles,
        "test_smiles":  test_smiles,
        "vocab":        vocab,
    }


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess ZINC SMILES → aligned seq/graph/hybrid datasets"
    )
    parser.add_argument("--input",  required=True,
                        help="Path to raw CSV or TXT file with SMILES")
    parser.add_argument("--output", default="data/processed",
                        help="Directory to save processed artefacts")
    parser.add_argument("--target", type=int, default=100_000,
                        help="Target number of clean molecules (default 100000)")
    parser.add_argument("--max_atoms", type=int, default=60,
                        help="Discard molecules with more heavy atoms (default 60)")
    parser.add_argument("--max_len",   type=int, default=100,
                        help="Tokenised sequence length (default 100)")
    args = parser.parse_args()

    build_full_dataset(
        input_path  = args.input,
        output_dir  = args.output,
        target_size = args.target,
        max_atoms   = args.max_atoms,
        max_len     = args.max_len,
    )
