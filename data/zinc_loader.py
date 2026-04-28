"""
data/zinc_loader.py

Loads raw SMILES strings from a local file.

Supported formats:
- CSV  : must have a column named "smiles" (case-insensitive)
- TXT  : one SMILES per line (plain text)

Usage:
    from data.zinc_loader import load_smiles
    smiles = load_smiles("data/raw/train_molecules.csv")
"""

import os
import csv


def load_smiles(filepath: str) -> list:
    """
    Load raw SMILES from a CSV or plain-text file.

    Args:
        filepath : absolute or relative path to the data file.

    Returns:
        List of raw SMILES strings (un-validated, un-cleaned).

    Raises:
        FileNotFoundError : if filepath does not exist.
        ValueError        : if a CSV has no recognisable 'smiles' column.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Data file not found: {filepath}")

    ext = os.path.splitext(filepath)[1].lower()

    # ------------------------------------------------------------------ CSV
    if ext == ".csv":
        smiles_list = []
        with open(filepath, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)

            # Locate the smiles column (case-insensitive)
            header = [c.lower() for c in reader.fieldnames]
            if "smiles" not in header:
                raise ValueError(
                    f"CSV '{filepath}' has no 'smiles' column. "
                    f"Found: {reader.fieldnames}"
                )
            col_idx = header.index("smiles")
            col_name = reader.fieldnames[col_idx]

            for row in reader:
                smi = row[col_name]
                if smi and smi.strip():
                    smiles_list.append(smi.strip())

        return smiles_list

    # ------------------------------------------------- Plain text (.txt / other)
    smiles_list = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            smi = line.strip()
            if smi:
                smiles_list.append(smi)

    return smiles_list
