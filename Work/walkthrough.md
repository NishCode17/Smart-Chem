# Molecule Generation Upgrade Walkthrough

## Goal
To improve the quality of generated molecules by implementing a tiered filtering system (Hard vs Soft) and a weighted scoring mechanism, ensuring only valid, high-quality candidates are presented to the user.

## Changes Implemented

### 1. Tiered Filtering Logic (`backend/main.py`)
I implemented a strict "Hard Filter" that outright rejects invalid molecules before they are scored.

```python
# --- STRICT FILTER ---
def is_valid_candidate(mol):
    if mol is None: return False
    
    # 1. Size Rule (STRICT)
    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms < 5 or num_atoms > 50: return False
    
    # 2. Composition Rule (NO PURE CARBON)
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    unique_atoms = set(atoms)
    if len(unique_atoms) == 1 and "C" in unique_atoms:
        return False # Reject Pure Carbon
        
    return True
```

### 2. Weighted Scoring System (`backend/chem_utils.py`)
I added a new scoring function that rewards drug-like features and penalizes undesirable structures.

| Feature | Condition | Impact |
| :--- | :--- | :--- |
| **Base Score** | - | `1.0` |
| **Long Chain** | Length > 7 | `-0.15` |
| **Bad LogP** | < -1 or > 5 | `-0.20` |
| **High QED** | > 0.5 | `+0.20` |
| **Aromatic Rings** | >= 1 | `+0.25` |

```python
def score_molecule(mol, props):
    score = 1.0
    if props['longest_chain'] > 7: score -= 0.15
    if props['logp'] < -1 or props['logp'] > 5: score -= 0.2
    if props['qed'] > 0.5: score += 0.2
    if props['aromatic_rings'] >= 1: score += 0.25
    return round(score, 3)
```

### 3. Robust Generation Flow (`backend/main.py`)
The targeted generator now:
1.  Generates a **massive pool** of candidates (300+).
2.  Applies **Hard Filters** to remove junk.
3.  Calculates **Soft Scores** for the survivors.
4.  **Sorts** by Score (Highest First).
5.  **Fallback Guarantee**: If targeted generation fails (rare), it automatically falls back to random generation to ensure the user always sees results.

## Verification
- **Hard Filters**: Verified that pure carbon chains and tiny/huge molecules are rejected.
- **Scoring**: Valid molecules are now ranked by their drug-likeness.
- **Fallback**: The system will never return an empty list.

## Next Steps
- Monitor the "Score" values in the UI (if exposed) to see how well the model is performing.
- Adjust weights (e.g., penalties for rings > 4) if needed based on chemist feedback.
