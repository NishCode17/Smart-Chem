# Molecule Generation Upgrade Plan

## Goal
Improve the quality of generated molecules by implementing a tiered filtering system (Hard vs Soft filters) and a weighted scoring mechanism. Ensure the targeted generator never returns empty results.

## User Review Required
> [!IMPORTANT]
> **Hard Filters**: Molecules with < 5 or > 50 atoms, or pure carbon compositions, will be **rejected outright**.
> **Soft Filters**: Long chains (>7) and extreme LogP will be penalized. High QED (>0.5) and aromatic rings (>=1) will be rewarded.

## Proposed Changes

### Backend (`chem_utils.py`)
#### [MODIFY] [chem_utils.py](file:///d:/Final%20Year%20Project/Smart%20Chem/backend/chem_utils.py)
- **New Helper Function**: `get_longest_carbon_chain_length(mol)` to calculate the length of the longest linear aliphatic carbon chain.
- **New Function**: `score_molecule(mol, stats)` to calculate the fitness score based on:
    - Base score: 1.0
    - Penalty: Long chain (> 7) -> -0.15
    - Penalty: LogP outside [-1, 5] -> -0.2
    - Reward: QED > 0.5 -> +0.2
    - Reward: Aromatic Rings >= 1 -> +0.25
- **Update**: `calculate_properties` to include `score`, `longest_chain`, and `num_aromatic_rings`.

### Backend (`main.py`)
#### [MODIFY] [main.py](file:///d:/Final%20Year%20Project/Smart%20Chem/backend/main.py)
- **Update `is_valid_candidate`**:
    - Enforce Hard Filters strictly:
        - 5 <= Heavy Atoms <= 50
        - Reject Pure Carbon
- **Update `generate_targeted` and `optimize_lead`**:
    - Increase internal batch size to ensure enough candidates pass hard filters.
    - Calculate scores for all valid candidates.
    - Sort results by `score` (descending) instead of just QED distance.
    - **Fallback Logic**: If hard filters reject too many, automatically retry with a larger batch or fallback to random generation if strictly necessary (though increased batch size should suffice).

## Verification Plan
### Automated Tests
- Run the targeted generator endpoint.
- Verify that resulting molecules:
    - Have mixed atom types (no pure carbon).
    - Favor structures with rings and reasonable logP.
    - Are never empty 
- Check logs to see scoring happening.
