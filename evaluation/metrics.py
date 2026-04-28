import random
import torch
from rdkit import Chem
from rdkit.Chem import QED, AllChem
from rdkit import DataStructs


def validity(smiles_list):
    valid = [s for s in smiles_list if Chem.MolFromSmiles(s) is not None]
    return len(valid) / len(smiles_list) if smiles_list else 0.0


def uniqueness(smiles_list):
    if not smiles_list:
        return 0.0
    unique_smiles = set(smiles_list)
    return len(unique_smiles) / len(smiles_list)


def novelty(generated_smiles, training_smiles_set):
    if not generated_smiles:
        return 0.0
    novel = [s for s in generated_smiles if s not in training_smiles_set]
    return len(novel) / len(generated_smiles)


def internal_diversity(smiles_list):
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    mols = [m for m in mols if m is not None]
    if len(mols) < 2:
        return 0.0

    # Sample subset for performance safety
    if len(mols) > 100:
        mols = random.sample(mols, 100)

    fps = [
        AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048)
        for m in mols
    ]

    similarities = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            similarities.append(sim)

    mean_similarity = sum(similarities) / len(similarities) if similarities else 0.0
    return 1.0 - mean_similarity


def property_scores(smiles_list, predictor, device):
    qed_scores = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            try:
                qed_scores.append(QED.qed(mol))
            except Exception:
                pass

    mean_qed = float(sum(qed_scores) / len(qed_scores)) if qed_scores else 0.0
    return {"mean_qed": mean_qed}


def success_rate(smiles_list, threshold=0.7):
    successes = 0
    total = 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            total += 1
            try:
                if QED.qed(mol) > threshold:
                    successes += 1
            except Exception:
                pass
    return successes / total if total > 0 else 0.0


def compute_all_metrics(smiles_list, training_smiles_set, predictor, device, qed_threshold=0.7):
    val = validity(smiles_list)
    valid_smiles = [s for s in smiles_list if Chem.MolFromSmiles(s) is not None]

    uniq = uniqueness(valid_smiles)
    nov = novelty(valid_smiles, training_smiles_set)
    props = property_scores(valid_smiles, predictor, device)
    sr = success_rate(valid_smiles, qed_threshold)
    div = internal_diversity(valid_smiles)

    return {
        "validity": val,
        "uniqueness": uniq,
        "novelty": nov,
        "mean_qed": props["mean_qed"],
        "success_rate": sr,
        "diversity": div,
    }
