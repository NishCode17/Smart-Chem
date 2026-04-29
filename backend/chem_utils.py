from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Draw, rdMolDescriptors, AllChem
from rdkit.Chem import FilterCatalog
from backend.admet import compute_admet
import selfies as sf
import base64
from io import BytesIO
import math

# Filter Catalogs (Singleton)
_filters_loaded = False
_catalog = None

def load_filters():
    global _filters_loaded, _catalog
    if not _filters_loaded:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
        _catalog = FilterCatalog.FilterCatalog(params)
        _filters_loaded = True
    return _catalog

def get_3d_mol_block(smiles):
    """
    Generate 3D MOL block.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        
        mol = Chem.AddHs(mol)
        # Generate 3D conformation
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)
        
        return Chem.MolToMolBlock(mol)
    except Exception as e:
        print(f"3D Error: {e}")
        return None


def get_mol_from_sequence(seq, mode="selfies"):
    mol = None
    try:
        if mode == "selfies":
            smiles = sf.decoder(seq)
            mol = Chem.MolFromSmiles(smiles)
        else:
            mol = Chem.MolFromSmiles(seq)
    except:
        return None
    return mol

def calculate_admet(mol):
    """
    Calculate ADMET properties.
    """
    if mol is None: return {}
    
    # 1. Lipinski Rule of 5 Parameters
    mw = Descriptors.MolWt(mol)           # Target: < 500
    logp = Descriptors.MolLogP(mol)       # Target: < 5
    hbd = Descriptors.NumHDonors(mol)     # Target: < 5
    hba = Descriptors.NumHAcceptors(mol)  # Target: < 10
    
    # 2. Absorption
    tpsa = Descriptors.TPSA(mol)          # Target: < 140 for good absorption
    
    # 3. Size / Complexity
    rotatable = Descriptors.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Check Lipinski Pass/Fail
    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    lipinski_status = "✅ Pass" if violations <= 1 else f"⚠️ {violations} Violations"

    return {
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable": rotatable,
        "rings": rings,
        "lipinski": lipinski_status
    }

def estimate_admet_scores(mol):
    """
    Returns ADMET summary scores.
    """
    if mol is None:
        return {}
    return compute_admet(mol).get("summary", {})

def check_toxicity_alerts(mol):
    """
    Check toxicity alerts using RDKit.
    """
    if mol is None: return {}
    catalog = load_filters()
    
    entries = catalog.GetMatches(mol)
    alerts = []
    
    has_pains = False
    has_brenk = False
    has_nih = False
    
    for entry in entries:
        desc = entry.GetDescription()
        alerts.append(desc)
        if "PAINS" in desc: has_pains = True
        if "brenk" in desc.lower(): has_brenk = True
        if "nih" in desc.lower(): has_nih = True
        
    return {
        "pains": has_pains,
        "brenk": has_brenk,
        "nih": has_nih,
        "alerts_count": len(alerts),
        "details": alerts[:5] # Top 5 alerts
    }

# Scoring
def get_longest_carbon_chain_length(mol):
    """
    Get longest carbon chain length.
    """
    try:
        # Match aliphatic carbons
        max_len = 0
        for length in range(20, 2, -1): # Check down from 20
            query = "[C;R0]" + ("[C;R0]" * (length - 1))
            patt = Chem.MolFromSmarts(query)
            if mol.HasSubstructMatch(patt):
                max_len = length
                break
        return max_len
    except:
        return 0

def score_molecule(mol, props):
    """
    Calculate fitness score.
    """
    score = 1.0
    
    # Extract needed props
    longest_chain = props.get('longest_chain', 0)
    logp = props.get('logp', 0.0)
    qed = props.get('qed', 0.0)
    rings = props.get('aromatic_rings', 0)
    
    # 1. Penalty: Long aliphatic chain
    if longest_chain > 7:
        score -= 0.15
        
    # 2. Penalty: LogP
    if logp < -1 or logp > 5:
        score -= 0.2
        
    # 3. Reward: QED
    if qed > 0.5:
        score += 0.2
        
    # 4. Reward: Aromatic rings
    if rings >= 1:
        score += 0.25
        
    return round(score, 3)

def calculate_properties(mol):
    """
    Calculate molecule properties.
    """
    if mol is None:
        return {"valid": False, "image": None, "status": "Error"}

    try:
        # Basic Props
        qed = round(QED.qed(mol), 3)
        logp = round(Descriptors.MolLogP(mol), 2)
        
        # Additional Structural Features
        longest_chain = get_longest_carbon_chain_length(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        
        # Filter extreme LogP
        if logp < -5 or logp > 8:
             return {"valid": False, "status": "Extreme LogP"}

        admet = calculate_admet(mol)
        
        # Scoring
        temp_props = {
            "longest_chain": longest_chain,
            "logp": logp,
            "qed": qed,
            "aromatic_rings": aromatic_rings
        }
        score = score_molecule(mol, temp_props)
        
        # Determine Status
        status = "🟢 Lead-Like" if score >= 1.2 else ("🟡 Promising" if score >= 0.9 else "🔴 Poor Fit")
        if qed < 0.4: status = "🔴 Poor Candidate"

        # Generate Image
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        admet_full = compute_admet(mol)

        return {
            "valid": True,
            "smiles": Chem.MolToSmiles(mol),
            "qed": qed,
            "logp": logp,
            "score": score,
            "longest_chain": longest_chain,
            "aromatic_rings": aromatic_rings,
            "admet":      admet_full.get("summary", {}),  # 5-score summary
            "admet_full": admet_full,                      # full admet profile
            "admet_props": calculate_admet(mol),           # raw properties
            "tox_alerts": check_toxicity_alerts(mol),
            "status": status,
            "image": f"data:image/png;base64,{img_str}"
        }
    except Exception as e:
        print(f"Prop Error: {e}")
        return {"valid": False}