from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Draw, rdMolDescriptors, AllChem
from rdkit.Chem import FilterCatalog
import selfies as sf
import base64
from io import BytesIO
import math

# --- Initialize Filter Catalogs (Singleton) ---
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
    Generates a 3D MOL block for a given SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        
        mol = Chem.AddHs(mol)
        # Generate 3D Conformation
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
    Calculates detailed ADMET-like properties for the 'Lab View'.
    """
    if mol is None: return {}
    
    # 1. Lipinski Rule of 5 Parameters
    mw = Descriptors.MolWt(mol)           # Target: < 500
    logp = Descriptors.MolLogP(mol)       # Target: < 5
    hbd = Descriptors.NumHDonors(mol)     # Target: < 5
    hba = Descriptors.NumHAcceptors(mol)  # Target: < 10
    
    # 2. Absorption / Permeability
    tpsa = Descriptors.TPSA(mol)          # Target: < 140 for good absorption
    
    # 3. Safety / Complexity
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
    Estimates 0-1 scores for ADMET properties based on functional heuristics.
    These are approximations, not clinical predictions.
    """
    if mol is None: return {}

    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    mw = Descriptors.MolWt(mol)
    
    # 1. Absorption (Human Intestinal Absorption)
    # TPSA < 140 is generally good. Sigmoid decay.
    absorption = 1.0 / (1.0 + math.exp((tpsa - 140) / 15))
    
    # 2. Distribution (BBB Penetration)
    # TPSA < 90 and LogP 2-4 is ideal. 
    # Gaussian-ish peak around LogP ~3
    bbb_score = math.exp(-0.5 * ((logp - 3.0) / 1.5)**2)
    if tpsa > 90: bbb_score *= 0.5
    
    # 3. Metabolism (CYP Stability - crude proxy)
    # High lipophilicity often correlates with clearance.
    metabolism = 1.0 / (1.0 + math.exp((logp - 4.0) / 1.0)) # Stable if hydrophilic
    
    # 4. Excretion (Half-life proxy)
    # Smaller molecules often cleared faster/renally.
    excretion = 1.0 / (1.0 + math.exp((mw - 400) / 100))
    
    # 5. Toxicity (General)
    # High LogP (>5) is toxic.
    toxicity_score = 1.0 / (1.0 + math.exp(-(logp - 5.0))) # Higher is more toxic
    
    return {
        "absorption": round(absorption, 2),
        "distribution": round(bbb_score, 2),
        "metabolism": round(metabolism, 2),
        "excretion": round(excretion, 2),
        "toxicity": round(toxicity_score, 2)
    }

def check_toxicity_alerts(mol):
    """
    Uses RDKit FilterCatalog to check for PAINS, Brenk, and NIH alerts.
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

def calculate_properties(mol):
    """
    Standard properties + Image for cards
    """
    if mol is None:
        return {"valid": False, "image": None, "status": "Error"}

    try:
        # Basic Props
        qed = round(QED.qed(mol), 3)
        logp = round(Descriptors.MolLogP(mol), 2)
        admet = calculate_admet(mol)
        
        # Determine Status
        status = "🟢 Lead-Like" if qed > 0.6 else "🟡 Needs Optimization"
        if qed < 0.4: status = "🔴 Poor Candidate"

        # Generate Image
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        return {
            "valid": True,
            "smiles": Chem.MolToSmiles(mol),
            "qed": qed,
            "logp": logp,
            "logp": logp,
            "admet": estimate_admet_scores(mol),         # Scored (0-1) for Charts
            "admet_props": calculate_admet(mol),         # Raw Properties (MW, TPSA)
            "tox_alerts": check_toxicity_alerts(mol),    # Structural Alerts
            "status": status,
            "image": f"data:image/png;base64,{img_str}"
        }
    except Exception as e:
        print(f"Prop Error: {e}")
        return {"valid": False}