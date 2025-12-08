from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Draw, rdMolDescriptors, AllChem
import selfies as sf
import base64
from io import BytesIO

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
            "admet": admet, # <--- NEW FIELD
            "status": status,
            "image": f"data:image/png;base64,{img_str}"
        }
    except Exception as e:
        print(f"Prop Error: {e}")
        return {"valid": False}