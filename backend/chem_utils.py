from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Draw
import selfies as sf
import base64
from io import BytesIO

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

def calculate_properties(mol):
    # Default Error State
    default_props = {
        "valid": False,
        "smiles": "Generation Failed",
        "qed": 0.0,
        "logp": 0.0,
        "status": "🔴 Error",
        "image": None
    }

    if mol is None: return default_props

    try:
        # 1. properties
        qed = round(QED.qed(mol), 3)
        logp = round(Descriptors.MolLogP(mol), 2)
        
        status = "🟢 High Potential" if qed > 0.6 else "🟡 Average"
        if qed < 0.4: status = "🔴 Low Drug-likeness"

        # 2. IMAGE GENERATION (Crucial Step)
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        return {
            "valid": True,
            "smiles": Chem.MolToSmiles(mol),
            "qed": qed,
            "logp": logp,
            "status": status,
            "image": f"data:image/png;base64,{img_str}"
        }
    except Exception as e:
        print(f"Prop Error: {e}")
        return default_props