import torch
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import DataLoader
from rdkit import Chem
from rdkit.Chem import rdchem


def atom_features(atom):
    # 1. Atomic number (one-hot over [6,7,8,9,16,17]) -> 6
    atomic_nums = [6, 7, 8, 9, 16, 17]
    atomic_one_hot = [int(atom.GetAtomicNum() == a) for a in atomic_nums]

    # 2. Degree (one-hot 0–4) -> 5
    degree_one_hot = [int(atom.GetDegree() == d) for d in range(5)]

    # 3. Hybridization (SP, SP2, SP3) -> 3
    hyb_map = {
        rdchem.HybridizationType.SP: 0,
        rdchem.HybridizationType.SP2: 1,
        rdchem.HybridizationType.SP3: 2,
    }
    hyb = atom.GetHybridization()
    hyb_one_hot = [0, 0, 0]
    if hyb in hyb_map:
        hyb_one_hot[hyb_map[hyb]] = 1

    # 4. Aromaticity -> 1
    aromaticity = [float(atom.GetIsAromatic())]

    # 5. Formal charge (scaled) -> 1
    formal_charge = [atom.GetFormalCharge() / 5.0]

    # 6. Number of hydrogens (scaled) -> 1
    num_hs = [atom.GetTotalNumHs() / 4.0]

    # Total: 6 + 5 + 3 + 1 + 1 + 1 = 17 features
    return atomic_one_hot + degree_one_hot + hyb_one_hot + aromaticity + formal_charge + num_hs


def bond_features(bond):
    # 1. Bond type (one-hot) -> 4
    bond_type_map = {
        rdchem.BondType.SINGLE: 0,
        rdchem.BondType.DOUBLE: 1,
        rdchem.BondType.TRIPLE: 2,
        rdchem.BondType.AROMATIC: 3,
    }
    bt = bond.GetBondType()
    bond_type_one_hot = [0, 0, 0, 0]
    if bt in bond_type_map:
        bond_type_one_hot[bond_type_map[bt]] = 1

    # 2. Conjugation -> 1
    conjugated = [float(bond.GetIsConjugated())]

    # 3. Ring membership -> 1
    in_ring = [float(bond.IsInRing())]

    # Total: 4 + 1 + 1 = 6 features
    return bond_type_one_hot + conjugated + in_ring


def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    node_features = []
    for atom in mol.GetAtoms():
        node_features.append(atom_features(atom))

    x = torch.tensor(node_features, dtype=torch.float)

    # Edge dimension
    edge_dim = 6

    bonds = mol.GetBonds()
    if mol.GetNumBonds() == 0:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, edge_dim), dtype=torch.float)
        return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

    edge_indices = []
    edge_attrs = []

    for bond in bonds:
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        feat = bond_features(bond)

        # i -> j
        edge_indices.append([i, j])
        edge_attrs.append(feat)

        # j -> i
        edge_indices.append([j, i])
        edge_attrs.append(feat)

    edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attrs, dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


class MoleculeDataset(Dataset):
    def __init__(self, smiles_list):
        super(MoleculeDataset, self).__init__()
        # Filtering invalid SMILES during initialization
        self.data_list = [smiles_to_graph(s) for s in smiles_list]
        self.data_list = [d for d in self.data_list if d is not None]

    def len(self):
        return len(self.data_list)

    def get(self, idx):
        return self.data_list[idx]
