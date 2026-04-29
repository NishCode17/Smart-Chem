from pydantic import BaseModel, ConfigDict
from typing import List, Dict, Optional, Any
from datetime import datetime

class MolProperties(BaseModel):
    logp: float
    qed: float
    mw: float
    tpsa: float
    hbd: int
    hba: int
    rot_bonds: int

class AdmetProperties(BaseModel):
    absorption: float
    distribution: float
    metabolism: float
    excretion: float
    toxicity: float

class MoleculeResponse(BaseModel):
    id: str
    project_id: str
    name: str
    smiles: str
    generated_by: str
    tags: List[str]
    properties: MolProperties
    admet: Optional[AdmetProperties]
    admet_full: Optional[Dict[str, Any]] = {}
    tox_alerts: Optional[Dict[str, Any]]
    created_at: datetime

m = {
    "_id": "60a7b4f5a3b9f3b3b4f5a3b9",
    "user_id": "test",
    "project_id": "test",
    "name": "test",
    "smiles": "C",
    "generated_by": "test",
    "tags": [],
    "properties": {
        "logp": 1,
        "qed": 1,
        "mw": 1,
        "tpsa": 1,
        "hbd": 1,
        "hba": 1,
        "rot_bonds": 1
    },
    "admet": None,
    "tox_alerts": None,
    "created_at": datetime.now()
}

try:
    resp = MoleculeResponse(**m, id=str(m["_id"]))
    print(resp.model_dump_json())
except Exception as e:
    print("Error:", e)
