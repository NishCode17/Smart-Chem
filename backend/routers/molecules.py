
from fastapi import APIRouter, HTTPException, Depends, status
from typing import List, Optional
from bson import ObjectId
from backend.database import db
from backend.models import MoleculeCreate, MoleculeResponse, MoleculeDB, UserResponse
from backend.routers.auth import get_current_user
from datetime import datetime

router = APIRouter(prefix="/molecules", tags=["Molecules"])

@router.post("/", response_model=MoleculeResponse)
async def create_molecule(molecule: MoleculeCreate, current_user: UserResponse = Depends(get_current_user)):
    # Verify project belongs to user
    print(">>> Incoming molecule payload:", molecule)
    print(">>> Incoming project_id:", molecule.project_id)
    print(">>> SAVE MOL PAYLOAD:", molecule)

    try:
        project_oid = ObjectId(molecule.project_id)
    except:
        raise HTTPException(status_code=400, detail="Invalid project ID format")

    project = await db.projects.find_one({"_id": project_oid})
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    molecule_db = MoleculeDB(
        user_id=current_user.id,
        project_id=molecule.project_id,
        name=molecule.name,
        smiles=molecule.smiles,
        generated_by=molecule.generated_by,
        tags=molecule.tags,
        properties=molecule.properties,
        admet=molecule.admet if molecule.admet else {},
        tox_alerts=molecule.tox_alerts if molecule.tox_alerts else {}
    )
    
    result = await db.molecules.insert_one(molecule_db.dict(by_alias=True))
    created_mol = await db.molecules.find_one({"_id": result.inserted_id})
    return MoleculeResponse(**created_mol, id=str(created_mol["_id"]))

@router.get("/", response_model=List[MoleculeResponse])
async def get_molecules(project_id: Optional[str] = None, current_user: UserResponse = Depends(get_current_user)):
    query = {"user_id": current_user.id}
    if project_id:
        if not ObjectId.is_valid(project_id):
            raise HTTPException(status_code=400, detail="Invalid project ID")
        query["project_id"] = project_id
        
    molecules = await db.molecules.find(query).sort("created_at", -1).to_list(100)
    return [MoleculeResponse(**m, id=str(m["_id"])) for m in molecules]

@router.get("/{molecule_id}", response_model=MoleculeResponse)
async def get_molecule(molecule_id: str, current_user: UserResponse = Depends(get_current_user)):
    if not ObjectId.is_valid(molecule_id):
         raise HTTPException(status_code=400, detail="Invalid molecule ID")
    
    molecule = await db.molecules.find_one({"_id": ObjectId(molecule_id), "user_id": current_user.id})
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
        
    return MoleculeResponse(**molecule, id=str(molecule["_id"]))

@router.delete("/{molecule_id}", status_code=204)
async def delete_molecule(molecule_id: str, current_user: UserResponse = Depends(get_current_user)):
    if not ObjectId.is_valid(molecule_id):
        raise HTTPException(status_code=400, detail="Invalid molecule ID")
        
    result = await db.molecules.delete_one({"_id": ObjectId(molecule_id), "user_id": current_user.id})
    if result.deleted_count == 0:
        raise HTTPException(status_code=404, detail="Molecule not found")
