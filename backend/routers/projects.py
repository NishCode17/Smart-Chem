from fastapi import APIRouter, HTTPException, Depends
from typing import List
from bson import ObjectId
from backend.database import db
from backend.models import ProjectCreate, ProjectResponse, ProjectDB, UserResponse
from backend.routers.auth import get_current_user

router = APIRouter(prefix="/projects", tags=["Projects"])

@router.post("/", response_model=ProjectResponse)
async def create_project(project: ProjectCreate, current_user: UserResponse = Depends(get_current_user)):
    project_db = ProjectDB(
        user_id=current_user.id,
        name=project.name,
        description=project.description
    )
    print(">>> DEBUG PROJECT INSERT:", project_db.model_dump(by_alias=True))

    result = await db.projects.insert_one(project_db.dict(by_alias=True))
    created = await db.projects.find_one({"_id": result.inserted_id})

    return ProjectResponse(
        id=str(created["_id"]),
        name=created["name"],
        description=created.get("description"),
        created_at=created["created_at"]
    )


@router.get("/", response_model=List[ProjectResponse])
async def get_projects(current_user: UserResponse = Depends(get_current_user)):
    projects = await db.projects.find({"user_id": current_user.id}).to_list(200)

    return [
        ProjectResponse(
            id=str(p["_id"]),
            name=p["name"],
            description=p.get("description"),
            created_at=p["created_at"],
        )
        for p in projects
    ]


@router.delete("/{project_id}", status_code=204)
async def delete_project(project_id: str, current_user: UserResponse = Depends(get_current_user)):

    if not ObjectId.is_valid(project_id):
        raise HTTPException(status_code=400, detail="Invalid project ID")

    result = await db.projects.delete_one({
        "_id": ObjectId(project_id),
        "user_id": current_user.id
    })

    if result.deleted_count == 0:
        raise HTTPException(status_code=404, detail="Project not found")
