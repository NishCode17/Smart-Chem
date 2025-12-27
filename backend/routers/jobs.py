from fastapi import APIRouter, HTTPException, Depends
from bson import ObjectId
from datetime import datetime
from typing import Dict, Any, Optional
from pydantic import BaseModel

from backend.database import db
from backend.models import JobDB, JobResponse, UserResponse
from backend.routers.auth import get_current_user

router = APIRouter(prefix="/jobs", tags=["Async Jobs"])

# ---------------------------------------------------------
# Request Models (Specific for /optimize)
# ---------------------------------------------------------
class OptimizeJobRequest(BaseModel):
    smiles: str
    target_qed: float = 0.8
    target_logp: float = 2.5
    
# ---------------------------------------------------------
# 1. POST /jobs/optimize
# ---------------------------------------------------------
@router.post("/optimize", response_model=JobResponse, status_code=202)
async def create_optimization_job(
    req: OptimizeJobRequest, 
    current_user: UserResponse = Depends(get_current_user)
):
    """
    Creates a new optimization job.
    Returns the job_id immediately. Status is PENDING.
    """
    # 1. Define Protocol Payload
    payload = {
        "smiles": req.smiles,
        "target_qed": req.target_qed,
        "target_logp": req.target_logp
    }

    # 2. Create Job Document (PENDING)
    job = JobDB(
        user_id=current_user.id,
        task_type="OPTIMIZE_LEAD",
        params=payload,
        status="PENDING", # Waits for external worker
        created_at=datetime.utcnow(),
        updated_at=datetime.utcnow()
    )
    
    # 3. Save to DB
    new_job = await db.jobs.insert_one(job.dict(by_alias=True))
    
    # 4. Return ID
    return JobResponse(
        job_id=str(new_job.inserted_id),
        status=job.status,
        created_at=job.created_at,
        updated_at=job.updated_at
    )

# ---------------------------------------------------------
# 2. GET /jobs/{job_id}
# ---------------------------------------------------------
@router.get("/{job_id}", response_model=JobResponse)
async def get_job_status(
    job_id: str, 
    current_user: UserResponse = Depends(get_current_user)
):
    """
    Polls the status of a job.
    """
    if not ObjectId.is_valid(job_id):
        raise HTTPException(400, "Invalid Job ID")
        
    job = await db.jobs.find_one({
        "_id": ObjectId(job_id),
        "user_id": current_user.id
    })
    
    if not job:
        raise HTTPException(404, "Job not found")
        
    return JobResponse(
        job_id=str(job["_id"]),
        status=job["status"],
        created_at=job["created_at"],
        updated_at=job["updated_at"],
        result=job.get("result"),
        error=job.get("error")
    )
