from fastapi import APIRouter, HTTPException, Depends, BackgroundTasks
from bson import ObjectId
from datetime import datetime
from backend.database import db
from backend.models import JobCreate, JobDB, JobResponse, UserResponse
from backend.routers.auth import get_current_user
import asyncio

router = APIRouter(prefix="/jobs", tags=["Jobs"])

async def process_job(job_id: str):
    """
    Background Task: Executes ML Logic (Generation/Optimization).
    Triggered by POST /jobs via BackgroundTasks.
    """
    print(f"🔄 [Worker] Starting Job {job_id}...")
    
    # --- Local Imports to prevent Circular Dependency & Access Globals ---
    # We need access to the initialized model, predictor, and config in main.py
    from backend.main import (
        model, predictor, DEVICE, MODE, 
        decode_tensor, smiles_to_latent, is_valid_candidate, 
        get_mol_from_sequence, calculate_properties, optimize_latent_vector
    )
    import torch

    try:
        # 1. Update Status -> PROCESSING
        await db.jobs.update_one(
            {"_id": ObjectId(job_id)},
            {"$set": {"status": "PROCESSING", "updated_at": datetime.utcnow()}}
        )
        
        # 2. Fetch Job Params
        job = await db.jobs.find_one({"_id": ObjectId(job_id)})
        if not job:
            print(f"❌ [Worker] Job {job_id} not found in DB execution.")
            return

        task_type = job.get("task_type")
        params = job.get("params", {})
        results = []
        
        # Checking Resource Availability
        if model is None:
            raise Exception("Model not loaded on server.")

        # =========================================================================
        # TASK: GENERATE RANDOM
        # =========================================================================
        if task_type == "GENERATE_RANDOM":
            num_molecules = params.get("num_molecules", 5)
            attempts = 0
            # Reuse logic from main.py /generate
            while len(results) < num_molecules and attempts < 100:
                attempts += 1
                with torch.no_grad():
                    # Generate 5x more samples than needed to buffer for invalid ones
                    indices = model.sample(num_molecules * 5, DEVICE, temperature=2.0)
                    for i in range(len(indices)):
                        seq = decode_tensor(indices[i])
                        mol = get_mol_from_sequence(seq, mode=MODE)
                        
                        if is_valid_candidate(mol):
                            # Dedup in current batch
                            if any(r['sequence'] == seq for r in results): 
                                continue
                            
                            props = calculate_properties(mol)
                            if props['valid'] and props['image']:
                                results.append({"sequence": seq, "properties": props})
                                if len(results) >= num_molecules: break
        
        # =========================================================================
        # TASK: GENERATE TARGETED
        # =========================================================================
        elif task_type == "GENERATE_TARGETED":
            if predictor is None: raise Exception("Predictor not loaded.")
            
            target_qed = params.get("target_qed", 0.8)
            target_logp = params.get("target_logp", 2.5)
            target_sas = params.get("target_sas", 3.0)
            num_molecules = params.get("num_molecules", 5)
            
            target_props = [target_qed, target_logp, target_sas]
            INTERNAL_BATCH = 300
            
            # 1. Optimization Step
            with torch.enable_grad():
                z = torch.randn(INTERNAL_BATCH, 128).to(DEVICE)
                z_opt = optimize_latent_vector(z, predictor, target_props)
            
            # 2. Decoding Step
            with torch.no_grad():
                indices = model.decode(z_opt, DEVICE, temperature=0.8)
                
                temp_results = []
                for i in range(len(indices)):
                    seq = decode_tensor(indices[i])
                    mol = get_mol_from_sequence(seq, mode=MODE)
                    
                    if is_valid_candidate(mol):
                        props = calculate_properties(mol)
                        if props['valid'] and props['image']:
                            score = props.get('score', 0)
                            props['status'] = f"🎯 Targeted (Score {score})"
                            temp_results.append({"sequence": seq, "properties": props, "score": score})
            
            # Sort & Trim
            temp_results.sort(key=lambda x: x['score'], reverse=True)
            results = temp_results[:num_molecules]

        # =========================================================================
        # TASK: OPTIMIZE LEAD
        # =========================================================================
        elif task_type == "OPTIMIZE_LEAD":
            smiles = params.get("smiles")
            if not smiles: raise Exception("Missing 'smiles' parameter.")
            
            z_lead = smiles_to_latent(smiles)
            if z_lead is None: raise Exception("Invalid input molecule (could not encode).")
            
            BATCH_SIZE = 200
            with torch.no_grad():
                z_batch = z_lead.repeat(BATCH_SIZE, 1)
                
                # Create Cloud
                noise_close = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.1
                noise_far = torch.randn(BATCH_SIZE // 2, 128).to(DEVICE) * 0.3
                
                z_batch[:BATCH_SIZE//2] += noise_close
                z_batch[BATCH_SIZE//2:] += noise_far
                
                # Decode
                indices = model.decode(z_batch, DEVICE, temperature=0.9)
                
                temp_results = []
                unique_smiles = set()
                unique_smiles.add(smiles) # Exclude self
                
                for i in range(BATCH_SIZE):
                    seq = decode_tensor(indices[i])
                    mol = get_mol_from_sequence(seq, mode=MODE)
                    
                    if is_valid_candidate(mol):
                        props = calculate_properties(mol)
                        if props['valid'] and props['image']:
                            score = props.get('score', 0)
                            s_out = props['smiles']
                            
                            if s_out not in unique_smiles:
                                unique_smiles.add(s_out)
                                props['status'] = f"✨ Optimized ({score})"
                                temp_results.append({"sequence": seq, "properties": props, "score": score})
            
            # Sort
            temp_results.sort(key=lambda x: x['score'], reverse=True)
            results = temp_results[:5]

        else:
            raise Exception(f"Unknown Task Type: {task_type}")

        # ---------------------------------------------------
        # 3. Update Status -> COMPLETED
        # ---------------------------------------------------
        await db.jobs.update_one(
            {"_id": ObjectId(job_id)},
            {"$set": {
                "status": "COMPLETED", 
                "result": results,
                "updated_at": datetime.utcnow()
            }}
        )
        print(f"✅ [Worker] Job {job_id} Completed. Found {len(results)} items.")

    except Exception as e:
        print(f"💥 [Worker] Job {job_id} Failed: {e}")
        import traceback
        traceback.print_exc()
        await db.jobs.update_one(
            {"_id": ObjectId(job_id)},
            {"$set": {
                "status": "FAILED", 
                "error": str(e),
                "updated_at": datetime.utcnow()
            }}
        )

@router.post("/", response_model=JobResponse, status_code=202)
async def create_job(
    job_req: JobCreate, 
    background_tasks: BackgroundTasks,
    current_user: UserResponse = Depends(get_current_user)
):
    """
    Async Entry Point for Long-Running Tasks.
    """
    # 1. Create Job Document
    job_db = JobDB(
        user_id=current_user.id,
        task_type=job_req.task_type,
        params=job_req.params,
        status="PENDING"
    )
    
    new_job = await db.jobs.insert_one(job_db.dict(by_alias=True))
    created_id = str(new_job.inserted_id)
    
    # 2. Schedule Background Task
    background_tasks.add_task(process_job, created_id)
    
    # 3. Return Immediate Response
    return JobResponse(
        job_id=created_id,
        status="PENDING",
        created_at=job_db.created_at,
        updated_at=job_db.updated_at
    )

@router.get("/{job_id}", response_model=JobResponse)
async def get_job_status(job_id: str, current_user: UserResponse = Depends(get_current_user)):
    if not ObjectId.is_valid(job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID")
        
    job = await db.jobs.find_one({"_id": ObjectId(job_id), "user_id": current_user.id})
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
        
    return JobResponse(
        job_id=str(job["_id"]),
        status=job["status"],
        created_at=job["created_at"],
        updated_at=job["updated_at"],
        result=job.get("result"),
        error=job.get("error")
    )
