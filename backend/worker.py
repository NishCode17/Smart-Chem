import asyncio
import os
import sys
import traceback
from datetime import datetime
from bson import ObjectId

# Add project root
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.database import db
import backend.ml_executor as ml_exec

# -----------------------------------------------------------------------------
# WORKER LOOP
# -----------------------------------------------------------------------------
async def worker_loop():
    # Load resources via shared executor
    ml_exec.load_resources()
    print("🚀 [Worker] Looking for jobs...")
    
    while True:
        try:
            # Atomic Claim
            job = await db.jobs.find_one_and_update(
                {"status": "PENDING"},
                {"$set": {"status": "PROCESSING", "updated_at": datetime.utcnow()}},
                sort=[("created_at", 1)],
                return_document=True
            )
            
            if not job:
                await asyncio.sleep(2)
                continue
                
            job_id = job["_id"]
            print(f"🔄 Processing Job {job_id} ({job['task_type']})")
            
            try:
                # ROUTE TASK
                result = None
                params = job['params']
                
                if job['task_type'] == "OPTIMIZE_LEAD":
                    # DELEGATE TO EXECUTOR
                    result = ml_exec.run_lead_optimization(smiles=params['smiles'])
                    
                else:
                    raise ValueError(f"Unknown Task: {job['task_type']}")
                    
                # SUCCESS
                await db.jobs.update_one(
                    {"_id": job_id},
                    {"$set": {
                        "status": "COMPLETED",
                        "result": result,
                        "updated_at": datetime.utcnow()
                    }}
                )
                print(f"✅ Job {job_id} Completed")
                
            except Exception as e:
                print(f"💥 Job {job_id} Failed: {e}")
                traceback.print_exc()
                await db.jobs.update_one(
                    {"_id": job_id},
                    {"$set": {
                        "status": "FAILED",
                        "error": str(e),
                        "updated_at": datetime.utcnow()
                    }}
                )
                
        except Exception as e:
            print(f"Worker Error: {e}")
            await asyncio.sleep(5)

if __name__ == "__main__":
    asyncio.run(worker_loop())
