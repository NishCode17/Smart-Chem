
import asyncio
import os
from motor.motor_asyncio import AsyncIOMotorClient
from dotenv import load_dotenv

load_dotenv()

async def reset_user():
    mongo_url = os.getenv("MONGO_URI")
    db_name = os.getenv("MONGO_DB")
    
    print(f"Connecting to MongoDB...")
    client = AsyncIOMotorClient(mongo_url)
    db = client[db_name]
    
    username = "nish"
    result = await db.users.delete_one({"username": username})
    
    if result.deleted_count > 0:
        print(f"Successfully deleted user: {username}")
    else:
        print(f"User {username} not found.")

if __name__ == "__main__":
    asyncio.run(reset_user())
